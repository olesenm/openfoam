/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "masterOFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::masterOFstream::checkWrite
(
    const fileName& fName,
    const char* str,
    std::streamsize len
)
{
    if (!str || !(len > 0))
    {
        // Can skip everything if there is nothing to write
        return;
    }

    Foam::mkDir(fName.path());

    OFstream os
    (
        atomic_,
        fName,
        IOstreamOption(IOstreamOption::BINARY, version(), compression_),
        append_
    );
    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Could not open file " << fName << nl
            << exit(FatalIOError);
    }

    // Write characters directly to std::ostream
    os.writeRaw(str, len);

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Failed writing to " << fName << nl
            << exit(FatalIOError);
    }
}


void Foam::masterOFstream::commit()
{
    // Take ownership of serialized content
    DynamicList<char> charData(OCharStream::release());

    if (!UPstream::parRun())
    {
        // Write (non-empty) data
        checkWrite(pathName_, charData);
    }
    else
    {
        // Ignore content if not writing
        if (!writeOnProc_)
        {
            charData.clear();
        }

        List<fileName> filePaths(UPstream::nProcs(comm_));
        filePaths[UPstream::myProcNo(comm_)] = pathName_;
        Pstream::gatherList(filePaths, UPstream::msgType(), comm_);

        // Test for identical output paths
        bool uniform =
        (
            UPstream::master(comm_)
         && fileOperation::uniformFile(filePaths)
        );

        Pstream::broadcast(uniform, comm_);

        if (uniform)
        {
            // Identical file paths - write on master
            if (UPstream::master(comm_) && writeOnProc_)
            {
                checkWrite(pathName_, charData);
            }
            return;
        }

        // Different files
        // ---------------
        //
        // Non-sparse (most ranks have writeOnProc_ == true),
        // so gather sizes first and use PEX-like handling,
        // with polling for when data becomes available.
        //
        // Could also consider double buffering + write to reduce
        // memory overhead.

        // Or int64_t
        const label dataSize =
        (
            (UPstream::is_subrank(comm_) && writeOnProc_)
          ? charData.size()
          : 0
        );

        const labelList recvSizes
        (
            UPstream::listGatherValues<label>(dataSize, comm_)
        );

        // Receive from these procs
        DynamicList<int> recvProcs;

        if (UPstream::master(comm_))
        {
            // Sorted by message size
            labelList order(Foam::sortedOrder(recvSizes));
            recvProcs.reserve_exact(order.size());

            // Want to receive large messages first. Ignore empty slots
            forAllReverse(order, i)
            {
                const label proci = order[i];

                // Ignore empty slots
                if (recvSizes[proci] > 0)
                {
                    recvProcs.push_back(proci);
                }
            }
        }

        // Non-blocking communication
        const auto startOfRequests = UPstream::nRequests();

        // Some unique tag for this read/write grouping (extra precaution)
        const int messageTag = (UPstream::msgType() + 256);

        if (UPstream::is_subrank(comm_) && dataSize > 0)
        {
            // Send to content to master
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::masterNo(),
                charData.cdata_bytes(),
                charData.size_bytes(),
                messageTag,
                comm_
            );
        }
        else if (UPstream::master(comm_))
        {
            // The receive slots
            List<List<char>> recvBuffers(UPstream::nProcs(comm_));

            // Receive from these procs (non-empty slots)
            for (const int proci : recvProcs)
            {
                auto& slot = recvBuffers[proci];
                slot.resize_nocopy(recvSizes[proci]);

                // Receive content
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    slot.data_bytes(),
                    slot.size_bytes(),
                    messageTag,
                    comm_
                );
            }

            if (writeOnProc_)
            {
                // Write non-empty master data
                checkWrite(pathName_, charData);
                charData.clear();
            }

            // Poll for completed receive requests and dispatch
            DynamicList<int> indices(recvProcs.size());
            while
            (
                UPstream::waitSomeRequests
                (
                    startOfRequests,
                    recvProcs.size(),
                   &indices
                )
            )
            {
                for (const int i : indices)
                {
                    const int proci = recvProcs[i];
                    auto& slot = recvBuffers[proci];

                    // Write non-empty sub-proc data
                    checkWrite(filePaths[proci], slot);

                    // Eager cleanup
                    slot.clear();
                }
            }
        }

        UPstream::waitRequests(startOfRequests);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterOFstream::masterOFstream
(
    IOstreamOption::atomicType atomic,
    const int communicator,
    const fileName& pathName,
    IOstreamOption streamOpt,
    IOstreamOption::appendType append,
    const bool writeOnProc
)
:
    OCharStream(streamOpt),
    pathName_(pathName),
    atomic_(atomic),
    compression_(streamOpt.compression()),
    append_(append),
    writeOnProc_(writeOnProc),
    comm_(communicator < 0 ? UPstream::worldComm : communicator)
{
    // Start with a slightly larger buffer
    OCharStream::reserve(4*1024);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterOFstream::~masterOFstream()
{
    commit();
}


// ************************************************************************* //
