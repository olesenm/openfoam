/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Mark Olesen
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

Application
    Test-parallel-file-write1

Description
    Simple test of writing with MPI/IO

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Switch.H"
#include "UPstreamFile.H"
#include "SpanStream.H"

using namespace Foam;

template<class IntType>
void zeropadded(std::ostream& os, IntType val, char fillch = '0')
{
    // set fill char and width
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    fillch = os.fill(fillch);
    os.width(std::numeric_limits<IntType>::digits10+1);
    os  << val;
    // restore fill char
    os.fill(fillch);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();

    #include "setRootCase.H"

    if (!UPstream::parRun())
    {
        Info<< "###############" << nl
            << "Not running in parallel. Stopping now" << nl
            << "###############" << endl;
        return 1;
    }

    Info<< "Create time (without controlDict)\n" << endl;

    auto runTimePtr = Time::New();
    auto& runTime = runTimePtr();

    const auto comm_ = UPstream::worldComm;

    const auto myProc = UPstream::myProcNo(comm_);
    const auto nProcs = UPstream::nProcs(comm_);

    // Some content
    OCharStream charset;
    for (int i = 0; i < 10; ++i)
    {
        charset<< char('A' + i);
    }

    // Content buffer
    OCharStream os(IOstreamOption::BINARY);

    {
        const auto v = charset.view();

        os << nl;
        os.beginBlock(word("rank" + Foam::name(myProc)));

        for (int repeat = 0;  repeat <= myProc; ++repeat)
        {
            os << indent << word("entry" + Foam::name(repeat))
                << ' ' << word("List<char>");
            // os << nl;
            os << ' ';
            os << label(v.size());
            os.write(v.data(), v.size());
            // os << nl;
            os.endEntry();
        }

        os.endBlock();
    }


    // All content now exists - commit to disk
    const std::string_view blockData(os.view());
    const int64_t blockSize(blockData.size());

    // The sizes (without footer!) gathered onto the master
    List<int64_t> sizes
    (
        UPstream::listGatherValues<int64_t>(blockSize, comm_)
    );

    // Overall header - most flexible to keep separate from block content
    OCharStream header;
    if (UPstream::master(comm_))
    {
        header
            << "Simple MPI/IO test with " << nProcs << " ranks" << nl << nl;

        ocharstream labelbuf;
        labelbuf.reserve_exact(32);

        // Position before writing a label
        auto labelBegin = header.tellp();

        header.beginBlock("meta");
        {
            header << indent << word("data.start") << ' ';

            labelBegin = header.tellp();

            // Add the start value (placeholder)
            {
                labelbuf.rewind();
                zeropadded(labelbuf, label(0));
                header.append(labelbuf.view());
            }

            header.endEntry();

            header << indent << word("data.sizes") << nl;
            sizes.writeList(header);  // flatOutput
            header.endEntry();
        }
        header.endBlock();

        header << nl;
        IOobject::writeDivider(header);

        // Now update with the correct size
        {
            labelbuf.rewind();
            zeropadded(labelbuf, label(header.view().size()));
            header.overwrite(labelBegin, labelbuf.view());
        }
    }

    // Overall footer.
    // Will be written by the last block, but format for everyone
    // so that the size is known
    OCharStream footer;
    {
        IOobject::writeEndDivider(footer);
    }


    // Calculate the begin offsets for each block and total file size

    int64_t totalSize(header.view().size());
    for (auto& val : sizes)
    {
        const auto count = val;
        val = totalSize;
        totalSize += count;
    }
    totalSize += int64_t(footer.view().size());


    // The file begin offset for my block
    const int64_t blockOffset = UPstream::listScatterValues(sizes, comm_);

    // Everyone needs to know this as well
    Pstream::broadcast(totalSize, comm_);

    Pout<< "write size=" << blockSize
        << " at=" << blockOffset << " total=" << totalSize << nl;

    // The last block also gets the footer data to write
    if (myProc == nProcs-1)
    {
        os.extend_exact(footer.view().size());
        os.append(footer.view());
    }

    {
        UPstream::File file;

        bool ok = file.open_write
        (
            comm_,
            runTime.globalPath()/"mpiio-test1.txt",
            IOstreamOption::ATOMIC
        );

        if (ok)
        {
            Info<< "writing: " << file.name() << nl;

            // header from master
            if (UPstream::master(comm_))
            {
                // A no-op for empty buffer
                ok = file.write_at(0, header.view());
            }

            // data from all - footer is already in the last block
            ok = file.write_at_all(blockOffset, os.view());
        }

        file.set_size(totalSize);
        file.close();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
