/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    Test-globalIndex3

Description
    Tests for globalIndex with node-wise splitting

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"
#include "globalMeshData.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IndirectList.H"
#include "IOstreams.H"
#include "Random.H"
#include "openfoam_mpi.H"

// pre-scan for "-split-size NUM"
int option_splitsize(int argc, char *argv[])
{
    int ivalue = -1;
    for (int argi = 1; argi < argc-1; ++argi)
    {
        if (strcmp(argv[argi], "-split-size") == 0)
        {
            ++argi;
            ivalue = atoi(argv[argi]);
        }
    }

    return ivalue;
}

using namespace Foam;

template<class T>
void printList(Ostream& os, const UList<T>& list)
{
    os << list.size() << " " << flatOutput(list) << nl;
}

void printGlobalIndex(Ostream& os, const globalIndex& gi)
{
    printList(os, gi.offsets());
}


template<class ProcIDsContainer, class Type>
void globalIndexGather
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& fld,
    UList<Type>& allFld,    // must be adequately sized on master
    const int tag,
    UPstream::commsTypes commsType,
    bool useWindow = false
)
{
    // low-level: no parRun guard
    const int masterProci = procIDs.size() ? procIDs[0] : 0;

    // Protection for disjoint calls
    if (FOAM_UNLIKELY(!UPstream::is_rank(comm)))
    {
        FatalErrorInFunction
            << "Calling with process not on the communicator"
            << Foam::abort(FatalError);
    }

    // Require contiguous data for non-blocking
    if constexpr (!is_contiguous_v<Type>)
    {
        if (commsType == UPstream::commsTypes::nonBlocking)
        {
            commsType = UPstream::commsTypes::scheduled;
        }
    }

    const auto startOfRequests = UPstream::nRequests();


    // Very hard-coded at the moment
    int returnCode = MPI_SUCCESS;
    const int nCmpts = pTraits<Type>::nComponents;

    MPI_Win win;
    MPI_Datatype dataType = MPI_DOUBLE;
    if (useWindow)
    {
        using cmptType = typename pTraits<Type>::cmptType;

        if (std::is_same_v<float, cmptType>)
        {
            dataType = MPI_FLOAT;
        }
        else if (std::is_same_v<double, cmptType>)
        {
            dataType = MPI_DOUBLE;
        }
        else
        {
            // Not supported
            useWindow = false;
        }
    }

    if (useWindow)
    {
        MPI_Comm mpiComm =
            PstreamUtils::Cast::to_mpi(UPstream::Communicator::lookup(comm));

        char commName[MPI_MAX_OBJECT_NAME];
        int nameLen = 0;

        if
        (
            MPI_COMM_NULL != mpiComm
         && MPI_SUCCESS == MPI_Comm_get_name(mpiComm, commName, &nameLen)
         && (nameLen > 0)
        )
        {
            Pout<< "window on " << commName << nl;
        }

        if (UPstream::myProcNo(comm) == masterProci || fld.empty())
        {
            // Collective
            returnCode = MPI_Win_create
            (
                nullptr,
                0,
                1,  // disp_units
                MPI_INFO_NULL,
                mpiComm,
               &win
            );
        }
        else
        {
            // Collective
            returnCode = MPI_Win_create
            (
                const_cast<char *>(fld.cdata_bytes()),
                fld.size_bytes(),
                sizeof(Type),  // disp_units
                MPI_INFO_NULL,
                mpiComm,
               &win
            );
        }

        if (MPI_SUCCESS != returnCode || MPI_WIN_NULL == win)
        {
            FatalErrorInFunction
                << "MPI_Win_create() failed"
                << Foam::abort(FatalError);
            // return nullptr;
        }
    }


    if (UPstream::myProcNo(comm) == masterProci)
    {
        const label total = off.back();  // == totalSize()

        if (allFld.size() < total)
        {
            FatalErrorInFunction
                << "[out] UList size=" << allFld.size()
                << " too small to receive " << total << nl
                << Foam::abort(FatalError);
        }


        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied.
        // Also handle the case where we have a slice of the full
        // list.
        {
            SubList<Type> dst(allFld, off[1]-off[0], off[0]);
            SubList<Type> src(fld, off[1]-off[0]);

            if (!dst.empty() && (dst.data() != src.data()))
            {
                dst = src;
            }
        }

        if (useWindow)
        {
            MPI_Win_lock_all(MPI_MODE_NOCHECK, win);
        }

        for (label i = 1; i < procIDs.size(); ++i)
        {
            SubList<Type> slot(allFld, off[i+1]-off[i], off[i]);

            if (slot.empty())
            {
                // Nothing to do
            }
            else if (useWindow)
            {
                returnCode = MPI_Get
                (
                    // origin
                    slot.data(),
                    slot.size()*(nCmpts),
                    dataType,

                    // target
                    procIDs[i],
                    0,  // displacement
                    slot.size()*(nCmpts),
                    dataType,
                    win
                );

                if (MPI_SUCCESS != returnCode)
                {
                    FatalErrorInFunction
                        << "MPI_Get failed"
                        << Foam::abort(FatalError);
                    // return nullptr;
                }
            }
            else if constexpr (is_contiguous_v<Type>)
            {
                UIPstream::read
                (
                    commsType,
                    procIDs[i],
                    slot,
                    tag,
                    comm
                );
            }
            else
            {
                IPstream::recv(slot, procIDs[i], tag, comm);
            }
        }

        if (useWindow)
        {
            MPI_Win_unlock_all(win);
        }
    }
    else if (!useWindow)
    {
        if (fld.empty())
        {
            // Nothing to do
        }
        else if constexpr (is_contiguous_v<Type>)
        {
            UOPstream::write
            (
                commsType,
                masterProci,
                fld,
                tag,
                comm
            );
        }
        else
        {
            OPstream::send(fld, commsType, masterProci, tag, comm);
        }
    }

    if (useWindow)
    {
        // Collective
        MPI_Win_free(&win);
    }

    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        // Wait for outstanding requests
        UPstream::waitRequests(startOfRequests);
    }
}


// Report inter-node/intra-node offsets
static void reportOffsets(const globalIndex& gi)
{
    labelList interNodeOffsets;
    labelList localNodeOffsets;
    labelRange nodeRange;

    const label numProc = UPstream::nProcs(UPstream::commConstWorld());

    gi.splitNodeOffsets
    (
        interNodeOffsets,
        localNodeOffsets,
        UPstream::worldComm
    );

    const auto interNodeComm = UPstream::commInterNode();

    // Only communicate to the node leaders
    labelList allOffsets;
    if (UPstream::is_rank(interNodeComm))
    {
        // Send top-level offsets to the node leaders
        if (UPstream::master(interNodeComm))
        {
            allOffsets = gi.offsets();
        }
        else  // ie, UPstream::is_subrank(interNodeComm)
        {
            allOffsets.resize_nocopy(numProc+1);
        }

        UPstream::broadcast
        (
            allOffsets.data(),
            allOffsets.size(),
            interNodeComm
        );
    }

    // Ranges (node leaders only)
    if (UPstream::is_rank(interNodeComm))
    {
        const auto& procIds = UPstream::procID(interNodeComm);
        const int ranki = UPstream::myProcNo(interNodeComm);

        // For reporting
        nodeRange.reset
        (
            procIds[ranki],
            (
                (ranki+1 < procIds.size() ? procIds[ranki+1] : numProc)
              - procIds[ranki]
            )
        );
    }

    Pout<< "node-range: " << nodeRange << nl;
    Pout<< "all-offset: "; printList(Pout, allOffsets);
    Pout<< "inter-offset: "; printList(Pout, interNodeOffsets);
    Pout<< "intra-offset: "; printList(Pout, localNodeOffsets);
}


template<class Type>
void globalIndexGather
(
    const globalIndex& gi,
    const UList<Type>& sendData,
    List<Type>& allData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm = UPstream::worldComm,
    bool useWindow = false
)
{
    if (!UPstream::parRun())
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }

    if (UPstream::master(comm))
    {
        allData.resize_nocopy(gi.offsets().back());  // == totalSize()
    }
    else
    {
        allData.clear();  // zero-size on non-master
    }


    const auto& offsets = gi.offsets();  // needed on master only

    Info<< "Using node-comms: " << UPstream::usingNodeComms(comm) << nl;

    const auto interNodeComm = UPstream::commInterNode();
    const auto localNodeComm = UPstream::commLocalNode();

    if (UPstream::usingNodeComms(comm))
    {
        // Stage 0 : The inter-node/intra-node offsets
        labelList interNodeOffsets;
        labelList localNodeOffsets;

        gi.splitNodeOffsets(interNodeOffsets, localNodeOffsets, comm);

        // The first node re-uses the output (allData) when collecting
        // content. All other nodes require temporary node-local storage.

        List<Type> tmpNodeData;
        if (UPstream::is_subrank(interNodeComm))
        {
            tmpNodeData.resize(localNodeOffsets.back());
        }

        List<Type>& nodeData =
        (
            UPstream::master(interNodeComm) ? allData : tmpNodeData
        );


        // Stage 1 : Gather data within the node
        {
            globalIndexGather
            (
                localNodeOffsets,  // (master only)
                localNodeComm,
                UPstream::allProcs(localNodeComm),
                sendData,
                nodeData,
                tag,
                commsType,
                useWindow
            );
        }

        // Stage 2 : Gather data between nodes
        if (UPstream::is_rank(interNodeComm))
        {
            globalIndexGather
            (
                interNodeOffsets,  // (master only)
                interNodeComm,
                UPstream::allProcs(interNodeComm),
                nodeData,
                allData,
                tag,
                commsType,
                useWindow
            );
        }
    }
    else
    {
        globalIndexGather
        (
            offsets,  // needed on master only
            comm,
            UPstream::allProcs(comm),   // All communicator ranks
            sendData,
            allData,
            tag,
            commsType,
            useWindow
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("Set UPstream::debug level");
    argList::addOption("split-size", "NUM", "split with ncores/node");
    argList::addBoolOption("builtin", "only use builtin globalIndex::gather");
    argList::addBoolOption("window", "get data via window");

    // Check -verbose before initialisation
    UPstream::debug = argList::verbose(argc, argv);

    // Check -split-size before initialisation
    {
        int splitSize = option_splitsize(argc, argv);

        if (splitSize >= 0)
        {
            UPstream::nodeCommsControl_ = splitSize;
        }
    }

    #include "setRootCase.H"

    const bool useLocalComms = UPstream::usingNodeComms(UPstream::worldComm);
    bool useWindow = args.found("window");
    bool useBuiltin = args.found("builtin");

    Info<< nl
        << "Getting local-comms: " << Switch::name(useLocalComms) << nl
        << "Getting data with window: " << Switch::name(useWindow) << nl
        << nl;

    if (useWindow && useBuiltin)
    {
        Info<< "Selected '-window' and '-builtin' : ignoring -builtin'"
            << nl;
        useBuiltin = false;
    }

    Random rng(31 + 2*UPstream::myProcNo());

    const label localSize = (5*rng.position<label>(1, 15));

    globalIndex globIndex
    (
        globalIndex::gatherOnly{},
        localSize,
        UPstream::commWorld()
    );

    Info<< "global-index: ";
    printGlobalIndex(Info, globIndex);
    reportOffsets(globIndex);

    Field<scalar> allData;
    Field<scalar> localFld(localSize, scalar(UPstream::myProcNo()));

    if (useBuiltin)
    {
        globIndex.gather
        (
            localFld,
            allData,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking,
            UPstream::commWorld()
        );
    }
    else
    {
        globalIndexGather
        (
            globIndex,
            localFld,
            allData,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking,
            UPstream::commWorld(),
            useWindow
        );
    }

    Pout<< "local: " << flatOutput(localFld) << nl;
    Info<< "field: " << flatOutput(allData) << nl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
