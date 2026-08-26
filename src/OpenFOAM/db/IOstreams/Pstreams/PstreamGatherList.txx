/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    values[UPstream::myProcNo(comm)].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use broadcast to distribute the data.

\*---------------------------------------------------------------------------*/

#include "contiguous.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// Implementation: gather data, keeping individual values separate.
// Output is only valid (consistent) on master
template<class T>
void gatherList_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    UList<T>& values,
    const int tag,
    const int communicator
)
{
    const bool printDebug = (Pstream::debug & 2);

    if (FOAM_UNLIKELY(!UPstream::is_parallel(communicator)))
    {
        // Nothing to do
        return;
    }
    else
    {
        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        // if (comms.empty()) return;  // extra safety?
        const auto myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];


        // Local buffer for send/recv of contiguous
        [[maybe_unused]] DynamicList<T> buffer;

        // Presize buffer
        if constexpr (is_contiguous_v<T>)
        {
            label maxCount = 0;

            for (const auto belowID : myComm.below())
            {
                auto count = comms[belowID].allBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            if (myComm.above() >= 0)
            {
                auto count = myComm.allBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            buffer.reserve(maxCount + 1);
        }


        // Receive from my downstairs neighbours
        for (const auto belowID : myComm.below())
        {
            const auto& leaves = comms[belowID].allBelow();

            if constexpr (is_contiguous_v<T>)
            {
                if (leaves.empty())
                {
                    // Receive directly into destination
                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        belowID,
                       &(values[belowID]),
                        1,
                        tag,
                        communicator
                    );
                }
                else
                {
                    // Receive via intermediate buffer
                    buffer.resize_nocopy(leaves.size() + 1);

                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        belowID,
                        buffer,
                        tag,
                        communicator
                    );

                    label recvIdx(0);
                    values[belowID] = buffer[recvIdx++];

                    for (const auto leafID : leaves)
                    {
                        values[leafID] = buffer[recvIdx++];
                    }
                }
            }
            else
            {
                IPstream fromBelow(belowID, tag, communicator);
                fromBelow >> values[belowID];

                if (printDebug)
                {
                    Perr<< " received through "
                        << belowID << " data from:" << belowID
                        << " data:" << values[belowID] << endl;
                }

                // Receive from all other processors below belowID
                for (const auto leafID : leaves)
                {
                    fromBelow >> values[leafID];

                    if (printDebug)
                    {
                        Perr<< " received through "
                            << belowID << " data from:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }

        // Send up from values:
        // - my own value first
        // - all belowLeaves next
        if (const auto above = myComm.above(); above >= 0)
        {
            const auto& leaves = myComm.allBelow();

            if (printDebug)
            {
                Perr<< " sending to " << above
                    << " data from me:" << myProci
                    << " data:" << values[myProci] << endl;
            }

            if constexpr (is_contiguous_v<T>)
            {
                if (leaves.empty())
                {
                    // Send directly
                    UOPstream::write
                    (
                        UPstream::commsTypes::scheduled,
                        above,
                       &(values[myProci]),
                        1,
                        tag,
                        communicator
                    );
                }
                else
                {
                    // Send via intermediate buffer
                    buffer.resize_nocopy(leaves.size() + 1);

                    label sendIdx(0);
                    buffer[sendIdx++] = values[myProci];

                    for (const auto leafID : leaves)
                    {
                        buffer[sendIdx++] = values[leafID];
                    }

                    UOPstream::write
                    (
                        UPstream::commsTypes::scheduled,
                        above,
                        buffer,
                        tag,
                        communicator
                    );
                }
            }
            else
            {
                OPstream toAbove(above, tag, communicator);
                toAbove << values[myProci];

                for (const auto leafID : leaves)
                {
                    if (printDebug)
                    {
                        Perr<< " sending to "
                            << above << " data from:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                    toAbove << values[leafID];
                }
            }
        }
    }
}


// Implementation: gather data, keeping individual values separate
// using topo algorithm
// Return: True if topo algorithm was applied
template<class T>
bool gatherList_topo_algorithm
(
    UList<T>& values,
    const int tag,
    const int communicator
)
{
    const bool withTopo =
    (
        UPstream::is_parallel(communicator)
     && UPstream::usingTopoControl(UPstream::topoControls::gatherList)
     && UPstream::usingNodeComms(communicator)
    );

    if (withTopo)
    {
        // Topological gathering

        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        // Overall node-wise offsets
        const auto& off = UPstream::interNode_offsets();

        // The per-node processor range
        const auto& nodeProcs = UPstream::localNode_parentProcs();

        // The per-node sub-section of values
        auto nodeValues = values.slice(nodeProcs.start(), nodeProcs.size());

        const auto oldWarnComm = UPstream::warnComm;
        const bool usesWarnComm = (communicator == oldWarnComm);

        // Stage 1: gather values within a node
        // - linear for local-node (assume communication is fast)

        if
        (
            const auto subComm = UPstream::commLocalNode();
            UPstream::is_parallel(subComm)
        )
        {
            constexpr bool linear(true);

            if (usesWarnComm) UPstream::commWarn(subComm);

            PstreamDetail::gatherList_algorithm<T>
            (
                UPstream::whichCommunication(subComm, linear),
                nodeValues,
                tag,
                subComm
            );
        }

        // Stage 2: gather between node leaders
        // - this unfortunately corresponds to a gatherv process
        //   (number of cores per node is not identical)
        // - code strongly resembles globalIndex::gather

        if
        (
            const auto subComm = UPstream::commInterNode();
            UPstream::is_parallel(subComm)
        )
        {
            if (usesWarnComm) UPstream::commWarn(subComm);

            if (UPstream::master(subComm))
            {
                for (auto proci : UPstream::subProcs(subComm))
                {
                    auto slot =
                        values.slice(off[proci], off[proci+1]-off[proci]);

                    // Probably not contiguous though,
                    // otherwise would have used mpiGather()

                    if constexpr (is_contiguous_v<T>)
                    {
                        UIPstream::read
                        (
                            UPstream::commsTypes::scheduled,
                            proci,
                            slot,
                            tag,
                            subComm
                        );
                    }
                    else
                    {
                        IPstream::recv(slot, proci, tag, subComm);
                    }
                }
            }
            else
            {
                if constexpr (is_contiguous_v<T>)
                {
                    UOPstream::write
                    (
                        UPstream::commsTypes::scheduled,
                        UPstream::masterNo(),
                        nodeValues,
                        tag,
                        subComm
                    );
                }
                else
                {
                    OPstream::send
                    (
                        nodeValues,
                        UPstream::commsTypes::scheduled,
                        UPstream::masterNo(),
                        tag,
                        subComm
                    );
                }
            }
        }

        UPstream::commWarn(oldWarnComm);  // Restore warning communicator
    }

    return withTopo;
}


// Implementation: inverse of gatherList_algorithm
template<class T>
void scatterList_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    UList<T>& values,
    const int tag,
    const int communicator
)
{
    const bool printDebug = (Pstream::debug & 2);

    if (FOAM_UNLIKELY(!UPstream::is_parallel(communicator)))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Apart from the additional size check, the only difference
        // between scatterList() and using broadcast(List<T>&) or a regular
        // scatter(List<T>&) is that processor-local data is skipped.

        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        // if (comms.empty()) return;  // extra safety?
        const auto myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];


        // Local buffer for send/recv of contiguous
        [[maybe_unused]] DynamicList<T> buffer;

        // Presize buffer
        if constexpr (is_contiguous_v<T>)
        {
            label maxCount = 0;

            if (myComm.above() >= 0)
            {
                auto count = myComm.allNotBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            for (const auto belowID : myComm.below())
            {
                auto count = comms[belowID].allNotBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            buffer.reserve(maxCount);
        }


        // Receive from up
        if (const auto above = myComm.above(); above >= 0)
        {
            const auto& leaves = myComm.allNotBelow();

            if constexpr (is_contiguous_v<T>)
            {
                buffer.resize_nocopy(leaves.size());

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    above,
                    buffer,
                    tag,
                    communicator
                );

                label recvIdx(0);
                for (const auto leafID : leaves)
                {
                    values[leafID] = buffer[recvIdx++];
                }
            }
            else
            {
                IPstream fromAbove(above, tag, communicator);

                for (const auto leafID : leaves)
                {
                    fromAbove >> values[leafID];

                    if (printDebug)
                    {
                        Perr<< " received through "
                            << above << " data for:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }

        // Send to my downstairs neighbours
        forAllReverse(myComm.below(), belowI)
        {
            const auto belowID = myComm.below()[belowI];
            const auto& leaves = comms[belowID].allNotBelow();

            if constexpr (is_contiguous_v<T>)
            {
                buffer.resize_nocopy(leaves.size());

                label sendIdx(0);
                for (const auto leafID : leaves)
                {
                    buffer[sendIdx++] = values[leafID];
                }

                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    buffer,
                    tag,
                    communicator
                );
            }
            else
            {
                OPstream toBelow(belowID, tag, communicator);

                // Send data destined for all other processors below belowID
                for (const auto leafID : leaves)
                {
                    toBelow << values[leafID];

                    if (printDebug)
                    {
                        Perr<< " sent through "
                            << belowID << " data for:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }
    }
}

} // End namespace PstreamDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::Pstream::gatherList
(
    UList<T>& values,
    [[maybe_unused]] const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (is_contiguous_v<T>)
    {
        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        // In-place gather for contiguous types - one element per rank
        // - all ranks:
        //   * send pointer is source location from within the list
        // - on master:
        //   * recv pointer is the full list (same address as first location)
        // - on rank:
        //   * recv pointer is irrelevant
        //
        // So can simply use identical pointers for send/recv

        auto* ptr = values.data() + UPstream::myProcNo(communicator);
        UPstream::mpiGather(ptr, ptr, 1, communicator);
    }
    else if
    (
        !PstreamDetail::gatherList_topo_algorithm
        (
            values,
            tag,
            communicator
        )
    )
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        PstreamDetail::gatherList_algorithm
        (
            commOrder, values, tag, communicator
        );
    }
}


template<class T>
void Foam::Pstream::allGatherList
(
    UList<T>& values,
    [[maybe_unused]] const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (is_contiguous_v<T>)
    {
        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values is too small:" << values.size()
                << " vs numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        // Allgather for contiguous types - one element per rank
        UPstream::mpiAllGather(values.data(), 1, communicator);
    }
    else
    {
        // IMPORTANT: always call the *_algorithm() versions here and
        // never the base versions [eg, Pstream::gatherList()] since
        // the communcation order must be absolutely identical
        // for the gatherList and scatterList, otherwise the results will
        // not replicate the allGather behaviour.
        //
        // This also means that we must avoid the gatherList_topo_algorithm()
        // as well, since this does not pair well with scatterList_algorithm()

        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        PstreamDetail::gatherList_algorithm
        (
            commOrder, values, tag, communicator
        );
        PstreamDetail::scatterList_algorithm
        (
            commOrder, values, tag, communicator
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Housekeeping

template<class T>
void Foam::Pstream::scatterList
(
    UList<T>& values,
    const int tag,
    const int communicator
)
{
    if (UPstream::is_parallel(communicator))
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        PstreamDetail::scatterList_algorithm
        (
            commOrder, values, tag, communicator
        );
    }
}


// ************************************************************************* //
