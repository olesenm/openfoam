/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
    Copyright (C) 2026 Keysight Technologies
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
    The gathered data will be a single value constructed from the values
    on individual processors using a user-specified operator.

Note
    Normal gather uses:
    - binary operator that returns a value.
      So assignment that return value to yield the new value

    Combine gather uses:
    - binary operator modifies its first parameter in-place

\*---------------------------------------------------------------------------*/

#include "contiguous.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// Implementation: gather (reduce) single element data onto master
template<class T, class BinaryOp, bool InplaceMode>
void gather_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    T& value,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    [[maybe_unused]]
    const bool printDebug = (Pstream::debug & 2);

    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // if (comms.empty()) return;  // extra safety?
        const auto myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];
        const auto& below = myComm.below();


        // Receive from my downstairs neighbours
        //
        // use element or component type (or byte-wise) for contiguous data

        for (const auto proci : below)
        {
            T received;

            if constexpr (is_contiguous_v<T>)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    &received,
                    std::streamsize(1),  //< count
                    tag,
                    communicator
                );
            }
            else
            {
                IPstream::recv(received, proci, tag, communicator);
            }

            if constexpr (InplaceMode)
            {
                if (printDebug)
                {
                    Perr<< " received from "
                        << proci << " data:" << received << endl;
                }
            }

            if constexpr (InplaceMode)
            {
                // In-place binary operation
                bop(value, received);
            }
            else
            {
                // Assign result of binary operation
                value = bop(value, received);
            }
        }

        // Send up value
        if (const auto above = myComm.above(); above >= 0)
        {
            if constexpr (InplaceMode)
            {
                if (printDebug)
                {
                    Perr<< " sending to " << above
                        << " data:" << value << endl;
                }
            }

            if constexpr (is_contiguous_v<T>)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    above,
                    &value,
                    std::streamsize(1),  //< count
                    tag,
                    communicator
                );
            }
            else
            {
                OPstream::send(value, above, tag, communicator);
            }
        }
    }
}


// Implementation: gather (reduce) single element data onto :masterNo
// using a topo algorithm.
// Return: True if topo algorithm was applied
template<class T, class BinaryOp, bool InplaceMode>
bool gather_topo_algorithm
(
    T& value,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    const bool withTopo =
    (
        UPstream::is_parallel(communicator)
     && UPstream::usingTopoControl(UPstream::topoControls::combine)
     && UPstream::usingNodeComms(communicator)
    );

    if (withTopo)
    {
        // Topological reduce
        // - linear for local-node (assume communication is fast)
        // - tree for inter-node (no assumption about speed)

        using control = std::pair<int, bool>;

        const auto oldWarnComm = UPstream::warnComm;
        const bool usesWarnComm = (communicator == oldWarnComm);

        for
        (
            auto [subComm, linear] :
            {
                // 1: within node
                control{ UPstream::commLocalNode(), true },
                // 2: between nodes
                control{ UPstream::commInterNode(), false }
            }
        )
        {
            if (usesWarnComm) UPstream::commWarn(subComm);

            if (UPstream::is_parallel(subComm))
            {
                PstreamDetail::gather_algorithm<T, BinaryOp, InplaceMode>
                (
                    UPstream::whichCommunication(subComm, linear),
                    value,
                    bop,
                    tag,
                    subComm
                );
            }
        }

        UPstream::commWarn(oldWarnComm);  // Restore warning communicator
    }

    return withTopo;
}

} // End namespace PstreamDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::gather
(
    T& value,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (!InplaceMode && UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiReduce
        (
            &value,
            1,
            UPstream_opType<BinaryOp>::opcode_id,
            communicator
        );
    }
    else if
    (
        !PstreamDetail::gather_topo_algorithm<T, BinaryOp, InplaceMode>
        (
            value,
            bop,
            tag,
            communicator
        )
    )
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        PstreamDetail::gather_algorithm<T, BinaryOp, InplaceMode>
        (
            commOrder,
            value,
            bop,
            tag,
            communicator
        );
    }
}


template<class T, class CombineOp>
void Foam::Pstream::combineGather
(
    T& value,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::gather<T, CombineOp, true>(value, cop, tag, comm);
}


template<class T, class CombineOp>
void Foam::Pstream::combineReduce
(
    T& value,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    if (UPstream::is_parallel(comm))
    {
        // In-place binary operation
        Pstream::gather<T, CombineOp, true>(value, cop, tag, comm);
        Pstream::broadcast(value, comm);
    }
}


// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// Implementation: gather (reduce) list element data onto master

template<class T, class BinaryOp, bool InplaceMode>
void listGather_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    UList<T>& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    [[maybe_unused]]
    const bool printDebug = (Pstream::debug & 2);

    if (!UPstream::is_parallel(communicator) || values.empty())
    {
        // Nothing to do
        return;
    }
    else
    {
        // if (comms.empty()) return;  // extra safety?
        const auto myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];
        const auto& below = myComm.below();

        // Same length on all ranks
        const label listLen = values.size();

        if (printDebug)
        {
            Perr<< "listGather: length=" << listLen << endl;
        }

        List<T> received;

        if (!below.empty())
        {
            // Pre-size for contiguous reading
            if constexpr (is_contiguous_v<T>)
            {
                received.resize_nocopy(listLen);
            }
        }

        // Receive from my downstairs neighbours
        for (const auto proci : below)
        {
            if (printDebug)
            {
                Perr<< "listGather: receive from " << proci << endl;
            }

            if constexpr (is_contiguous_v<T>)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    received,
                    tag,
                    communicator
                );
            }
            else
            {
                received.clear();  // extra safety?
                IPstream::recv(received, proci, tag, communicator);
            }

            if (printDebug)
            {
                Perr<< "listGather: received from "
                    << proci << " data:" << received << endl;
            }

            for (label i = 0; i < listLen; ++i)
            {
                if constexpr (InplaceMode)
                {
                    // In-place binary operation
                    bop(values[i], received[i]);
                }
                else
                {
                    // Assign result of binary operation
                    values[i] = bop(values[i], received[i]);
                }
            }
        }

        // Send up values
        if (const auto above = myComm.above(); above >= 0)
        {
            if (printDebug)
            {
                Perr<< "listGather: sending to " << above
                    << " data:" << values << endl;
            }

            if constexpr (is_contiguous_v<T>)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    above,
                    values,
                    tag,
                    communicator
                );
            }
            else
            {
                OPstream::send(values, above, tag, communicator);
            }
        }
    }
}

} // End namespace PstreamDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::listGather
(
    UList<T>& values,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator) || values.empty())
    {
        // Nothing to do
        return;
    }
    else if constexpr (!InplaceMode && UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiReduce
        (
            values.data(),
            values.size(),  // Same length on all ranks
            UPstream_opType<BinaryOp>::opcode_id,
            communicator
        );
    }
    else if (values.size() == 1)
    {
        // Single value - optimized version
        Pstream::gather<T, BinaryOp, InplaceMode>
        (
            values[0],
            bop,
            tag,
            communicator
        );
    }
    else
    {
        // Multiple values - same length on all ranks
        const label listLen = values.size();

        // Chunk length (== listLen if no chunking is required)
        const label chunkLen = PstreamDetail::listMaxChunkSize(values);

        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        // Operate on chunks, single pass when no chunking is needed
        for (label start = 0; start < listLen; start += chunkLen)
        {
            const auto len = Foam::min(chunkLen, listLen - start);

            auto slice = values.slice(start, len);

            PstreamDetail::listGather_algorithm<T, BinaryOp, InplaceMode>
            (
                commOrder,
                slice,
                bop,
                tag,
                communicator
            );
        }
    }
}


template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::listReduce
(
    UList<T>& values,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int comm
)
{
    if (!UPstream::is_parallel(comm) || values.empty())
    {
        // Nothing to do
    }
    else if constexpr (!InplaceMode && UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiAllReduce
        (
            values.data(),
            values.size(),  // Same length on all ranks
            UPstream_opType<BinaryOp>::opcode_id,
            comm
        );
    }
    else if (values.size() == 1)
    {
        // Single value - optimized version
        Pstream::gather<T, BinaryOp, InplaceMode>(values[0], bop, tag, comm);
        Pstream::broadcast(values[0], comm);
    }
    else
    {
        // Multiple values - same length on all ranks
        Pstream::listGather<T, BinaryOp, InplaceMode>(values, bop, tag, comm);
        PstreamDetail::UList_chunked_broadcast(values, comm);
    }
}


template<class T, class CombineOp>
void Foam::Pstream::listCombineGather
(
    UList<T>& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::listGather<T, CombineOp, true>(values, cop, tag, comm);
}


template<class T, class CombineOp>
void Foam::Pstream::listCombineReduce
(
    UList<T>& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::listReduce<T, CombineOp, true>(values, cop, tag, comm);
}


// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// Implementation: gather (reduce) Map/HashTable containers onto master

template<class Container, class BinaryOp, bool InplaceMode>
void mapGather_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    Container& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    [[maybe_unused]]
    const bool printDebug = (Pstream::debug & 2);

    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // if (comms.empty()) return;  // extra safety?
        const auto myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];
        const auto& below = myComm.below();


        // Receive from my downstairs neighbours
        for (const auto proci : below)
        {
            // Map/HashTable: non-contiguous
            Container received;
            IPstream::recv(received, proci, tag, communicator);

            if constexpr (InplaceMode)
            {
                if (printDebug)
                {
                    Perr<< " received from "
                        << proci << " data:" << received << endl;
                }
            }

            const auto last = received.end();

            for (auto iter = received.begin(); iter != last; ++iter)
            {
                auto slot = values.find(iter.key());

                if (slot.good())
                {
                    // Combine with existing entry

                    if constexpr (InplaceMode)
                    {
                        // In-place binary operation
                        bop(slot.val(), iter.val());
                    }
                    else
                    {
                        // Assign result of binary operation
                        slot.val() = bop(slot.val(), iter.val());
                    }
                }
                else
                {
                    // Create a new entry
                    values.emplace(iter.key(), std::move(iter.val()));
                }

            }
        }

        // Send up values
        if (const auto above = myComm.above(); above >= 0)
        {
            if constexpr (InplaceMode)
            {
                if (printDebug)
                {
                    Perr<< " sending to " << above
                        << " data:" << values << endl;
                }
            }

            OPstream::send(values, above, tag, communicator);
        }
    }
}


// Implementation gather (reduce) Map/HashTable containers onto master
// using a topo algorithm.
// Return: True if topo algorithm was applied
template<class Container, class BinaryOp, bool InplaceMode>
bool mapGather_topo_algorithm
(
    Container& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    const bool withTopo =
    (
        UPstream::is_parallel(communicator)
     && UPstream::usingTopoControl(UPstream::topoControls::mapGather)
     && UPstream::usingNodeComms(communicator)
    );

    if (withTopo)
    {
        // Topological reduce
        // - linear for local-node (assume communication is fast)
        // - tree for inter-node (no assumption about speed)

        using control = std::pair<int, bool>;

        const auto oldWarnComm = UPstream::warnComm;
        const bool usesWarnComm = (communicator == oldWarnComm);

        for
        (
            auto [subComm, linear] :
            {
                // 1: within node
                control{ UPstream::commLocalNode(), true },
                // 2: between nodes
                control{ UPstream::commInterNode(), false }
            }
        )
        {
            if (usesWarnComm) UPstream::commWarn(subComm);

            if (UPstream::is_parallel(subComm))
            {
                PstreamDetail::mapGather_algorithm
                <Container, BinaryOp, InplaceMode>
                (
                    UPstream::whichCommunication(subComm, linear),
                    values,
                    bop,
                    tag,
                    subComm
                );
            }
        }

        UPstream::commWarn(oldWarnComm);  // Restore warning communicator
    }

    return withTopo;
}

} // End namespace PstreamDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class BinaryOp, bool InplaceMode>
void Foam::Pstream::mapGather
(
    Container& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if
    (
        !PstreamDetail::mapGather_topo_algorithm
        <Container, BinaryOp, InplaceMode>
        (
            values,
            bop,
            tag,
            communicator
        )
    )
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        PstreamDetail::mapGather_algorithm
        <Container, BinaryOp, InplaceMode>
        (
            commOrder,
            values,
            bop,
            tag,
            communicator
        );
    }
}


template<class Container, class BinaryOp, bool InplaceMode>
void Foam::Pstream::mapReduce
(
    Container& values,
    BinaryOp bop,
    const int tag,
    const int comm
)
{
    Pstream::mapGather<Container, BinaryOp, InplaceMode>
    (
        values, bop, tag, comm
    );
    Pstream::broadcast(values, comm);
}


template<class Container, class CombineOp>
void Foam::Pstream::mapCombineGather
(
    Container& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::mapGather<Container, CombineOp, true>
    (
        values, cop, tag, comm
    );
}


template<class Container, class CombineOp>
void Foam::Pstream::mapCombineReduce
(
    Container& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::mapReduce<Container, CombineOp, true>
    (
        values, cop, tag, comm
    );
}


// ************************************************************************* //
