/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"
#include <functional>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class IntType>
Foam::List<IntType>
Foam::globalIndex::calcRecvSizes(IntType localSize, int communicator)
{
    static_assert(std::is_integral_v<IntType>, "Integral required");

    List<IntType> count(UPstream::listGatherValues(localSize, communicator));
    // Replace element [0] with the max receive count
    if (count.size() > 1)
    {
        count[0] = *(std::max_element(count.begin()+1, count.end()));
    }
    return count;
}


// Helpers

template<class Addr>
Foam::labelList
Foam::globalIndex::calcOffsets
(
    const IndirectListBase<label, Addr>& counts,
    const bool checkOverflow
)
{
    labelList values;

    if (const label len = counts.size(); len > 0)
    {
        values.resize(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];

            values[i] = start;
            start += count;

            if (checkOverflow && start < values[i])
            {
                reportOverflowAndExit(i, values[i], count);
            }
        }
        values[len] = start;
    }

    return values;
}


template<class SubListType>
Foam::labelList
Foam::globalIndex::calcListOffsets
(
    const List<SubListType>& lists,
    const bool checkOverflow
)
{
    labelList values;

    if (const label len = lists.size(); len > 0)
    {
        values.resize(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = lists[i].size();

            values[i] = start;
            start += count;

            if (checkOverflow && start < values[i])
            {
                reportOverflowAndExit(i, values[i], count);
            }
        }
        values[len] = start;
    }

    return values;
}


// Low-level

template<class ProcIDsContainer, class Type>
Foam::List<Type> Foam::globalIndex::listGatherValues
(
    const label comm,
    const ProcIDsContainer& procIDs,
    const Type& localValue,
    const int tag,
    UPstream::commsTypes commsType
)
{
    // low-level: no parRun guard?
    const int masterProci = (procIDs.empty() ? 0 : procIDs[0]);

    // if (!UPstream::is_parallel(comm))
    // {
    //     List<Type> allValues(1);
    //     allValues[0] = localValue;
    //     return allValues;
    // }

    List<Type> allValues;

    // Cannot use non-blocking for non-contiguous data
    if constexpr (!is_contiguous_v<Type>)
    {
        commsType = UPstream::commsTypes::scheduled;
    }


    const auto startOfRequests = UPstream::nRequests();

    if (UPstream::myProcNo(comm) == masterProci)
    {
        allValues.resize_nocopy(procIDs.size());
        allValues[0] = localValue;

        for (label i = 1; i < procIDs.size(); ++i)
        {
            if constexpr (is_contiguous_v<Type>)
            {
                UIPstream::read
                (
                    commsType,
                    procIDs[i],
                    reinterpret_cast<char*>(&allValues[i]),
                    sizeof(Type),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream::recv(allValues[i], procIDs[i], tag, comm);
            }
        }
    }
    else
    {
        if constexpr (is_contiguous_v<Type>)
        {
            UOPstream::write
            (
                commsType,
                masterProci,
                reinterpret_cast<const char*>(&localValue),
                sizeof(Type),
                tag,
                comm
            );
        }
        else
        {
            OPstream::send(localValue, commsType, masterProci, tag, comm);
        }
    }

    // Process sync
    UPstream::waitRequests(startOfRequests);

    return allValues;
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::gather
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& fld,
    UList<Type>& allFld,    // must be adequately sized on master
    const int tag,
    UPstream::commsTypes commsType
)
{
    // low-level: no parRun guard
    const int masterProci = (procIDs.empty() ? 0 : procIDs[0]);

    // Cannot use non-blocking for non-contiguous data
    if constexpr (!is_contiguous_v<Type>)
    {
        commsType = UPstream::commsTypes::scheduled;
    }

    const auto startOfRequests = UPstream::nRequests();

    if (UPstream::myProcNo(comm) == masterProci)
    {
        if (FOAM_UNLIKELY(allFld.size() < off.back()))  // ie, totalSize()
        {
            FatalErrorInFunction
                << "[out] UList size=" << allFld.size()
                << " too small to receive " << off.back() << nl
                << Foam::abort(FatalError);
        }

        for (label i = 1; i < procIDs.size(); ++i)
        {
            SubList<Type> slot(allFld, off[i+1]-off[i], off[i]);

            if (slot.empty())
            {
                // Nothing to do
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
    }
    else
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

    // Process sync
    UPstream::waitRequests(startOfRequests);
}


template<class ProcIDsContainer, class Type, class Addr>
void Foam::globalIndex::gather
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const IndirectListBase<Type, Addr>& fld,
    UList<Type>& allFld,    // must be adequately sized on master
    const int tag,
    UPstream::commsTypes commsType
)
{
    // low-level: no parRun guard
    const int masterProci = (procIDs.empty() ? 0 : procIDs[0]);

    if constexpr (is_contiguous_v<Type>)
    {
        if (commsType == UPstream::commsTypes::nonBlocking)
        {
            // Contiguous data and requested nonBlocking.
            //
            // Flatten list (locally) so that we can benefit from using
            // direct read/write of contiguous data

            List<Type> flattened(fld);

            gather
            (
                off,
                comm,
                procIDs,
                flattened,
                allFld,
                tag,
                commsType
            );
            return;
        }
    }


    // Non-contiguous is always non-blocking

    if (UPstream::myProcNo(comm) == masterProci)
    {
        if (FOAM_UNLIKELY(allFld.size() < off.back()))  // ie, totalSize()
        {
            FatalErrorInFunction
                << "[out] UList size=" << allFld.size()
                << " too small to receive " << off.back() << nl
                << Foam::abort(FatalError);
        }

        for (label i = 1; i < procIDs.size(); ++i)
        {
            SubList<Type> slot(allFld, off[i+1]-off[i], off[i]);

            if (slot.empty())
            {
                // Nothing to do
            }
            else
            {
                IPstream::recv(slot, procIDs[i], tag, comm);
            }
        }

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied
        {
            SubList<Type> dst(allFld, off[1]-off[0], off[0]);

            if (!dst.empty() && (dst.size() == fld.size()))
            {
                dst.deepCopy(fld);
            }
        }
    }
    else
    {
        if (fld.empty())
        {
            // Nothing to do
        }
        else
        {
            OPstream::send(fld, masterProci, tag, comm);
        }
    }
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::gatherInplace
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    List<Type>& fld,
    const int tag,
    UPstream::commsTypes commsType
)
{
    if (!UPstream::is_parallel(comm))
    {
        // Serial: (no-op)
        return;
    }

    const bool master =
    (
        UPstream::myProcNo(comm) == (procIDs.empty() ? 0 : procIDs[0])
    );

    List<Type> allData;
    if (master)
    {
        allData.resize_nocopy(off.back());  // == totalSize()
    }

    globalIndex::gather(off, comm, procIDs, fld, allData, tag, commsType);

    if (master)
    {
        fld = std::move(allData);
    }
    else
    {
        fld.clear();  // zero-size on non-master
    }
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::gather
(
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& fld,
    List<Type>& allData,
    const int tag,
    UPstream::commsTypes commsType
) const
{
    if (!UPstream::is_parallel(comm))
    {
        // Serial: (no-op)
        return;
    }

    if (UPstream::myProcNo(comm) == (procIDs.empty() ? 0 : procIDs[0]))
    {
        // presize => totalSize()
        allData.resize_nocopy(totalSize());
    }
    else
    {
        allData.clear();  // zero-size on non-master
    }

    globalIndex::gather(offsets(), comm, procIDs, fld, allData, tag, commsType);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::globalIndex::gather
(
    const UList<Type>& sendData,
    List<Type>& allData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
) const
{
    if (!UPstream::parRun())
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }

    if (UPstream::master(comm))
    {
        allData.resize_nocopy(totalSize());
    }
    else
    {
        allData.clear();  // zero-size on non-master
    }

    if (!UPstream::usingNodeComms(comm))
    {
        globalIndex::gather
        (
            offsets(),                  // needed on master only
            comm,
            UPstream::allProcs(comm),   // All communicator ranks
            sendData,
            allData,
            tag,
            commsType
        );
    }
    else
    {
        // Using node-based hierarchy

        // Using comm-world and have node communication active
        const auto interNodeComm = UPstream::commInterNode();
        const auto localNodeComm = UPstream::commLocalNode();

        // Stage 0 : The inter-node/intra-node offsets
        labelList interNodeOffsets;
        labelList localNodeOffsets;
        this->splitNodeOffsets(interNodeOffsets, localNodeOffsets, comm);

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

        const auto oldWarnComm = UPstream::warnComm;
        const bool usesWarnComm = (comm == oldWarnComm);

        // Stage 1 : Gather data within the node
        {
            if (usesWarnComm) UPstream::commWarn(localNodeComm);

            globalIndex::gather
            (
                localNodeOffsets,  // (master only)
                localNodeComm,
                UPstream::allProcs(localNodeComm),
                sendData,
                nodeData,  // node-local dest (or the allData parameter)
                tag,
                commsType
            );
        }

        // Stage 2 : Gather data between nodes
        if (UPstream::is_rank(interNodeComm))
        {
            if (usesWarnComm) UPstream::commWarn(interNodeComm);

            globalIndex::gather
            (
                interNodeOffsets,  // (master only)
                interNodeComm,
                UPstream::allProcs(interNodeComm),
                nodeData,
                allData,
                tag,
                commsType
            );
        }

        UPstream::commWarn(oldWarnComm);  // Restore warning communicator
    }
}


template<class Type, class Addr>
void Foam::globalIndex::gather
(
    const IndirectListBase<Type, Addr>& sendData,
    List<Type>& allData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
) const
{
    if (!UPstream::parRun())
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }
    else if constexpr (is_contiguous_v<Type>)
    {
        if (commsType == UPstream::commsTypes::nonBlocking)
        {
            // Contiguous data and requested nonBlocking.
            //
            // Flatten list (locally) so that we can benefit from using
            // direct read/write of contiguous data

            List<Type> flattened(sendData);

            this->gather
            (
                flattened,
                allData,
                tag,
                commsType,
                comm
            );
            return;
        }
    }

    if (UPstream::master(comm))
    {
        allData.resize_nocopy(totalSize());
    }
    else
    {
        allData.clear();  // zero-size on non-master
    }

    {
        globalIndex::gather
        (
            offsets(),                  // needed on master only
            comm,
            UPstream::allProcs(comm),   // All communicator ranks
            sendData,
            allData,
            tag,
            commsType
        );
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::gather
(
    const UList<Type>& sendData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
) const
{
    OutputContainer allData;
    this->gather(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type, class Addr, class OutputContainer>
OutputContainer Foam::globalIndex::gather
(
    const IndirectListBase<Type, Addr>& sendData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
) const
{
    OutputContainer allData;
    this->gather(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type>
void Foam::globalIndex::gatherInplace
(
    List<Type>& fld,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
) const
{
    if (UPstream::parRun())
    {
        List<Type> allData;
        this->gather(fld, allData, tag, commsType, comm);

        if (UPstream::master(comm))
        {
            fld = std::move(allData);
        }
        else
        {
            fld.clear();  // zero-size on non-master
        }
    }
    // Serial: (no-op)
}


template<class Type, class OutputContainer>
void Foam::globalIndex::mpiGather
(
    const UList<Type>& sendData,
    OutputContainer& allData,
    const label comm,
    [[maybe_unused]] UPstream::commsTypes commsType,
    [[maybe_unused]] const int tag
) const
{
    if (!UPstream::is_parallel(comm))
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }

    if (UPstream::master(comm))
    {
        allData.resize_nocopy(totalSize());
    }
    else
    {
        allData.clear();  // zero-size on non-master
    }

    if constexpr (UPstream_dataType<Type>::value)
    {
        // Restrict to basic (or aliased) MPI types
        // - simplifies calculating counts/offsets
        //   and can call UPstream::mpiGatherv directly

        // The sizing factor is constexpr
        constexpr std::streamsize count = UPstream_dataType<Type>::size(1);

        static_assert
        (
            (count == 1),
            "Code does not (yet) work with aggregate types"
        );

        List<int> recvCounts;
        List<int> recvOffsets;

        if (UPstream::master(comm))
        {
            // Must be same as Pstream::nProcs(comm), at least on master!
            // if (UPstream::nProcs(comm) != this->nProcs()) ...

            recvCounts.resize(offsets().size()-1);
            recvOffsets.resize(offsets().size());

            // Copy offsets
            std::copy(offsets().begin(), offsets().end(), recvOffsets.begin());

            // Calculate sizes
            std::transform
            (
                offsets().begin()+1, offsets().end(),
                offsets().begin(), recvCounts.begin(),
                std::minus<>{}
            );

            // FUTURE .. fix sizes and offsets by the element factor...
            // if constexpr (UPstream_basic_dataType<Type>::size(1) > 1)
        }

        int sendSize = static_cast<int>(sendData.size());

        // Note we let MPI_Gatherv copy back the local data as well...

        UPstream::mpiGatherv
        (
            sendData.cdata(),
            sendSize,
            allData.data(),
            recvCounts,
            recvOffsets,
            comm
        );
    }
    else
    {
        // Regular (manual) gathering
        globalIndex::gather
        (
            offsets(),                  // needed on master only
            comm,
            UPstream::allProcs(comm),   // All communicator ranks
            sendData,
            allData,
            tag,
            commsType
        );
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::mpiGather
(
    const UList<Type>& sendData,
    const label comm,
    UPstream::commsTypes commsType,
    const int tag
) const
{
    OutputContainer allData;
    mpiGather(sendData, allData, comm, commsType, tag);
    return allData;
}


template<class Type>
void Foam::globalIndex::mpiGatherInplace
(
    List<Type>& fld,
    const label comm,
    UPstream::commsTypes commsType,
    const int tag
) const
{
    if (UPstream::parRun())
    {
        List<Type> allData;
        mpiGather(fld, allData, comm, commsType, tag);

        if (UPstream::master(comm))
        {
            fld = std::move(allData);
        }
        else
        {
            fld.clear();  // zero-size on non-master
        }
    }
    // Serial: (no-op)
}


template<class Type, class OutputContainer>
void Foam::globalIndex::mpiGatherOp
(
    const UList<Type>& sendData,
    OutputContainer& allData,
    const label comm,
    UPstream::commsTypes commsType,
    const int tag
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, sendData.size(), comm)
            .mpiGather(sendData, allData, comm, commsType, tag);
    }
    else
    {
        // Serial: direct copy
        allData = sendData;
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::mpiGatherOp
(
    const UList<Type>& sendData,
    const label comm,
    UPstream::commsTypes commsType,
    const int tag
)
{
    OutputContainer allData;
    mpiGatherOp(sendData, allData, comm, commsType, tag);
    return allData;
}


template<class Type>
void Foam::globalIndex::mpiGatherInplaceOp
(
    List<Type>& fld,
    const label comm,
    UPstream::commsTypes commsType,
    const int tag
)
{
    if (UPstream::parRun())
    {
        List<Type> allData;
        mpiGatherOp(fld, allData, comm, commsType, tag);

        if (UPstream::master(comm))
        {
            fld = std::move(allData);
        }
        else
        {
            fld.clear();  // zero-size on non-master
        }
    }
    // Serial: (no-op)
}


template<class Type>
void Foam::globalIndex::gatherOp
(
    const UList<Type>& sendData,
    List<Type>& allData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, sendData.size(), comm)
            .gather(sendData, allData, tag, commsType, comm);
    }
    else
    {
        // Serial: direct copy
        allData = sendData;
    }
}


template<class Type, class Addr>
void Foam::globalIndex::gatherOp
(
    const IndirectListBase<Type, Addr>& sendData,
    List<Type>& allData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, sendData.size(), comm)
            .gather(sendData, allData, tag, commsType, comm);
    }
    else
    {
        // Serial: direct copy
        allData = List<Type>(sendData);
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::gatherOp
(
    const UList<Type>& sendData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
)
{
    OutputContainer allData;
    gatherOp(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type, class Addr, class OutputContainer>
OutputContainer Foam::globalIndex::gatherOp
(
    const IndirectListBase<Type, Addr>& sendData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
)
{
    OutputContainer allData;
    gatherOp(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type>
void Foam::globalIndex::gatherInplaceOp
(
    List<Type>& fld,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, fld.size(), comm)
            .gatherInplace(fld, tag, commsType, comm);
    }
    // Serial: (no-op)
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::scatter
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& allFld,
    UList<Type>& fld,
    const int tag,
    UPstream::commsTypes commsType
)
{
    // low-level: no parRun guard
    const int masterProci = (procIDs.empty() ? 0 : procIDs[0]);

    // Cannot use non-blocking for non-contiguous data
    if constexpr (!is_contiguous_v<Type>)
    {
        commsType = UPstream::commsTypes::scheduled;
    }


    const auto startOfRequests = UPstream::nRequests();

    if (UPstream::myProcNo(comm) == masterProci)
    {
        for (label i = 1; i < procIDs.size(); ++i)
        {
            const SubList<Type> slot(allFld, off[i+1]-off[i], off[i]);

            if (slot.empty())
            {
                // Nothing to do
            }
            else if constexpr (is_contiguous_v<Type>)
            {
                UOPstream::write
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
                OPstream::send(slot, commsType, procIDs[i], tag, comm);
            }
        }

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied.
        // Also handle the case where we have a slice of the full
        // list.

        {
            SubList<Type> dst(fld, off[1]-off[0]);
            SubList<Type> src(allFld, off[1]-off[0], off[0]);

            if (!dst.empty() && (dst.data() != src.data()))
            {
                dst = src;
            }
        }
    }
    else
    {
        // Note: we are receiving into UList, so sizes MUST match or we
        // have a problem. Can therefore reasonably assume that a zero-sized
        // send matches a zero-sized receive, and we can skip that.

        if (fld.empty())
        {
            // Nothing to do
        }
        else if constexpr (is_contiguous_v<Type>)
        {
            UIPstream::read
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
            IPstream::recv(fld, masterProci, tag, comm);
        }
    }

    // Process sync
    UPstream::waitRequests(startOfRequests);
}


template<class Type>
void Foam::globalIndex::scatter
(
    const UList<Type>& allData,
    UList<Type>& localData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
) const
{
    if (UPstream::parRun())
    {
        scatter
        (
            offsets(),                 // needed on master only
            comm,
            UPstream::allProcs(comm),  // All communicator ranks
            allData,
            localData,
            tag,
            commsType
        );
    }
    else
    {
        // Serial: direct copy
        // - fails miserably if incorrectly dimensioned!
        localData.deepCopy(allData);
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::scatter
(
    const UList<Type>& allData,
    const int tag,
    UPstream::commsTypes commsType,
    const label comm
) const
{
    if (UPstream::parRun())
    {
        // The globalIndex might be correct on master only,
        // so scatter local sizes to ensure consistency

        const label count
        (
            UPstream::listScatterValues<label>(this->localSizes(), comm)
        );

        OutputContainer localData;
        localData.resize(count);
        this->scatter(allData, localData, tag, commsType, comm);

        return localData;
    }
    else
    {
        // Serial: direct copy
        return OutputContainer(allData);
    }
}


template<class Type, class CombineOp>
void Foam::globalIndex::get
(
    List<Type>& allFld,
    const labelUList& globalIds,
    const CombineOp& cop,
    const label comm,
    const int tag
) const
{
    allFld.resize_nocopy(globalIds.size());

    if (globalIds.size())
    {
        // Sort according to processor
        labelList order;
        DynamicList<label> validBins(Pstream::nProcs());

        CompactListList<label> bins
        (
            bin(offsets(), globalIds, order, validBins)
        );

        // Send local indices to individual processors as local index
        PstreamBuffers sendBufs(comm, tag);

        for (const auto proci : validBins)
        {
            labelList localIDs(bins[proci]);

            for (label& val : localIDs)
            {
                val = toLocal(proci, val);
            }

            UOPstream os(proci, sendBufs);
            os << localIDs;
        }
        sendBufs.finishedSends();


        PstreamBuffers returnBufs(comm, tag);

        for (const int proci : sendBufs.allProcs())
        {
            if (sendBufs.recvDataCount(proci))
            {
                UIPstream is(proci, sendBufs);
                labelList localIDs(is);

                // Collect entries
                List<Type> fld(localIDs.size());
                cop(fld, localIDs);

                UOPstream os(proci, returnBufs);
                os << fld;
            }
        }
        returnBufs.finishedSends();

        // Slot back
        for (const auto proci : validBins)
        {
            label start = bins.offsets()[proci];
            const SubList<label> es
            (
                order,
                bins.offsets()[proci+1]-start,  // start
                start
            );
            UIPstream is(proci, returnBufs);
            List<Type> fld(is);

            UIndirectList<Type>(allFld, es) = fld;
        }
    }
}


// ************************************************************************* //
