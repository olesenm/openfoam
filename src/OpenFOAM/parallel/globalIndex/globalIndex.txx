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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class IntType>
void Foam::GlobalIndex<IntType>::reset_nonParallel
(
    label localSize,
    int communicator
)
{
    auto& starts = addr_type::offsets();

    if
    (
        label numProc = 0;
        (communicator >= 0 && (numProc = UPstream::nProcs(communicator)) > 0)
    )
    {
        // localSize on-proc, zero elsewhere:
        // offsets are 0 up to and including proci, then localSize
        starts.resize_fill(numProc+1, 0);

        if
        (
            label proci = UPstream::myProcNo(communicator);
            (proci >= 0 && proci < numProc)
        )
        {
            starts.slice(proci+1) = localSize;
        }
    }
    else
    {
        starts.clear();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
void Foam::GlobalIndex<IntType>::reset
(
    label localSize,
    int communicator,
    bool parallel
)
{
    if (parallel && UPstream::is_parallel(communicator))
    {
        constexpr bool checkOverflow = (sizeof(IntType) < sizeof(int64_t));

        reset
        (
            UPstream::allGatherValues(localSize, communicator),
            checkOverflow
        );
    }
    else
    {
        reset_nonParallel(localSize, communicator);
    }
}


template<class IntType>
void Foam::GlobalIndex<IntType>::reset
(
    globalIndexBase::gatherOnly,
    label localSize,
    int communicator,
    bool parallel
)
{
    if (!parallel)
    {
        addr_type::single(localSize);
        return;
    }

    // one-sided gather only
    reset(UPstream::listGatherValues(localSize, communicator));
}


template<class IntType>
void Foam::GlobalIndex<IntType>::reset
(
    globalIndex::gatherTotal,
    label localSize,
    int communicator,
    bool parallel
)
{
    if (!parallel)
    {
        addr_type::single(localSize);
        return;
    }

    // one-sided gather only
    reset(UPstream::listGatherValues(localSize, communicator));

    if (UPstream::is_parallel(communicator))
    {
        IntType total = addr_type::totalSize();

        UPstream::broadcast(&total, 1, communicator);

        if (UPstream::is_subrank(communicator))
        {
            // total value only (no content for rank 0)
            addr_type::sentinel(total);
        }
    }
}


template<class IntType>
void Foam::GlobalIndex<IntType>::reset
(
    globalIndexBase::gatherOffset,
    label localSize,
    int communicator,
    bool parallel
)
{
    if (!parallel)
    {
        addr_type::single(localSize);
        return;
    }

    if (UPstream::is_parallel(communicator))
    {
        const auto numProc = UPstream::nProcs(communicator);

        // Number of gatherOnly counts (master: numProc, sub-rank: 0)
        const int numCounts
        (
            UPstream::is_master(communicator) ? static_cast<int>(numProc) : 0
        );

        // Work space : allocated with double space! (still 0 on non-master).
        // - managed as an "uint64_t" element type to ensure good alignment
        //   whilst avoiding aligned allocator (fragile on mingw)

        std::unique_ptr<uint64_t[]> workBuffer;

        if (numCounts)
        {
            workBuffer = std::make_unique<uint64_t[]>(2*numCounts);
        }

        // The full content as <IntType> for scattering (start, total) tuples
        UList<IntType> tuples
        (
            reinterpret_cast<IntType*>(workBuffer.get()),
            (2*numCounts)
        );

        // The initial slice as <label> for getting/setting the counts
        labelUList counts
        (
            reinterpret_cast<label*>(workBuffer.get()),
            numCounts
        );

        // Gather the local sizes (master only)
        UPstream::mpiGather(&localSize, counts.data(), 1, communicator);

        // On the master: from local sizes to offsets and subsequently
        // pack (start, total) per rank for scattering

        if (UPstream::is_master(communicator))
        {
            // Calculate offsets from the collected local sizes
            reset(counts);

            // Populate the work array with (offset,total) tuples for
            // scattering to sub-ranks.
            // Note: also scattering copies of 'total' everywhere instead
            // using an additional broadcast.

            const IntType total = totalSize();

            auto iter = offsets().cbegin();

            for (label proci = 0; proci < numProc; ++proci)
            {
                tuples[2*proci]   = *iter++;
                tuples[2*proci+1] = total;
            }
        }

        IntType dims[2];  // (start, total) tuple

        UPstream::mpiScatter(tuples.cdata(), dims, 2, communicator);

        if (UPstream::is_subrank(communicator))
        {
            // Store either as (beg len total) tuple, or as (beg len) tuple
            // if the is no trailing content
            addr_type::single(dims[0], localSize, dims[1]);
        }
    }
    else
    {
        reset_nonParallel(localSize, communicator);
    }
}


// * * * * * * * * * * * * * * General Aggregation * * * * * * * * * * * * * //

template<class IntType>
Foam::List<IntType>
Foam::globalIndexBase::calcRecvSizes(IntType localSize, int communicator)
{
    static_assert(std::is_integral_v<IntType>, "Integral required");

    List<IntType> count
    (
        UPstream::listGatherValues<IntType>(localSize, communicator)
    );

    // Replace element [0] with the max receive count
    if (count.size() == 1)
    {
        // Nothing to receive from self
        count[0] = 0;
    }
    else if (count.size() > 1)
    {
        count[0] = *(std::max_element(count.begin()+1, count.end()));
    }

    return count;
}


template<class OffsetType, class IntType>
Foam::List
<
    std::conditional_t
    <std::is_void_v<OffsetType>, IntType, OffsetType>
>
Foam::globalIndexBase::calcOffsets
(
    const UList<IntType>& counts,
    bool checkOverflow
)
{
    using offset_type =
        std::conditional_t
        <std::is_void_v<OffsetType>, IntType, OffsetType>;

    // Could check this, but leave as undocumented requirement for now
    // static_assert(sizeof(offset_type) >= sizeof(IntType))
    // also assume that they both use consistent signed/unsigned

    List<offset_type> starts;

    if (const label len = counts.size(); len > 0)
    {
        starts.resize(len+1);

        offset_type total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = static_cast<offset_type>(counts[i]);

            // Ignore negative values as possible placeholders,
            // and 0 has no effect on total anyhow.
            if (counts[i] < 1)
            {
                starts[i] = total;
                continue;
            }

            // No overflow check for int64.
            if constexpr (sizeof(offset_type) < sizeof(int64_t))
            {
                if (checkOverflow && isOverflow(total, count))
                {
                    reportOverflowAndExit(i, total, count);
                }
            }

            starts[i] = total;
            total += count;
        }
        starts[len] = total;
    }

    return starts;
}


template<class OffsetType, class IntType, class Addr>
Foam::List
<
    std::conditional_t
    <std::is_void_v<OffsetType>, IntType, OffsetType>
>
Foam::globalIndexBase::calcOffsets
(
    const IndirectListBase<IntType, Addr>& counts,
    bool checkOverflow
)
{
    using offset_type =
        std::conditional_t
        <std::is_void_v<OffsetType>, IntType, OffsetType>;

    // Could check this, but leave as undocumented requirement for now
    // static_assert(sizeof(offset_type) >= sizeof(IntType))
    // also assume that they both use consistent signed/unsigned

    List<offset_type> starts;

    if (const label len = counts.size(); len > 0)
    {
        starts.resize(len+1);

        offset_type total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = static_cast<offset_type>(counts[i]);

            // Ignore negative values as possible placeholders,
            // and 0 has no effect on total anyhow.
            if (counts[i] < 1)
            {
                starts[i] = total;
                continue;
            }

            // No overflow check for int64.
            if constexpr (sizeof(offset_type) < sizeof(int64_t))
            {
                if (checkOverflow && isOverflow(total, count))
                {
                    reportOverflowAndExit(i, total, count);
                }
            }

            starts[i] = total;
            total += count;
        }
        starts[len] = total;
    }

    return starts;
}


template<class OffsetType, class SubListType>
Foam::List
<
    std::conditional_t
    <std::is_void_v<OffsetType>, Foam::label, OffsetType>
>
Foam::globalIndexBase::calcListOffsets
(
    const UList<SubListType>& lists,
    bool checkOverflow
)
{
    using offset_type =
        std::conditional_t
        <std::is_void_v<OffsetType>, Foam::label, OffsetType>;

    // Could check this, but leave as undocumented requirement for now
    // static_assert(sizeof(offset_type) >= sizeof(IntType))
    // also assume that they both use consistent signed/unsigned

    List<offset_type> starts;

    if (const label len = lists.size(); len > 0)
    {
        starts.resize(len+1);

        offset_type total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = static_cast<offset_type>(lists[i].size());

            // Zero-sized lists have no effect on total
            if (lists[i].size() < 1)
            {
                starts[i] = total;
                continue;
            }

            // No overflow check for int64.
            // Note: also checks endOffset of last range

            if constexpr (sizeof(offset_type) < sizeof(int64_t))
            {
                if (checkOverflow && isOverflow(total, count))
                {
                    reportOverflowAndExit(i, total, count);
                }
            }

            starts[i] = total;
            total += count;
        }
        starts[len] = total;
    }

    return starts;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        allData.resize_nocopy(offsets().back());
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
        allData.resize_nocopy(offsets().back());  // == totalSize()
    }
    else
    {
        allData.clear();  // zero-size on non-master
    }

    if (!UPstream::usingNodeComms(comm))
    {
        globalIndex::gather
        (
            offsets(),  // needed on master only
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
        allData.resize_nocopy(offsets().back());  // == totalSize()
    }
    else
    {
        allData.clear();  // zero-size on non-master
    }

    {
        globalIndex::gather
        (
            offsets(),  // needed on master only
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
        allData.resize_nocopy(offsets().back());  // == totalSize()
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

            // Calculate sizes. Currently without std::minus <functional>
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
            offsets(),  // needed on master only
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
            offsets(),  // needed on master only
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
