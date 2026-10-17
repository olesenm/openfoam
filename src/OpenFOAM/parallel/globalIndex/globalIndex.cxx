/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
    Copyright (C) 2026 Mark Olesen
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
#include "Pair.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::globalIndexBase::reportOverflowAndExit
(
    label idx,
    int64_t prevOffset,
    int64_t count
)
{
    if (idx < 0)
    {
        // No overflow tagged
        return;
    }

    FatalErrorInFunction
        << "Overflow : sum of sizes exceeds labelMax ("
        << labelMax << ") after index " << idx;

    if (prevOffset >= 0 && count >= 0)
    {
        FatalError
            << " while trying to add (" << count
            << ") to offset (" << prevOffset << ")";
    }

    FatalError
        << nl
        << "Please recompile with larger datatype for label." << nl
        << exit(FatalError);
}


Foam::labelList
Foam::globalIndex::calcOffsets
(
    const labelUList& counts,
    const bool checkOverflow
)
{
    labelList starts;

    if (const label len = counts.size(); len > 0)
    {
        starts.resize(len+1);

        label total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = counts[i];

            if (checkOverflow && isOverflow(total, count))
            {
                reportOverflowAndExit(i, total, count);
            }

            starts[i] = total;
            total += count;
        }
        starts[len] = total;
    }

    return starts;
}


Foam::List<Foam::labelRange>
Foam::globalIndex::calcRanges
(
    const labelUList& counts,
    const bool checkOverflow
)
{
    List<labelRange> values;

    if (const label len = counts.size(); len > 0)
    {
        values.resize(len);

        label total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = counts[i];

            // Note: do not care about endOffset of last range

            if (checkOverflow && (i < len-1) && isOverflow(total, count))
            {
                reportOverflowAndExit(i, total, count);
            }

            values[i].reset(total, count);
            total += count;
        }
    }

    return values;
}


bool Foam::globalIndex::splitNodeOffsets
(
    labelList& interNodeOffsets,
    labelList& localNodeOffsets,
    int communicator,
    const bool absoluteLocalNodeOffsets
) const
{
    // Require const-world as the starting point
    if (!UPstream::parRun() || communicator != UPstream::commConstWorld())
    {
        interNodeOffsets.clear();
        localNodeOffsets.clear();
        return false;
    }

    const auto interNodeComm = UPstream::commInterNode();

    // Only generate information on the node leaders
    if (!UPstream::is_rank(interNodeComm))
    {
        interNodeOffsets.clear();
        localNodeOffsets.clear();
        return true;  // Not involved, but return true to match others...
    }

    const label numProc = UPstream::nProcs(UPstream::commConstWorld());
    const auto& procIds = UPstream::procID(interNodeComm);
    const int ranki = UPstream::myProcNo(interNodeComm);

    if (FOAM_UNLIKELY(procIds.empty()))
    {
        // Should not happen...
        interNodeOffsets.clear();
        localNodeOffsets.clear();
        return true;  // Return true to match others...
    }

    // The inter-node offsets from the node-specific segment of the
    // overall offsets, but avoiding MPI_Scatterv (slow, doesn't
    // handle overlaps) and using MPI_Bcast() instead.

    // Send top-level offsets to the node leaders.
    // Could also be a mutable operation and use offsets_ directly.
    //
    // - number of overall offsets is always (nProc+1) [worldComm]
    labelList allOffsets;
    if (UPstream::master(interNodeComm))
    {
        allOffsets = offsets();
    }
    else  // ie, UPstream::is_subrank(interNodeComm)
    {
        allOffsets.resize_nocopy(numProc+1);
    }

    UPstream::broadcast(allOffsets.data(), allOffsets.size(), interNodeComm);


    if (FOAM_UNLIKELY(allOffsets.empty()))
    {
        // Should not happen...
        interNodeOffsets.clear();
        localNodeOffsets.clear();
        return true;  // Return true to match others...
    }

    // The local node span
    const label firstProc = procIds[ranki];
    const label lastProc =
    (
        (ranki+1 < procIds.size())
      ? procIds[ranki+1]
      : numProc
    );

    // Offsets (within a node)
    localNodeOffsets = allOffsets.slice
    (
        firstProc,
        (lastProc - firstProc) + 1  // +1 since offsets
    );

    if (!absoluteLocalNodeOffsets && !localNodeOffsets.empty())
    {
        const auto start0 = localNodeOffsets.front();
        for (auto& val : localNodeOffsets)
        {
            val -= start0;
        }
    }

    // Offsets (between nodes)
    interNodeOffsets.resize_nocopy(procIds.size()+1);
    {
        forAll(procIds, i)
        {
            interNodeOffsets[i] = allOffsets[procIds[i]];
        }
        interNodeOffsets.back() = allOffsets.back();
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndex::globalIndex(Istream& is)
{
    is >> offsets_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::CompactListList<Foam::label>
Foam::globalIndex::bin
(
    const labelUList& offsets,
    const labelUList& globalIds,
    labelList& order,
    DynamicList<label>& validBins
)
{
    Foam::sortedOrder(globalIds, order);
    validBins.clear();

    CompactListList<label> bins;

    if (!globalIds.empty())
    {
        labelList& binOffsets = bins.offsets();
        binOffsets.resize(offsets.size(), 0);

        labelList& binValues = bins.values();
        binValues = UIndirectList<label>(globalIds, order);

        const label id = binValues[0];
        label proci = findLower(offsets, id+1);

        validBins.push_back(proci);
        label binSize = 1;

        for (label i = 1; i < order.size(); i++)
        {
            const label id = binValues[i];

            if (id < offsets[proci+1])
            {
                ++binSize;
            }
            else
            {
                // Not local. Reset proci
                const label oldProci = proci;
                proci = findLower(offsets, id+1);

                // Set offsets
                for (label j = oldProci+1; j < proci; ++j)
                {
                    binOffsets[j] = binOffsets[oldProci]+binSize;
                }
                binOffsets[proci] = i;
                validBins.push_back(proci);
                binSize = 1;
            }
        }

        for (label j = proci+1; j < binOffsets.size(); ++j)
        {
            binOffsets[j] = binOffsets[proci]+binSize;
        }
    }

    return bins;
}


void Foam::globalIndex::resize(label n)
{
    if (n < 1)
    {
        offsets_.clear();
    }
    else
    {
        offsets_.resize(n+1, end_value());
    }
}


void Foam::globalIndex::reset_nonParallel
(
    label localSize,
    int communicator
)
{
    // localSize on-proc, zero elsewhere:
    // offsets are 0 up to and including proci, then localSize
    if
    (
        label len = 0;
        (communicator >= 0 && (len = UPstream::nProcs(communicator)) > 0)
    )
    {
        offsets_.resize_fill(len+1, 0);
        if (label proci = UPstream::myProcNo(communicator); proci >= 0)
        {
            for (label i = proci+1; i <= len; ++i)
            {
                offsets_[i] = localSize;
            }
        }
    }
    else
    {
        offsets_.clear();
    }
}


void Foam::globalIndex::reset
(
    const label localSize,
    int communicator,
    bool parallel
)
{
    if (parallel && UPstream::is_parallel(communicator))
    {
        constexpr bool checkOverflow = true;

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


void Foam::globalIndex::reset
(
    const labelUList& counts,
    const bool checkOverflow
)
{
    if (const label len = counts.size(); len > 0)
    {
        offsets_.resize_nocopy(len+1);

        label total = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];

            if (checkOverflow && isOverflow(total, count))
            {
                reportOverflowAndExit(i, total, count);
            }

            offsets_[i] = total;
            total += count;
        }
        offsets_[len] = total;
    }
    else
    {
        offsets_.clear();
    }
}


void Foam::globalIndex::reset
(
    globalIndex::gatherOffset,
    label localSize,
    int communicator,
    bool parallel
)
{
    if (!parallel)
    {
        single(localSize);
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

        // Work array : allocated with double space! (still 0 on non-master)
        labelList work(2*numCounts);

        // Slice into the initial part of the work array
        // for getting/setting the counts.

        labelUList counts(work.data(), numCounts);


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

            const auto total = totalSize();

            auto iter = offsets().cbegin();

            for (label proci = 0; proci < numProc; ++proci)
            {
                work[2*proci]   = *iter++;
                work[2*proci+1] = total;
            }
        }

        label dims[2];  // (start, total) tuple

        UPstream::mpiScatter(work.cdata(), dims, 2, communicator);

        if (UPstream::is_subrank(communicator))
        {
            // Store either as (beg len total) tuple or if there is no
            // no trailing content simply as a (beg len) tuple
            single(dims[0], localSize, dims[1]);
        }
    }
    else
    {
        reset_nonParallel(localSize, communicator);
    }
}


void Foam::globalIndex::setLocalSize(const label proci, const label len)
{
    if (proci >= 0 && proci+1 < offsets_.size() && len >= 0)
    {
        const label delta = (len - (offsets_[proci+1] - offsets_[proci]));

        // TBD: additional overflow check
        if (delta)
        {
            for (label i = proci+1; i < offsets_.size(); ++i)
            {
                offsets_[i] += delta;
            }
        }
    }
}


Foam::labelList Foam::globalIndex::localSizes() const
{
    labelList values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label proci=0; proci < len; ++proci)
    {
        values[proci] = offsets_[proci+1] - offsets_[proci];
    }

    return values;
}


Foam::List<Foam::labelRange>
Foam::globalIndex::ranges() const
{
    List<labelRange> values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label proci=0; proci < len; ++proci)
    {
        values[proci].reset
        (
            offsets_[proci],
            (offsets_[proci+1] - offsets_[proci])
        );
    }

    return values;
}


Foam::label Foam::globalIndex::maxNonLocalSize(const label proci) const
{
    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return 0;
    }

    label maxLen = 0;

    for (label i=0; i < len; ++i)
    {
        if (i != proci)
        {
            const label count = (offsets_[i+1] - offsets_[i]);
            maxLen = std::max(maxLen, count);
        }
    }

    return maxLen;
}


Foam::labelRange Foam::globalIndex::front() const
{
    return
    (
        (offsets_.size() < 2)
      ? labelRange()
      : labelRange(Pair<label>(offsets_[0], offsets_[1]))
    );
}


Foam::labelRange Foam::globalIndex::back() const
{
    return
    (
        (offsets_.size() < 2)
      ? labelRange()
      : labelRange
        (
            Pair<label>
            (
                offsets_[offsets_.size()-2],
                offsets_[offsets_.size()-1]
            )
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, globalIndex& gi)
{
    return is >> gi.offsets();
}


Foam::Ostream& Foam::operator<<(Ostream& os, const globalIndex& gi)
{
    return os << gi.offsets();
}


// ************************************************************************* //
