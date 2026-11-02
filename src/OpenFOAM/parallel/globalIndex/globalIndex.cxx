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
        << "Overflow : sum of sizes exceeds integral limit after index " << idx;

    if (prevOffset >= 0 && count >= 0)
    {
        FatalError
            << " while trying to add (" << count
            << ") to offset (" << prevOffset << ")";
    }

    FatalError
        << nl
        << "Please consider a larger integral datatype." << nl
        << Foam::exit(FatalError);
}


Foam::List<Foam::labelRange>
Foam::globalIndex::calcRanges
(
    const labelUList& counts,
    bool checkOverflow
)
{
    List<labelRange> values;

    if (const label len = counts.size(); len > 0)
    {
        values.resize(len);

        label total = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];

            // No overflow check for int64...
            // Note: do not care about endOffset of last range
            // if constexpr (sizeof(IntType) < sizeof(int64_t))
            {
                if
                (
                    checkOverflow && (i < len-1)
                 && globalIndexBase::isOverflow(total, count)
                )
                {
                    globalIndexBase::reportOverflowAndExit(i, total, count);
                }
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


Foam::List<Foam::labelRange>
Foam::globalIndex::ranges() const
{
    // CSRAddressing delivers IntRange<label>, we want labelRange

    List<labelRange> values;

    const label len = (offsets().size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label i = 0; i < len; ++i)
    {
        values[i].reset
        (
            offsets()[i],
            (offsets()[i+1] - offsets()[i])
        );
    }

    return values;
}


// ************************************************************************* //
