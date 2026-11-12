/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "simpleGeomDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "globalIndex.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleGeomDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        simpleGeomDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// assignToProcessorGroup : given nCells cells and nProcGroup processor
// groups to share them, how do we share them out? Answer : each group
// gets nCells/nProcGroup cells, and the first few get one
// extra to make up the numbers. This should produce almost
// perfect load balancing

namespace Foam
{

static void assignToProcessorGroup
(
    labelUList& processorGroup,
    const label nProcGroup
)
{
    const label jump = processorGroup.size()/nProcGroup;
    const label jumpb = jump + 1;
    const label fstProcessorGroup = processorGroup.size() - jump*nProcGroup;

    label ind = 0;
    label j = 0;

    // assign cells to the first few processor groups (those with
    // one extra cell each
    for (j=0; j<fstProcessorGroup; j++)
    {
        for (label k=0; k<jumpb; k++)
        {
            processorGroup[ind++] = j;
        }
    }

    // and now to the `normal' processor groups
    for (; j<nProcGroup; j++)
    {
        for (label k=0; k<jump; k++)
        {
            processorGroup[ind++] = j;
        }
    }
}


static void assignToProcessorGroup
(
    labelUList& processorGroup,
    const label nProcGroup,
    const labelUList& indices,
    const UList<scalar>& weights,
    const scalar summedWeights
)
{
    // This routine gets the sorted points.
    // Easiest to explain with an example.
    // E.g. 400 points, sum of weights : 513.
    // Now with number of divisions in this direction (nProcGroup) : 4
    // gives the split at 513/4 = 128
    // So summed weight from 0..128 goes into bin 0,
    //     ,,              128..256 goes into bin 1
    //   etc.
    // Finally any remaining ones go into the last bin (3).

    const scalar jump = summedWeights/nProcGroup;
    const label nProcGroupM1 = nProcGroup - 1;

    scalar sumWeights = 0;
    label ind = 0;
    label j = 0;

    // assign cells to all except last group.
    for (j=0; j<nProcGroupM1; j++)
    {
        const scalar limit = jump*scalar(j + 1);
        while (sumWeights < limit)
        {
            sumWeights += weights[indices[ind]];
            processorGroup[ind++] = j;
        }
    }
    // Ensure last included.
    processorGroup.slice(ind) = nProcGroupM1;
}


// As above, but split the band so that bin j receives the target fraction
// fractions[j] of the cells by count. Used when processorWeights requests an
// uneven banding.
static void assignToProcessorGroup
(
    labelUList& processorGroup,
    const label nProcGroup,
    const UList<scalar>& fractions  // per-band target fractions (sum 1)
)
{
    const label n = processorGroup.size();
    const label nProcGroupM1 = nProcGroup - 1;  // last bin ID

    label ind = 0;
    scalar cumFrac = 0;

    // Visit every bin except the last; the last bin receives the remainder
    for (label j = 0; j < nProcGroupM1; ++j)
    {
        cumFrac += fractions[j];
        const label limit = Foam::min(n, label(cumFrac*scalar(n) + 0.5));
        while (ind < limit)
        {
            processorGroup[ind++] = j;
        }
    }
    processorGroup.slice(ind) = nProcGroupM1;
}


// As above (cell-weighted), but split so that bin j receives the target
// fraction fractions[j] of the total cell weight.
static void assignToProcessorGroup
(
    labelUList& processorGroup,
    const label nProcGroup,
    const labelUList& indices,
    const UList<scalar>& weights,
    const scalar summedWeights,
    const UList<scalar>& fractions  // per-band target fractions (sum 1)
)
{
    const label n = processorGroup.size();
    const label nProcGroupM1 = nProcGroup - 1;  // last bin ID

    scalar sumWeights = 0;
    scalar cumFrac = 0;
    label ind = 0;

    // Visit all bins except the last
    for (label j = 0; j < nProcGroupM1; ++j)
    {
        cumFrac += fractions[j];
        const scalar limit = cumFrac*summedWeights;
        while (ind < n && sumWeights < limit)
        {
            sumWeights += weights[indices[ind]];
            processorGroup[ind++] = j;
        }
    }
    processorGroup.slice(ind) = nProcGroupM1;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::simpleGeomDecomp::decomposeOneProc
(
    const pointField& points
) const
{
    // construct a list for the final result
    labelList finalDecomp(points.size());

    labelList processorGroups(points.size());

    labelList pointIndices(Foam::identity(points.size()));

    const pointField rotatedPoints(adjustPoints(points));

    // For each direction (X, then Y, then Z) assign the processors to
    // groups of processors labelled 0..nX (or nY, nZ) to give a banded
    // structure on the mesh. The actual processor number is constructed
    // by treating the X group as the units, the Y group as a multiple of
    // nX (columns) and the Z group as a multiple of nX*nY (layers).

    label processed = 1;

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const label count = n_[cmpt];

        Foam::sort
        (
            pointIndices,
            [&rotatedPoints,cmpt](label a, label b)
            {
                return rotatedPoints[a][cmpt] < rotatedPoints[b][cmpt];
            }
        );

        if (procWeights_.empty())
        {
            assignToProcessorGroup(processorGroups, count);
        }
        else
        {
            assignToProcessorGroup
            (
                processorGroups,
                count,
                axisWeightFractions(cmpt)
            );
        }

        if (cmpt == 0)
        {
            // Assign on the first pass
            forAll(points, i)
            {
                finalDecomp[pointIndices[i]] = processorGroups[i];
            }
        }
        else
        {
            // Accumulate on subsequent passes
            forAll(points, i)
            {
                finalDecomp[pointIndices[i]] += processed*processorGroups[i];
            }
        }

        processed *= count;  // next
    }

    return finalDecomp;
}


Foam::labelList Foam::simpleGeomDecomp::decomposeOneProc
(
    const pointField& points,
    const UList<scalar>& weights
) const
{
    double totalWeight = 0;

    bool hadNegative = false;

    for (const auto& val : weights)
    {
        if (val < 0)
        {
            hadNegative = true;
        }
        else
        {
            totalWeight += double(val);
        }
    }

    bool unweighted = (weights.empty() || hadNegative);

    if (hadNegative)
    {
        WarningInFunction
            << "Negative weights, ignoring weights" << nl;
    }
    else if (!weights.empty())
    {
        const double magTotal = std::fabs(totalWeight);

        if (magTotal < double(VSMALL))
        {
            // Near zero weights for this scalar (float/double) type.
            // Treat as unweighted without a message,
            // it could be a legitimate input for unweighted

            unweighted = true;
        }
        else if (double(VGREAT) < magTotal)
        {
            // Excessive total sum for scalar (float/double) type.
            // This type of overflow will almost never occur
            // since the regular inputs will be "normal" values
            // like 0,1, 1000, etc. Never really anything silly like 1E+9...

            // Report this rare case, and treat as unweighted so that
            // things will work (not crash).

            unweighted = true;
            WarningInFunction
                << "Sum of weights exceeds scalar range, ignoring weights"
                << nl;
        }
    }

    if (unweighted)
    {
        return decomposeOneProc(points);
    }

    // Already checked that totalWeight was in range, non-negative etc so
    // can just use it now.
    const scalar summedWeights = scalar(totalWeight);

    // construct a list for the final result
    labelList finalDecomp(points.size());

    labelList processorGroups(points.size());

    labelList pointIndices(Foam::identity(points.size()));

    const pointField rotatedPoints(adjustPoints(points));

    // For each direction (X, then Y, then Z) assign the processors to
    // groups of processors labelled 0..nX (or nY, nZ) to give a banded
    // structure on the mesh. The actual processor number is constructed
    // by treating the X group as the units, the Y group as a multiple of
    // nX (columns) and the Z group as a multiple of nX*nY (layers).

    label processed = 1;

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const label count = n_[cmpt];

        Foam::sort
        (
            pointIndices,
            [&rotatedPoints,cmpt](label a, label b)
            {
                return rotatedPoints[a][cmpt] < rotatedPoints[b][cmpt];
            }
        );

        if (procWeights_.empty())
        {
            assignToProcessorGroup
            (
                processorGroups,
                count,
                pointIndices,
                weights,
                summedWeights
            );
        }
        else
        {
            assignToProcessorGroup
            (
                processorGroups,
                count,
                pointIndices,
                weights,
                summedWeights,
                axisWeightFractions(cmpt)
            );
        }

        if (cmpt == 0)
        {
            // Assign on the first pass
            forAll(points, i)
            {
                finalDecomp[pointIndices[i]] = processorGroups[i];
            }
        }
        else
        {
            // Accumulate on subsequent passes
            forAll(points, i)
            {
                finalDecomp[pointIndices[i]] += processed*processorGroups[i];
            }
        }

        processed *= count;  // next
    }

    return finalDecomp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleGeomDecomp::simpleGeomDecomp(const Vector<label>& divisions)
:
    geomDecomp(divisions)
{}


Foam::simpleGeomDecomp::simpleGeomDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    geomDecomp(typeName, decompDict, regionName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::simpleGeomDecomp::decompose
(
    const pointField& points,
    const scalarField& weights
) const
{
    // May want to handle different communicators in the future...
    const int communicator = UPstream::worldComm;
    const int messageTag = (UPstream::msgType() + 1812);
    const auto commsType = globalIndex::defaultCommsType();

    // Uniform weighting if weights are empty or poorly sized
    const bool isWeighted = returnReduceAnd
    (
        (points.size() == weights.size()),
        communicator
    );


    labelList finalDecomp;

    if (UPstream::is_parallel(communicator))
    {
        // Parallel
        const globalIndex globalNumbers(points.size(), communicator);

        if
        (
            globalNumbers.empty()
         || (globalNumbers.localSize(0) == globalNumbers.totalSize())
        )
        {
            // Only have values on the master (ie, decomposing)
            if (UPstream::master(communicator))
            {
                if (isWeighted)
                {
                    finalDecomp = decomposeOneProc(points, weights);
                }
                else
                {
                    finalDecomp = decomposeOneProc(points);
                }
            }
        }
        else
        {
            // Gather all, decompose on master, scatter back

            pointField allPoints
            (
                globalNumbers.gather
                (
                    points,
                    messageTag, commsType, communicator
                )
            );

            labelList allDecomp;
            if (isWeighted)
            {
                // Non-uniform weighting
                scalarField allWeights
                (
                    globalNumbers.gather
                    (
                        weights,
                        messageTag, commsType, communicator
                    )
                );

                if (UPstream::master(communicator))
                {
                    allDecomp = decomposeOneProc(allPoints, allWeights);
                }
            }
            else
            {
                if (UPstream::master(communicator))
                {
                    allDecomp = decomposeOneProc(allPoints);
                }
            }

            allPoints.clear();

            finalDecomp = globalNumbers.scatter
            (
                allDecomp,
                messageTag, commsType, communicator
            );
        }
    }
    else
    {
        // Non-parallel (or single-proc multi-world)

        if (isWeighted)
        {
            finalDecomp = decomposeOneProc(points, weights);
        }
        else
        {
            finalDecomp = decomposeOneProc(points);
        }
    }

    return finalDecomp;
}


// ************************************************************************* //
