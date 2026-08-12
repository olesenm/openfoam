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

namespace Foam
{

// A list compare binary predicate for normal sort by vector component
struct vectorLessOp
{
    const UList<vector>& values;
    direction sortCmpt;

    vectorLessOp(const UList<vector>& list, direction cmpt = vector::X)
    :
        values(list),
        sortCmpt(cmpt)
    {}

    void setComponent(direction cmpt)
    {
        sortCmpt = cmpt;
    }

    bool operator()(const label a, const label b) const
    {
        return values[a][sortCmpt] < values[b][sortCmpt];
    }
};

} // End namespace Foam


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
    const scalarField& weights,
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
    const scalarUList& fractions    // per-band target fractions (sum 1)
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
    const scalarField& weights,
    const scalar summedWeights,
    const scalarUList& fractions    // per-band target fractions (sum 1)
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

    labelList pointIndices(identity(points.size()));

    const pointField rotatedPoints(adjustPoints(points));

    vectorLessOp sorter(rotatedPoints);

    // and one to take the processor group id's. For each direction.
    // we assign the processors to groups of processors labelled
    // 0..nX to give a banded structure on the mesh. Then we
    // construct the actual processor number by treating this as
    // the units part of the processor number.

    sorter.setComponent(vector::X);
    Foam::sort(pointIndices, sorter);

    if (procWeights_.empty())
    {
        assignToProcessorGroup(processorGroups, n_.x());
    }
    else
    {
        assignToProcessorGroup
        (
            processorGroups,
            n_.x(),
            axisWeightFractions(vector::X)
        );
    }

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] = processorGroups[i];
    }


    // now do the same thing in the Y direction. These processor group
    // numbers add multiples of nX to the proc. number (columns)

    sorter.setComponent(vector::Y);
    Foam::sort(pointIndices, sorter);

    if (procWeights_.empty())
    {
        assignToProcessorGroup(processorGroups, n_.y());
    }
    else
    {
        assignToProcessorGroup
        (
            processorGroups,
            n_.y(),
            axisWeightFractions(vector::Y)
        );
    }

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*processorGroups[i];
    }


    // finally in the Z direction. Now we add multiples of nX*nY to give
    // layers

    sorter.setComponent(vector::Z);
    Foam::sort(pointIndices, sorter);

    if (procWeights_.empty())
    {
        assignToProcessorGroup(processorGroups, n_.z());
    }
    else
    {
        assignToProcessorGroup
        (
            processorGroups,
            n_.z(),
            axisWeightFractions(vector::Z)
        );
    }

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*n_.y()*processorGroups[i];
    }

    return finalDecomp;
}


Foam::labelList Foam::simpleGeomDecomp::decomposeOneProc
(
    const pointField& points,
    const scalarField& weights
) const
{
    if (weights.empty())
    {
        return decomposeOneProc(points);
    }

    // construct a list for the final result
    labelList finalDecomp(points.size());

    labelList processorGroups(points.size());

    labelList pointIndices(identity(points.size()));

    const pointField rotatedPoints(adjustPoints(points));

    vectorLessOp sorter(rotatedPoints);

    // and one to take the processor group id's. For each direction.
    // we assign the processors to groups of processors labelled
    // 0..nX to give a banded structure on the mesh. Then we
    // construct the actual processor number by treating this as
    // the units part of the processor number.

    sorter.setComponent(vector::X);
    Foam::sort(pointIndices, sorter);

    const scalar summedWeights = sum(weights);
    if (procWeights_.empty())
    {
        assignToProcessorGroup
        (
            processorGroups,
            n_.x(),
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
            n_.x(),
            pointIndices,
            weights,
            summedWeights,
            axisWeightFractions(vector::X)
        );
    }

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] = processorGroups[i];
    }


    // now do the same thing in the Y direction. These processor group
    // numbers add multiples of nX to the proc. number (columns)

    sorter.setComponent(vector::Y);
    Foam::sort(pointIndices, sorter);

    if (procWeights_.empty())
    {
        assignToProcessorGroup
        (
            processorGroups,
            n_.y(),
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
            n_.y(),
            pointIndices,
            weights,
            summedWeights,
            axisWeightFractions(vector::Y)
        );
    }

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*processorGroups[i];
    }


    // finally in the Z direction. Now we add multiples of nX*nY to give
    // layers

    sorter.setComponent(vector::Z);
    Foam::sort(pointIndices, sorter);

    if (procWeights_.empty())
    {
        assignToProcessorGroup
        (
            processorGroups,
            n_.z(),
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
            n_.z(),
            pointIndices,
            weights,
            summedWeights,
            axisWeightFractions(vector::Z)
        );
    }

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*n_.y()*processorGroups[i];
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
    // Uniform weighting if weights are empty or poorly sized
    const bool hasWeights = returnReduceAnd(points.size() == weights.size());

    if (!UPstream::parRun())
    {
        if (hasWeights)
        {
            return decomposeOneProc(points, weights);
        }
        else
        {
            return decomposeOneProc(points);
        }
    }
    else
    {
        // Parallel
        const globalIndex globalNumbers(points.size());

        labelList allDecomp;
        pointField allPoints(globalNumbers.gather(points));
        scalarField allWeights;

        if (hasWeights)
        {
            // Non-uniform weighting
            allWeights = globalNumbers.gather(weights);
        }

        if (UPstream::master())
        {
            if (hasWeights)
            {
                allDecomp = decomposeOneProc(allPoints, allWeights);
            }
            else
            {
                allDecomp = decomposeOneProc(allPoints);
            }
            allPoints.clear();  // Not needed anymore
            allWeights.clear(); // ...
        }

        return globalNumbers.scatter(allDecomp);
    }
}


// ************************************************************************* //
