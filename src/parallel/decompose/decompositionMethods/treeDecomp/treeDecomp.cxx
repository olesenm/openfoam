/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Notes
    The split index is chosen to so that each side of the split receive point
    counts proportional to its number of target domains. For example,
    splitting into 9 domains gives 4/9 of the points to the front bin and
    5/9 to the back bin.

    When processor weights are specified, the bin split and the point-split
    fraction are instead based on the relative weight of the
    domains on each side, rather than on domain counts. A domain with a
    weight of zero is deliberately starved of points, even when points
    remain plentiful.

    Recursion continues, alternating the split direction as dictated by the
    local point spread, until every domain corresponds to a group of
    points. Domains are only left empty when there are fewer points than
    requested domains, or (with processorWeights) when a domain has been
    assigned a weight of zero; otherwise every domain receives at least one
    point.

\*---------------------------------------------------------------------------*/

#include "treeDecomp.H"
#include "treeBoundBox.H"
#include "bitSet.H"
#include "globalIndex.H"
#include "addToRunTimeSelectionTable.H"
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(treeDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        treeDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Copy integral weights as scalar weights, treat -ve values like 0.
// Leaves weights unnormalized (except in fulldebug)
static void copy_weights
(
    const label numDomains,
    const labelUList& weights,
    List<scalar>& procWeights
)
{
    if (numDomains <= 1 || weights.empty())
    {
        procWeights.clear();
        return;
    }

    // Silently adjust range
    const label overlap = Foam::min(numDomains, weights.size());

    bool uniform = false;
    label nWeights = 0;
    label firstValue = -1;
    double totalWeight = 0;

    for (label i = 0; i < overlap; ++i)
    {
        if (label val = weights[i]; val > 0)
        {
            ++nWeights;
            totalWeight += val;

            if (firstValue < 0)
            {
                firstValue = val;
                uniform = true;
            }
            else if (uniform && (firstValue != val))
            {
                uniform = false;
            }
        }
    }

    if
    (
        (nWeights == 0)
     || (uniform && (nWeights == numDomains))
     || (totalWeight <= doubleScalarSMALL)
    )
    {
        // Identical weights everywhere,
        // or degenerate case of zero weights everywhere
        procWeights.clear();
    }
    else
    {
        // Appears to have valid and non-uniform weights
        // - can leave in non-normalized state but copy label -> scalar

        procWeights.resize_fill(numDomains, 0);

        for (label i = 0; i < overlap; ++i)
        {
            if (label val = weights[i]; val > 0)
            {
                procWeights[i] = val;
                #ifdef FULLDEBUG
                // Normalized value (for debugging)
                procWeights[i] /= totalWeight;
                #endif
            }
        }
    }
}


// Processor weights based on boolean selection.
// Leaves weights unnormalized (always)
template<class BoolListType>
static void select_weights
(
    const label numDomains,
    const BoolListType& selection,
    List<scalar>& procWeights
)
{
    if (numDomains <= 1)
    {
        procWeights.clear();
        return;
    }

    procWeights.resize_fill(numDomains, 0);

    // Can walk the full numDomains since .test() is failsafe for out-of-range

    label nWeights = 0;
    for (label i = 0; i < numDomains; ++i)
    {
        if (selection.test(i))
        {
            ++nWeights;
            procWeights[i] = 1;
        }
    }

    if ((nWeights == 0) || (nWeights == numDomains))
    {
        // Identical weights everywhere,
        // or degenerate case of zero weights everywhere
        procWeights.clear();
    }
    else
    {
        // Appears to have valid and non-uniform weights
        // - can leave in non-normalized state
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::treeDecomp::splitPoints
(
    const UList<point>& points,
    labelUList& indices,
    const labelRange& span,
    const labelRange& bins,
    const UList<scalar>& procWeights,
    labelUList& finalDecomp
)
{
    const bool hasProcWeights = !procWeights.empty();

    // Nothing to split - either a single target domain, or no more points
    // to distribute.
    //
    // Assign what remains to a domain of this sub-tree. Normally this is
    // simply the first domain, but if multiple domains remain unresolved
    // (points exhausted before the domain range collapsed to a single entry)
    // and weights are in use, prefer the domain with the largest weight.

    if (bins.size() <= 1 || span.size() <= 1)
    {
        label target = bins.begin_value();

        if (bins.size() > 1 && hasProcWeights)
        {
            scalar maxWeight = procWeights[target];

            for (label domain : bins)
            {
                if (maxWeight < procWeights[domain])
                {
                    maxWeight = procWeights[domain];
                    target = domain;
                }
            }
        }

        for (label i : span)
        {
            finalDecomp[indices[i]] = target;
        }
        return;
    }

    // Total weight assigned to the domains in this bins range.
    // Only meaningful (and used) when procWeights are supplied
    // and have at least one non-zero entry.

    double totalWeight = 0;

    if (hasProcWeights)
    {
        for (label domain : bins)
        {
            totalWeight += procWeights[domain];
        }
    }

    const bool useProcWeights =
    (
        hasProcWeights && (totalWeight > doubleScalarVSMALL)
    );

    // Determine the coordinate direction with the largest span
    // over the current group of points

    direction cmpt(vector::X);
    {
        auto subIndices = indices.slice(span);
        cmpt = treeBoundBox(points, subIndices).maxDir();
    }

    // Split the target domains as evenly as possible.
    // For an odd number of bins (eg, 9), this gives 4 to the front
    // half and 5 to the back half.
    label nLower = bins.size()/2;

    if (useProcWeights)
    {
        // Split the domains at the (approximate) weighted median:
        // walk the cumulative weight until at least half of the total
        // has been accounted for.

        const double half = 0.5*totalWeight;
        double cumulative = 0;

        nLower = 0;
        for (label domain : bins)
        {
            cumulative += procWeights[domain];
            ++nLower;

            if (half <= cumulative)
            {
                break;
            }
        }

        // Ensure the bin split makes progress on both sides
        nLower = Foam::clamp(nLower, label(1), bins.size()-1);
    }

    const labelRange binsLower
    (
        bins.begin_value(),
        nLower
    );
    const labelRange binsUpper
    (
        binsLower.end_value(),
        (bins.size() - binsLower.size())
    );

    // Choose the split index so that each half receives (as closely as
    // possible) the same fraction of points as the fraction of domains
    // (or, when weighted, the fraction of summed weight) it has been
    // assigned.

    double fraction =
    (
        double(binsLower.size())/double(bins.size())
    );

    // Pivot bounds: by default, both sides get at least one point
    // (progress guarantee). Only relaxed below for a weighted side
    // with (effectively) zero summed weight, which is deliberately
    // starved of points instead.
    label minPivot = (span.begin_value() + 1);
    label maxPivot = (span.end_value() - 1);

    if (useProcWeights)
    {
        double weightLower = 0;
        for (label domain : binsLower)
        {
            weightLower += procWeights[domain];
        }

        fraction = double(weightLower)/double(totalWeight);

        if (weightLower <= doubleScalarVSMALL)
        {
            minPivot = span.begin_value();
        }
        if ((totalWeight - weightLower) <= doubleScalarVSMALL)
        {
            maxPivot = span.end_value();
        }
    }

    label pivot =
    (
        span.begin_value()
      + label(std::lround(double(span.size())*fraction))
    );

    pivot = Foam::clamp(pivot, minPivot, maxPivot);

    // Partition (not fully sort) the points about the pivot location
    // along the chosen direction.
    std::nth_element
    (
        indices.begin() + span.begin_value(),
        indices.begin() + pivot,
        indices.begin() + span.end_value(),

        [&points,cmpt](label a, label b)
        {
            return (points[a][cmpt] < points[b][cmpt]);
        }
    );

    // The lower/upper point-index spans,
    // using the same pivot location determined above
    const labelRange spanLower
    (
        span.begin_value(),
        (pivot - span.begin_value())
    );
    const labelRange spanUpper
    (
        pivot,
        (span.end_value() - pivot)
    );

    // Recurse into the front/back halves, using the same proportional
    // split of the (target) number of domains (or summed weight)
    splitPoints
    (
        points,
        indices,
        spanLower,
        binsLower,
        procWeights,
        finalDecomp
    );
    splitPoints
    (
        points,
        indices,
        spanUpper,
        binsUpper,
        procWeights,
        finalDecomp
    );
}


Foam::labelList Foam::treeDecomp::treeDecompose_impl
(
    const label numDomains,
    const UList<point>& points,
    const UList<scalar>& procWeights
)
{
    const auto numPoints = points.size();

    // The final result
    labelList finalDecomp(numPoints, Foam::zero{});

    if (numDomains > 1)
    {
        // The sort indices
        labelList indices(Foam::identity(numPoints));

        const labelRange span(numPoints);   // span of point indices
        const labelRange bins(numDomains);  // range of target domain ids

        splitPoints(points, indices, span, bins, procWeights, finalDecomp);
    }

    return finalDecomp;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::treeDecomp::readCoeffs(const dictionary& coeffs)
{
    // Optional (cartesian) coordinate system transform
    auto csysPtr = coordinateSystem::NewIfPresent(coeffs, "transform");

    if (csysPtr)
    {
        csys_ = *csysPtr;
    }
    else
    {
        csys_.clear();  // Reset to identity
    }

    if
    (
        coeffs.readIfPresent("processorWeights", procWeights_)
     && (procWeights_.size() > 1)
    )
    {
        if (procWeights_.size() < nDomains_)
        {
            FatalIOErrorInFunction(coeffs)
                << "processorWeights (" << procWeights_.size()
                << ") less than numberOfSubdomains (" << nDomains_
                << ")\n"
                << exit(FatalIOError);
        }
        else if (procWeights_.size() > nDomains_)
        {
            WarningInFunction
                << "processorWeights (" << procWeights_.size()
                << ") greater than numberOfSubdomains (" << nDomains_
                << ") ... truncating" << endl;

            procWeights_.resize(nDomains_, 0);
        }

        sanitizeWeights(nDomains_, procWeights_);
    }
    else
    {
        // Also treat size 1 like uniform weighted
        procWeights_.clear();
    }
}


Foam::labelList Foam::treeDecomp::treeDecompose
(
    const pointField& points
) const
{
    // May want to handle different communicators in the future...
    const int communicator = UPstream::worldComm;
    const int messageTag = (UPstream::msgType() + 1812);
    const auto commsType = globalIndex::defaultCommsType();


    // The points, after any geometric transformations
    tmp<pointField> trotatedPoints(points);

    if (csys_.good() && !csys_.R().is_identity())
    {
        trotatedPoints = csys_.localPosition(points);
    }
    const auto& rotatedPoints = trotatedPoints.cref();

    const label numDomains = nDomains_;
    const auto& weights = procWeights_;


    labelList finalDecomp;

    if (UPstream::is_parallel(communicator))
    {
        // Parallel
        const globalIndex globalNumbers(rotatedPoints.size(), communicator);

        if
        (
            globalNumbers.empty()
         || (globalNumbers.localSize(0) == globalNumbers.totalSize())
        )
        {
            // Only have values on the master (ie, decomposing)
            if (UPstream::master(communicator))
            {
                finalDecomp =
                    treeDecompose_impl(numDomains, rotatedPoints, weights);
            }
        }
        else
        {
            // Gather all, decompose on master, scatter back

            pointField allPoints
            (
                globalNumbers.gather
                (
                    rotatedPoints,
                    messageTag, commsType, communicator
                )
            );

            labelList allDecomp;
            if (UPstream::master(communicator))
            {
                allDecomp = treeDecompose_impl(numDomains, allPoints, weights);
                allPoints.clear();
            }

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

        finalDecomp = treeDecompose_impl(numDomains, rotatedPoints, weights);
    }

    return finalDecomp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDecomp::treeDecomp(label numDomains)
:
    decompositionMethod(numDomains)
{}


Foam::treeDecomp::treeDecomp(label numDomains, const labelUList& weights)
:
    decompositionMethod(numDomains)
{
    copyWeights(numDomains, weights, procWeights_);
}


Foam::treeDecomp::treeDecomp(label numDomains, const bitSet& selection)
:
    decompositionMethod(numDomains)
{
    selectWeights(numDomains, selection, procWeights_);
}


Foam::treeDecomp::treeDecomp(label numDomains, const boolUList& selection)
:
    decompositionMethod(numDomains)
{
    selectWeights(numDomains, selection, procWeights_);
}


Foam::treeDecomp::treeDecomp(label numDomains, const coordinateSystem& csys)
:
    decompositionMethod(numDomains),
    csys_(csys)   // Check good() and !is_identity() later...
{}


Foam::treeDecomp::treeDecomp
(
    const dictionary& decompDict,
    const word& regionName,
    int select
)
:
    decompositionMethod(decompDict, regionName),
    csys_()
{
    const dictionary& coeffs = findCoeffsDict(typeName + "Coeffs", select);

    readCoeffs(coeffs);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::treeDecomp::setWeights(const labelUList& weights)
{
    procWeights_.clear();
    copy_weights(nDomains_, weights, procWeights_);
}


void Foam::treeDecomp::setWeights(const bitSet& selection)
{
    procWeights_.clear();
    selectWeights(nDomains_, selection, procWeights_);
}


void Foam::treeDecomp::setWeights(const boolUList& selection)
{
    procWeights_.clear();
    selectWeights(nDomains_, selection, procWeights_);
}


// ************************************************************************* //
