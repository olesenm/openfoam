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

Note
    When building the tree indexing (partitioning), a strict depth-first
    order is used which means that the addressing for the resulting leaves
    can simply be appended to a flat list of offsets.

\*---------------------------------------------------------------------------*/

#include "AABBTree.H"
#include <algorithm>
#include <cmath>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::AABBTreeBase::autosize
(
    control& ctrl,
    label total,
    int maxLevel,
    int minLeafSize
)
{
    // Absolute lower depth limit
    constexpr int minDepth = 1;

    // Small safety margin
    constexpr int depthSafety = 1;

    if (total <= 0)
    {
        return;
    }

    const bool autoDepth = (maxLevel <= 0);
    const bool autoWidth = (minLeafSize <= 0);

    // Special "automatic" heuristics
    if (autoDepth)
    {
        int depth = -1;

        if (autoWidth)
        {
            // Fully automatic:
            // Choose a depth that gives approximately targetOccupancy
            // objects per leaf in a balanced binary tree.

            depth = static_cast<int>
            (
                std::ceil(std::log2(double(total)/double(targetOccupancy)))
            );

            depth = std::clamp(depth, minDepth, maxDepthLimit);

            // The corresponding leaf occupancy threshold for this depth
            ctrl.minLeafSize = static_cast<int>
            (
                std::ceil(double(total) / std::exp2(depth))
            );
        }
        else  // (minLeafSize > 0)
        {
            // Leaf size explicitly specified:
            // derive an appropriate depth from it.
            depth = static_cast<int>
            (
                std::ceil(std::log2(double(total)/double(minLeafSize)))
            );

            depth = std::clamp(depth, minDepth, maxDepthLimit);
        }

        if (depth > 0)
        {
            // With small safety margin...
            ctrl.maxLevel = depth + depthSafety;
        }
    }
    else if (autoWidth)
    {
        // Depth explicitly specified:
        // derive the corresponding leaf occupancy.
        ctrl.minLeafSize = static_cast<int>
        (
            std::ceil(double(total)/std::exp2(maxLevel))
        );
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::AABBTreeBase::inflate(scalar factor)
{
    // Inflate box in case geometry reduces to 2-D
    if (factor > 0)
    {
        for (auto& bb : leafBoxes_)
        {
            bb.inflate(factor);
        }
    }
}


void Foam::AABBTreeBase::build_impl
(
    const label first,
    const label last,
    const UList<point>& points,
    DynamicList<label>& offsetsToBuild,
    const int currLevel,
    const control& ctrl
)
{
    #ifdef FULLDEBUG
    if (first < 0 || last > indices_.size() || first >= last)
    {
        FatalErrorInFunction
            << "Invalid range [" << first << ", " << last << ")"
            << " for indices of size " << indices_.size()
            << abort(FatalError);
    }

    if (currLevel == 0)
    {
        if (first != 0)
        {
            FatalErrorInFunction
                << "Partitioning " << first << ".." << last
                << " at level=0, should have first == 0 !!"
                << abort(FatalError);
        }
        else if (!offsetsToBuild.empty())
        {
            FatalErrorInFunction
                << "Partitioning " << first << ".." << last
                << " at level=0, should have an empty offset list!!"
                << abort(FatalError);
        }
    }
    else
    {
        if (offsetsToBuild.empty())
        {
            FatalErrorInFunction
                << "Offsets not initialized at level=" << currLevel
                << ", should have been called at level=0"
                << abort(FatalError);
        }
        else if (label prev = offsetsToBuild.back(); first != prev)
        {
            FatalErrorInFunction
                << "Partitioning " << first << ".." << last
                << " at level=" << currLevel
                << " but previous partitioning stopped at " << prev
                << abort(FatalError);
        }
    }
    #endif

    if (offsetsToBuild.empty() && currLevel == 0)
    {
        offsetsToBuild.push_back(0);
    }

    const label count = (last - first);

    const int minCount = std::max(1, ctrl.minLeafSize);

    if (currLevel >= ctrl.maxLevel || count/2 < minCount)
    {
        // End condition reached == this is a leaf
        offsetsToBuild.push_back(last);

        // Sort the per-leaf indices. Not completely necessary...
        std::sort(indices_.begin() + first, indices_.begin() + last);
    }
    else
    {
        treeBoundBox extentBb;
        for (auto i : indices_.slice(first, count))
        {
            extentBb.add(points[i]);
        }

        // Split direction: the longest side of the current range
        const direction cmpt = extentBb.maxDir();

        // Default is median split (equal counts)
        label pivot = (first + count/2);

        if (ctrl.algo == algorithm::geometric_split)
        {
            // Split at the geometric middle of the bounding box
            // (instead of the median element used for median_split).

            const scalar centreValue = extentBb.centre(cmpt);

            auto splitIter = std::partition
            (
                indices_.begin() + first,
                indices_.begin() + last,
                [&points, cmpt, centreValue](label a)
                {
                    return (points[a][cmpt] < centreValue);
                }
            );

            pivot = label(splitIter - indices_.begin());

            // Guard against a degenerate/skewed partitioning
            // - both children must still satisfy the minimum leaf size

            label repart_first = 0, repart_last = 0;

            if (label lowerLimit = (first + minCount); pivot < lowerLimit)
            {
                // Clamped partition - moves to the right
                repart_first = pivot;
                repart_last = last;
                pivot = lowerLimit;
            }
            else if (label upperLimit = (last - minCount); pivot > upperLimit)
            {
                // Clamped partition - moves to the left
                repart_first = first;
                repart_last = pivot;
                pivot = upperLimit;
            }

            if (repart_first != repart_last)
            {
                std::nth_element
                (
                    indices_.begin() + repart_first,
                    indices_.begin() + pivot,
                    indices_.begin() + repart_last,
                    [&points, cmpt](label a, label b)
                    {
                        return (points[a][cmpt] < points[b][cmpt]);
                    }
                );
            }
        }
        else  // algorithm::median_split
        {
            pivot = (first + count/2);

            std::nth_element
            (
                indices_.begin() + first,
                indices_.begin() + pivot,
                indices_.begin() + last,
                [&points, cmpt](label a, label b)
                {
                    return (points[a][cmpt] < points[b][cmpt]);
                }
            );
        }

        // Recurse left and then right. Order is essential!
        const int nextLevel(currLevel+1);

        build_impl(first, pivot, points, offsetsToBuild, nextLevel, ctrl);
        build_impl(pivot, last, points, offsetsToBuild, nextLevel, ctrl);
    }
}


void Foam::AABBTreeBase::clear()
{
    leafBoxes_.clear();
    clearAddressing();
}


void Foam::AABBTreeBase::clearAddressing()
{
    indices_.clear();
    offsets_.clear();
    nodes_.reset(nullptr);
}


void Foam::AABBTreeBase::clearNodes()
{
    nodes_.reset(nullptr);
}


void Foam::AABBTreeBase::reset
(
    const label nObjects,
    const UList<point>& points,
    bool medianSplit,
    int maxLevel,
    int minLeafSize
)
{
    clear();

    if (nObjects > 0)
    {
        control ctrl;
        ctrl.maxLevel = std::max(1, maxLevel);
        ctrl.minLeafSize = std::max(1, minLeafSize);
        ctrl.algo =
        (
            medianSplit
          ? algorithm::median_split
          : algorithm::geometric_split
        );

        // Special "automatic" heuristics
        autosize(ctrl, nObjects, maxLevel, minLeafSize);

        indices_ = Foam::identity(nObjects);

        DynamicList<label> offsetsToBuild;

        build_impl
        (
            0,                 // first
            indices_.size(),   // last
            points,
            offsetsToBuild,
            0,                 // starting at top level
            ctrl
        );

        // Top-level is done - transfer into permanent storage

        offsets_ = std::move(offsetsToBuild);
        const label numLeafs = nLeafs();

        if (numLeafs == 0)
        {
            // Degenerate case - probably never reached
            offsets_.clear();
            leafBoxes_.clear();
            return;
        }

        leafBoxes_.resize_nocopy(numLeafs);

        for (label leafi = 0; leafi < numLeafs; ++leafi)
        {
            auto& bb = leafBoxes_[leafi];
            bb.reset();

            for (auto i : indices(leafi))
            {
                bb.add(points[i]);
            }
        }
    }
}


void Foam::AABBTreeBase::calcNodeAddressing() const
{
    nodes_.reset(nullptr);
    nodes_ = std::make_unique<List<node>>();

    if (nLeafs() > 0)
    {
        DynamicList<node> nodesToBuild;

        reindex_impl(0, nLeafs(), nodesToBuild);

        *nodes_ = std::move(nodesToBuild);
    }
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::AABBTreeBase::AABBTreeBase
(
    const AABBTreeBase& other
)
:
    leafBoxes_(other.leafBoxes_),
    indices_(other.indices_),
    offsets_(other.offsets_)
{}


Foam::AABBTreeBase::AABBTreeBase
(
    AABBTreeBase&& other
)
:
    AABBTreeBase()
{
    transfer(other);
}


Foam::AABBTreeBase::AABBTreeBase
(
    const pointField& points,
    bool medianSplit,
    int maxLevel,
    int minLeafSize
)
{
    reset(points.size(), points, medianSplit, maxLevel, minLeafSize);
    inflate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::AABBTreeBase::pointInside(const point& pt) const
{
    if (leafBoxes_.empty())
    {
        return false;
    }

    if (cacheNodes() && !hasNodes())
    {
        calcNodeAddressing();
    }

    if (hasNodes())
    {
        return found_first_of
        (
            0,
            [&pt](const treeBoundBox& bb)
            {
                return bb.contains(pt);
            }
        );
    }
    else
    {
        // Flat searching
        for (const auto& bb : leafBoxes_)
        {
            if (bb.contains(pt))
            {
                return true;
            }
        }
    }

    return false;
}


bool Foam::AABBTreeBase::overlaps(const boundBox& queryBb) const
{
    if (leafBoxes_.empty() || !queryBb.good())
    {
        return false;
    }

    if (cacheNodes() && !hasNodes())
    {
        calcNodeAddressing();
    }

    if (hasNodes())
    {
        return found_first_of
        (
            0,
            [&queryBb](const treeBoundBox& bb)
            {
                return bb.overlaps(queryBb);
            }
        );
    }
    else
    {
        // Flat search
        for (const auto& bb : leafBoxes_)
        {
            if (bb.overlaps(queryBb))
            {
                return true;
            }
        }
    }

    return false;
}


void Foam::AABBTreeBase::writeOBJ(Ostream& os) const
{
    constexpr bool writeLinesOnly = false;

    label vertIndex(0);

    for (const auto& bb : leafBoxes_)
    {
        bb.writeOBJ(os, vertIndex, writeLinesOnly);
    }
}


void Foam::AABBTreeBase::transfer
(
    AABBTreeBase& other
)
{
    if (this == &other)
    {
        // no self-assignment
        return;
    }

    leafBoxes_ = std::move(other.leafBoxes_);
    indices_ = std::move(other.indices_);
    offsets_ = std::move(other.offsets_);
    nodes_ = std::move(other.nodes_);
}


void Foam::AABBTreeBase::reindex()
{
    calcNodeAddressing();
}


const Foam::List<Foam::AABBTreeBase::node>&
Foam::AABBTreeBase::nodes() const
{
    if (!nodes_ || (nodes_->empty() && nLeafs() > 0))
    {
        calcNodeAddressing();
    }

    return *nodes_;
}


Foam::treeBoundBox Foam::AABBTreeBase::bb() const
{
    treeBoundBox overallBb;

    for (const auto& leafBb : leafBoxes_)
    {
        overallBb.add(leafBb);
    }

    return overallBb;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::AABBTreeBase::operator=(const AABBTreeBase& other)
{
    if (this != &other)
    {
        // No self-assignment
        clear();

        leafBoxes_ = other.leafBoxes_;
        indices_ = other.indices_;
        offsets_ = other.offsets_;
    }
}


void Foam::AABBTreeBase::operator=(AABBTreeBase&& other)
{
    if (this != &other)
    {
        // No self-assignment
        transfer(other);
    }
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const AABBTreeBase& tree)
{
    os  << tree.leafBoxes()  // content
        << tree.offsets() << tree.indices();  // compact addressing

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, AABBTreeBase& tree)
{
    // Invalidate any navigation
    tree.clearNodes();

    is  >> tree.leafBoxes()  // content
        >> tree.offsets() >> tree.indices();  // compact addressing

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
