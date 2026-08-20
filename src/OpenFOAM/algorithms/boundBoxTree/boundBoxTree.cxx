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

\*---------------------------------------------------------------------------*/

#include "boundBoxTree.H"
#include "predicates.H"
#include <algorithm>
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace
{
    // Default constructed has addressing (0, -1) -> an empty leaf,
    // and the default box is inverted box
    const Foam::boundBoxTree::node emptyLeaf_singleton_;
}

const Foam::boundBoxTree::node& Foam::boundBoxTree::empty_leaf()
{
    return emptyLeaf_singleton_;
}

static constexpr bool debuggingOn = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// The distance of an outside point to the nearest box corner.
// Is zero for an inside point
Foam::scalar pointToBoxDistSqr
(
    const Foam::boundBox& bbox,
    const Foam::point& p
)
{
    Foam::scalar distSqr = 0;

    for
    (
        Foam::direction cmpt = 0;
        cmpt < Foam::vector::nComponents;
        ++cmpt
    )
    {
        Foam::scalar d = 0;

        if (p[cmpt] < bbox.min()[cmpt])
        {
            // An outside point : distance to min corner
            d = bbox.min()[cmpt] - p[cmpt];
        }
        else if (p[cmpt] > bbox.max()[cmpt])
        {
            // An outside point : distance to max corner
            d = p[cmpt] - bbox.max()[cmpt];
        }

        distSqr += Foam::sqr(d);
    }

    return distSqr;
}

} // End namespace Foam


// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * //

void Foam::boundBoxTree::autosize
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
                std::ceil(double(total)/std::exp2(depth))
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

    if constexpr (debuggingOn)
    {
        InfoErr
            << "boundBoxTree: autosize =>"
            << " algo:" << int(ctrl.algo)
            << " maxLevel:" << ctrl.maxLevel
            << " minLeaf:" << ctrl.minLeafSize << endl;
    }
}


void Foam::boundBoxTree::writeOBJ
(
    const label nodeId,
    Ostream& os,
    label& vertIndex,
    bool leavesOnly,
    bool writeLinesOnly
) const
{
    const node& nd = nodes_[nodeId];

    if (nd.is_leaf())
    {
        if (nd.size() <= 1)
        {
            // The containing bound as solid (faces)
            nd.bb().writeOBJ(os, vertIndex, false);
        }
        else
        {
            // The containing bound as lines
            nd.bb().writeOBJ(os, vertIndex, true);

            for (auto i : nd.slice(indices_))
            {
                // Individual boxes as solid (faces)
                leafBoxes_[i].writeOBJ(os, vertIndex, false);
            }
        }
    }
    else
    {
        writeOBJ(nd.left(), os, vertIndex, leavesOnly, writeLinesOnly);
        writeOBJ(nd.right(), os, vertIndex, leavesOnly, writeLinesOnly);
    }
}


Foam::label Foam::boundBoxTree::build_impl
(
    const label first,
    const label last,
    DynamicList<node>& nodesToBuild,
    const int currLevel
)
{
    #ifdef FULLDEBUG
    // Internal method, so first/last should never violate these conditions.
    // If they do, it is a programming error
    if (first < 0 || last > indices_.size() || first >= last)
    {
        FatalErrorInFunction
            << "Invalid range [" << first << ", " << last
            << ") for indices of size " << indices_.size()
            << abort(FatalError);
    }
    #endif

    const label count = (last - first);

    // Create this node's slot before recursing etc.
    const label nodeId = nodesToBuild.size();
    nodesToBuild.emplace_back();

    const int minCount = std::max(1, ctrl_.minLeafSize);

    if (currLevel >= ctrl_.maxLevel || (count < 2*minCount))
    {
        // End condition reached == this is a leaf
        auto& leaf = nodesToBuild[nodeId];

        // Bounding box enclosing all entries in [first, last)
        treeBoundBox leafBb;
        for (auto i : indices_.slice(first, count))
        {
            leafBb.add(leafBoxes_[i]);
        }

        // Sorting the per-leaf indices is not completely necessary...
        std::sort(indices_.begin() + first, indices_.begin() + last);

        leaf.set_bb(leafBb);
        leaf.set_leaf(first, count);

        return nodeId;
    }

    // Default is median split (equal counts)
    label pivot = (first + count/2);

    // The split plane, initialized with something invalid
    direction splitPlane(8);

    switch (ctrl_.algo)
    {
        case algorithm::median_split :
        {
            // Use the box centres
            treeBoundBox overallBb;
            for (auto i : indices_.slice(first, count))
            {
                overallBb.add(leafBoxes_[i].centre());
            }

            splitPlane = overallBb.maxDir();
            break;
        }

        case algorithm::box_split :
        {
            // Use the entire box extents
            treeBoundBox overallBb;
            for (auto i : indices_.slice(first, count))
            {
                overallBb.add(leafBoxes_[i]);
            }

            splitPlane = overallBb.maxDir();
            break;
        }

        case algorithm::geometric_split :
        {
            // Split at the geometric middle of the bounding box
            // (instead of the median element used for median_split).

            // Use the box centres
            treeBoundBox overallBb;
            for (auto i : indices_.slice(first, count))
            {
                overallBb.add(leafBoxes_[i].centre());
            }

            // Split direction: the longest side of the current range
            const direction cmpt = overallBb.maxDir();

            const scalar centreValue = overallBb.centre(cmpt);

            auto splitIter = std::partition
            (
                indices_.begin() + first,
                indices_.begin() + last,
                [this, cmpt, centreValue](label a)
                {
                    return (leafBoxes_[a].centre(cmpt) < centreValue);
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
                    [this, cmpt](label a, label b)
                    {
                        return
                        (
                            leafBoxes_[a].centre(cmpt)
                          < leafBoxes_[b].centre(cmpt)
                        );
                    }
                );
            }
            break;
        }

        case algorithm::pre_balanced :
        {
            // Trust the caller-supplied grouping - simply bisect the range.
            break;
        }

        default :
        {
            FatalErrorInFunction
                << "Unknown splitting algorithm: " << int(ctrl_.algo)
                << abort(FatalError);
            break;
        }
    }

    // Partition in the split plane by the box centres
    if (splitPlane <= vector::Z)
    {
        const direction cmpt = splitPlane;

        std::nth_element
        (
            indices_.begin() + first,
            indices_.begin() + pivot,
            indices_.begin() + last,
            [this, cmpt](label a, label b)
            {
                return
                (
                    leafBoxes_[a].centre(cmpt)
                  < leafBoxes_[b].centre(cmpt)
                );
            }
        );
    }

    // Recurse left and then right. Order is essential!
    const int nextLevel(currLevel+1);

    const auto lefti  = build_impl(first, pivot, nodesToBuild, nextLevel);
    const auto righti = build_impl(pivot, last, nodesToBuild, nextLevel);

    // The interior node bounding box spans both children
    treeBoundBox nodeBb
    (
        nodesToBuild[lefti].bb(),
        nodesToBuild[righti].bb()
    );

    nodesToBuild[nodeId].set_bb(nodeBb);
    nodesToBuild[nodeId].set_interior(lefti, righti);

    return nodeId;
}


Foam::label Foam::boundBoxTree::countElements(label nodei) const
{
    if (nodei < 0 || nodei >= nodes_.size())
    {
        // Protect against bad input (exposed method)
        return 0;
    }

    return count_if(nodei, predicates::always{});
}


Foam::label Foam::boundBoxTree::countLeafs(label nodei) const
{
    if (nodei < 0 || nodei >= nodes_.size())
    {
        // Protect against bad input (exposed method)
        return 0;
    }

    const node& nd = nodes_[nodei];

    if (nd.is_leaf())
    {
        return nd.size();
    }
    else  // is_interior
    {
        return
        (
            countLeafs(nd.left())
          + countLeafs(nd.right())
        );
    }
}


Foam::List<Foam::labelRange>
Foam::boundBoxTree::ranges() const
{
    List<labelRange> values;

    if (const label total = nLeafs(); total > 0)
    {
        values.resize(total);
        label leafi = 0;

        visit_leaf
        (
            [&](const node& nd)
            {
                values[leafi++] = labelRange(nd.start(), nd.size());
            }
        );
    }

    return values;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundBoxTree::setBoxes
(
    List<treeBoundBox>&& boxes
)
{
    clear();
    leafBoxes_ = std::move(boxes);
}


void Foam::boundBoxTree::setBoxes
(
    const UList<treeBoundBox>& boxes
)
{
    clear();
    leafBoxes_ = boxes;
}


void Foam::boundBoxTree::setBoxes
(
    const UList<boundBox>& boxes
)
{
    clear();
    leafBoxes_.resize_nocopy(boxes.size());
    std::copy(boxes.begin(), boxes.end(), leafBoxes_.begin());
}


void Foam::boundBoxTree::rebuild()
{
    // Effectively clearAddressing() but reclaiming indices_ if possible
    nodes_.clear();
    indices_.resize_nocopy(leafBoxes_.size());

    // Identity map, but only consider good bounding boxes
    label count = 0;
    forAll(leafBoxes_, i)
    {
        if (leafBoxes_[i].good())
        {
            indices_[count] = i;
            ++count;
        }
    }
    indices_.resize(count);

    if (count)
    {
        DynamicList<node> nodesToBuild(2*count);

        // If input order is trusted to already reflect the desired
        // spatial grouping, indices_ is an identity permutation and we
        // can build by simply splitting the indices

        const auto oldAlgo = ctrl_.algo;

        if
        (
            (oldAlgo == algorithm::pre_balanced)
         && (count < leafBoxes_.size())
        )
        {
            // The requirements of "pre_balanced" did not hold,
            // use the default algorithm instead
            ctrl_.algo = boundBoxTree::defaultAlgorithm();
        }

        constexpr int currLevel = 0;

        // NB: use count (number of valid boxes), not leafBoxes_.size() !!
        // since invalid boxes have been excluded above
        build_impl(0, count, nodesToBuild, currLevel);

        nodes_.transfer(nodesToBuild);

        // Restore algorithm if needed
        ctrl_.algo = oldAlgo;
    }
}


void Foam::boundBoxTree::clearAddressing()
{
    nodes_.clear();
    indices_.clear();
}


void Foam::boundBoxTree::clear()
{
    leafBoxes_.clear();
    clearAddressing();
}


Foam::label Foam::boundBoxTree::nInterior() const
{
    label total = 0;

    for (const auto& nd : nodes_)
    {
        if (nd.is_interior())
        {
            ++total;
        }
    }

    return total;
}


Foam::label Foam::boundBoxTree::nLeafs() const
{
    label total = 0;

    for (const auto& nd : nodes_)
    {
        if (nd.is_leaf())
        {
            ++total;
        }
    }

    return total;
}


// * * * * * * * * * * * * * * * Search Functions  * * * * * * * * * * * * * //

Foam::label Foam::boundBoxTree::findSphere
(
    const point& centre,
    const scalar radiusSqr,
    DynamicList<label>& hits
) const
{
    hits.clear();

    if (nodes_.empty())
    {
        return 0;
    }

    DynamicList<label> stack(64);
    stack.push_back(root());

    while (!stack.empty())
    {
        // Process nodes as LIFO
        label nodei = stack.back();
        stack.pop_back();

        const node& nd = nodes_[nodei];
        const auto& bb = nd.bb();

        if (!bb.overlaps(centre, radiusSqr))
        {
            continue;
        }

        if (nd.is_leaf())
        {
            for (auto i : nd.slice(indices_))
            {
                if (leafBoxes_[i].overlaps(centre, radiusSqr))
                {
                    hits.push_back(i);
                }
            }
        }
        else  // is_interior
        {
            label node0 = nd.left();
            label node1 = nd.right();

            if
            (
                treeBoundBox::distanceCmp
                (
                    nodes_[node0].bb(),
                    nodes_[node1].bb(),
                    centre
                ) > 0
            )
            {
                // visit right then left
                stack.push_back(node0); stack.push_back(node1);
            }
            else
            {
                // visit left then right
                stack.push_back(node1); stack.push_back(node0);
            }
        }
    }


    // TBD: indices sorted by nearest to sphere
    if constexpr (false)
    {
        std::sort
        (
            hits.begin(),
            hits.end(),

            [&centre,this](label a, label b)
            {
                return
                (
                    leafBoxes_[a].minDistSqr(centre)
                  < leafBoxes_[b].minDistSqr(centre)
                );
            }
        );
    }

    return hits.size();
}


Foam::Tuple2<Foam::scalar, Foam::volumeType>
Foam::boundBoxTree::findNearest
(
    const point& centre,
    const scalar radiusSqr,
    DynamicList<label>& hits
) const
{
    hits.clear();

    if (nodes_.empty())
    {
        // Empty tree - nothing is known
        return { Foam::sqr(GREAT), volumeType::UNKNOWN };
    }

    // We want to capture all possible contains() and overlaps(),
    // but if there are contains() we discard any overlaps() and only
    // retain the contained values.

    // The current scanning limits, volumeType also represents the state engine
    label bestIndex = -1;
    volumeType state = volumeType::OUTSIDE;
    scalar minDistSqr = radiusSqr;
    scalar maxDistSqr = radiusSqr;

    // Special handling for zero radius - restrict to contains() only
    // (no overlap checking). If nothing was actually found, bestIndex
    // will still be negative

    const bool degenerate = (radiusSqr <= VSMALL);
    if (degenerate)
    {
        // Toggle into INSIDE (contains) state to prevent overlaps checks.
        state = volumeType::INSIDE;
        minDistSqr = Foam::sqr(GREAT);
        maxDistSqr = Foam::sqr(GREAT);
    }

    DynamicList<label> stack(64);
    stack.push_back(root());

    while (!stack.empty())
    {
        // Process nodes as LIFO
        label nodei = stack.back();
        stack.pop_back();

        const node& nd = nodes_[nodei];

        // Prune: nothing under this node can beat the current best
        if
        (
            const auto& bbox = nd.bb();
            (
                state == volumeType::INSIDE
              ? !bbox.contains(centre)
              : !bbox.overlaps(centre, maxDistSqr)
            )
        )
        {
            continue;
        }

        if (nd.is_leaf())
        {
            for (auto i : nd.slice(indices_))
            {
                const auto& bbox = leafBoxes_[i];

                if (bbox.contains(centre))
                {
                    if (state != volumeType::INSIDE)
                    {
                        // Toggle permanently into INSIDE (contains) state
                        hits.clear();
                        bestIndex = i;
                        state = volumeType::INSIDE;
                        minDistSqr = Foam::sqr(GREAT);
                        maxDistSqr = Foam::sqr(GREAT);
                    }
                    hits.push_back(i);

                    // The smallest sphere that encloses this box, thus using
                    // the maxDistSqr() and not minDistSqr() method
                    if
                    (
                        scalar d2 = bbox.maxDistSqr(centre);
                        (d2 < minDistSqr)
                    )
                    {
                        bestIndex = i;
                        minDistSqr = d2;
                    }
                }
                else if
                (
                    // Check overlaps while not yet in INSIDE (contains) state
                    (state != volumeType::INSIDE)
                 && bbox.overlaps(centre, maxDistSqr)
                )
                {
                    state = volumeType::MIXED;
                    hits.push_back(i);

                    // min/max distances squared
                    const Pair<scalar> dsqr = bbox.distSqr(centre);

                    // Use the true point-to-AABB minimum distance for
                    // candidate selection and hit ordering. Ensure an
                    // inclusive first overlap can establish the nearest
                    // candidate at an equal search bound.

                    const scalar d2 = pointToBoxDistSqr(bbox, centre);

                    if
                    (
                        (bestIndex < 0)
                     || (d2 < minDistSqr)
                    )
                    {
                        bestIndex = i;
                        minDistSqr = d2;
                        maxDistSqr = Foam::min(maxDistSqr, dsqr.second());
                    }
                }
            }
        }
        else  // is_interior
        {
            label node0 = nd.left();
            label node1 = nd.right();

            // Visiting the closer child first (tightens minDistSqr sooner)
            if
            (
                treeBoundBox::distanceCmp
                (
                    nodes_[node0].bb(),
                    nodes_[node1].bb(),
                    centre
                ) > 0
            )
            {
                stack.push_back(node0); stack.push_back(node1);
            }
            else
            {
                stack.push_back(node1); stack.push_back(node0);
            }
        }
    }


    // Reprune/filter and sort by distance.
    // If we only had overlaps (no contains) also need to re-prune according
    // to the smallest maxDistSqr

    if (bestIndex < 0)
    {
        hits.clear();
        state = volumeType::OUTSIDE;
    }
    else if (hits.size() > 1)  // also implies that (bestIndex >= 0)
    {
        // Reuse stack for sorting etc.
        stack.clear();
        auto& idxs = stack;
        DynamicList<scalar> vals(hits.size());

        if (state == volumeType::INSIDE)
        {
            for (auto i : hits)
            {
                const auto& bbox = leafBoxes_[i];

                idxs.push_back(i);
                vals.push_back(bbox.maxDistSqr(centre));
            }
        }
        else
        {
            for (auto i : hits)
            {
                const auto& bbox = leafBoxes_[i];

                if
                (
                    (i == bestIndex)
                 || bbox.overlaps(centre, maxDistSqr)
                )
                {
                    idxs.push_back(i);
                    vals.push_back(pointToBoxDistSqr(bbox, centre));
                }
            }
        }

        // Temporarily use 'hits' for the sorting indices
        auto& sortOrder = hits;
        sortOrder.resize_nocopy(vals.size());
        Foam::identity(sortOrder);

        std::stable_sort
        (
            sortOrder.begin(),
            sortOrder.end(),
            ListOps::less(vals)
        );

        // Inplace transcribe in the correct order
        for (auto& x : sortOrder)
        {
            x = idxs[x];
        }
    }

    // The tightened search bound (radius squared) that still guarantees a hit:
    // - INSIDE -> smallest sphere enclosing a containing box
    // - MIXED  -> far distance of the nearest box. Never above the input.

    scalar boundSqr = radiusSqr;

    if (degenerate)
    {
        boundSqr = 0;
    }
    else if (bestIndex >= 0)
    {
        boundSqr =
        (
            (state == volumeType::INSIDE)
          ? std::min(radiusSqr, minDistSqr)
          : maxDistSqr
        );
    }

    return { boundSqr, state };
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundBoxTree::writeOBJ(Ostream& os) const
{
    constexpr bool leavesOnly = true;
    constexpr bool writeLinesOnly = false;

    if (!nodes_.empty())
    {
        label vertIndex(0);

        writeOBJ(root(), os, vertIndex, leavesOnly, writeLinesOnly);
    }
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBoxTree& tree)
{
    os << tree.boxes() << endl;
    os << flatOutput(tree.indices()) << nl;

    os << "root:" << tree.root() << endl;
    os << "bb:" << tree.bb() << endl;

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<boundBoxTree>& iproxy
)
{
    const auto& tree = *iproxy;

    // Visitor collecting depth and leaf occupancy statistics
    struct statsVisitor
    {
        // Remaining unvisited children for each open interior node
        DynamicList<int> pending;

        label nLeafs = 0;
        label nEmpty = 0;
        label nItems = 0;
        label minOccupancy = labelMax;
        label maxOccupancy = 0;
        label minDepth = labelMax;
        label maxDepth = 0;

        bool interior(const treeBoundBox&)
        {
            pending.push_back(2);
            return true;
        }

        void leaf(const treeBoundBox&, const labelUList& indices)
        {
            const label depth = pending.size();
            const label n = indices.size();

            ++nLeafs;
            nItems += n;
            if (!n) ++nEmpty;
            minOccupancy = Foam::min(minOccupancy, n);
            maxOccupancy = Foam::max(maxOccupancy, n);
            minDepth = Foam::min(minDepth, depth);
            maxDepth = Foam::max(maxDepth, depth);

            // Leaf complete: close any parents with all children done
            while (!pending.empty())
            {
                if (--pending.back() > 0) break;
                pending.pop_back();
            }
        }
    };

    statsVisitor stats;
    if (tree.hasNodes())
    {
        tree.visit(stats);
    }

    // Control parameters
    os  << "boxes: " << tree.boxes().size()
        << " indices:" << tree.indices().size()
        << " nodes:" << tree.nNodes()
        << " interior:" << tree.nInterior()
        << " leafs:" << tree.nLeafs()
        << nl;

    os << "control { max-level=" << tree.maxLevel()
        << " min-leaf-size=" << tree.minLeafSize()
        << " }" << nl;

    if (stats.nLeafs)
    {
        os  << "stats { depth(min/max)=" << stats.minDepth
            << '/' << stats.maxDepth
            << " occupancy(min/max/avg)=" << stats.minOccupancy
            << '/' << stats.maxOccupancy
            << '/' << (scalar(stats.nItems)/stats.nLeafs)
            << " empty-leafs=" << stats.nEmpty
            << " }" << nl;
    }

    return os;
}


// ************************************************************************* //
