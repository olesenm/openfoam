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

#include "AABBTree.H"
#include "bitSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Predicate>
bool Foam::AABBTreeBase::found_first_of
(
    label nodei,
    const Predicate& pred
) const
{
    const node& nd = (*nodes_)[nodei];

    if (!pred(nd.bb()))
    {
        return false;
    }

    if (nd.is_leaf())
    {
        return true;
    }
    else  // interior:
    {
        // Do left first, right second
        if
        (
            found_first_of(nd.left(), pred)
         || found_first_of(nd.right(), pred)
        )
        {
            return true;
        }
    }

    return false;
}


template<class NodeType>
Foam::label Foam::AABBTreeBase::reindex_impl
(
    const label first,
    const label last,
    DynamicList<NodeType>& nodesToBuild
) const
{
    constexpr bool has_box =
    (
        std::is_base_of_v<Foam::boundBox, typename NodeType::box_type>
    );

    #ifdef FULLDEBUG
    if (first < 0 || last > nLeafs() || first >= last)
    {
        FatalErrorInFunction
            << "Invalid range [" << first << ", " << last << ")"
            << " for leaves of size " << nLeafs()
            << abort(FatalError);
    }
    #endif

    const label count = (last - first);

    // Create this node's slot before recursing etc.
    const label nodeId = nodesToBuild.size();
    nodesToBuild.emplace_back();

    if (count == 1)
    {
        // Leaf reached
        auto& leaf = nodesToBuild[nodeId];

        // The leaf start/size simply shadows the offsets_ values:
        leaf.set_leaf(offsets_[first], (offsets_[first+1]-offsets_[first]));

        if constexpr (has_box)
        {
            // 'first' corresponds to the leaf index
            leaf.set_bb(leafBoxes_[first]);
        }

        return nodeId;
    }
    else
    {
        // median split
        const label pivot = (first + count/2);

        // Recurse left and then right. Order is essential!
        const auto lefti  = reindex_impl(first, pivot, nodesToBuild);
        const auto righti = reindex_impl(pivot, last, nodesToBuild);

        nodesToBuild[nodeId].set_interior(lefti, righti);

        if constexpr (has_box)
        {
            // The interior node bounding box spans both children
            treeBoundBox nodeBb
            (
                nodesToBuild[lefti].bb(),
                nodesToBuild[righti].bb()
            );

            nodesToBuild[nodeId].set_bb(nodeBb);
        }
    }

    return nodeId;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::AABBTreeBase::object_inflate
(
    const UList<Type>& objects,
    const UList<point>& objectPoints,
    scalar factor
)
{
    const label numLeafs = nLeafs();

    for (label leafi = 0; leafi < numLeafs; ++leafi)
    {
        auto& bb = leafBoxes_[leafi];

        for (auto i : indices(leafi))
        {
            // Use the point references from the object

            for (label pointi : objects[i])
            {
                bb.add(objectPoints[pointi]);
            }
        }
    }

    // Inflate all final leaves
    this->inflate(factor);
}


template<class Type>
Foam::tmp<Foam::pointField>
Foam::AABBTreeBase::objectCentres
(
    const UList<Type>& objects,
    const pointField& objectPoints
)
{
    if constexpr (nonIndexable<Type>())
    {
        // Cannot use indirect point references
        // - the "object" already is the point itself
        return objectPoints;
    }
    else
    {
        // The object uses point references - define a representative
        // centre value for each

        auto tcentres = tmp<pointField>::New(objects.size(), Foam::zero{});
        auto& centres = tcentres.ref();

        bitSet pointsUsed(objectPoints.size());
        DynamicList<label> currPoints;

        label objecti = 0;
        for (const auto& obj : objects)
        {
            auto& centre = centres[objecti];
            ++objecti;

            currPoints.clear();

            for (label pointi : obj)
            {
                if (pointsUsed.set(pointi))
                {
                    currPoints.push_back(pointi);
                    centre += objectPoints[pointi];
                }
            }

            if (const label nObjPoints = currPoints.size(); nObjPoints > 0)
            {
                centre /= nObjPoints;

                // clear out bits before the next loop
                pointsUsed.unset(currPoints);
            }
        }

        return tcentres;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AABBTreeBase::AABBTreeBase
(
    const UList<Type>& objects,
    const pointField& objectCentres,
    const pointField& objectPoints,
    const bool medianSplit,
    int maxLevel,
    int minLeafSize
)
{
    reset(objects.size(), objectCentres, medianSplit, maxLevel, minLeafSize);
    object_inflate(objects, objectPoints);
}


template<class Type>
Foam::AABBTreeBase::AABBTreeBase
(
    const UList<Type>& objects,
    const pointField& objectPoints,
    const bool medianSplit,
    int maxLevel,
    int minLeafSize
)
{
    const tmp<pointField> tcentres(objectCentres(objects, objectPoints));

    reset(objects.size(), tcentres(), medianSplit, maxLevel, minLeafSize);

    if constexpr (nonIndexable<Type>())
    {
        inflate();
    }
    else
    {
        object_inflate(objects, objectPoints);
    }
}


template<class Type>
Foam::AABBTree<Type>::AABBTree
(
    const pointField& points,
    const bool medianSplit,
    int maxLevel,
    int minLeafSize
)
:
    AABBTreeBase(points, medianSplit, maxLevel, minLeafSize)
{}


template<class Type>
Foam::AABBTree<Type>::AABBTree
(
    const UList<Type>& objects,
    const pointField& objectCentres,
    const pointField& objectPoints,
    const bool medianSplit,
    int maxLevel,
    int minLeafSize
)
:
    AABBTreeBase
    (
        objects,
        objectCentres,
        objectPoints,
        medianSplit, maxLevel, minLeafSize
    )
{}


template<class Type>
Foam::AABBTree<Type>::AABBTree
(
    const UList<Type>& objects,
    const pointField& objectPoints,
    const bool medianSplit,
    int maxLevel,
    int minLeafSize
)
:
    AABBTreeBase
    (
        objects,
        objectPoints,
        medianSplit, maxLevel, minLeafSize
    )
{}


// ************************************************************************* //
