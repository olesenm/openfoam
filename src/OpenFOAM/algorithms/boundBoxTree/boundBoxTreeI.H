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

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

inline Foam::boundBoxTree::boundBoxTree
(
    boundBoxTree::algorithm algo,
    int maxLevel,
    int minLeafSize
) noexcept
{
    setAlgorithm(algo);
    setMaxLevel(maxLevel);
    setMinLeafSize(minLeafSize);
}


inline Foam::boundBoxTree::boundBoxTree
(
    boundBoxTree::algorithm algo,
    List<treeBoundBox>&& boxes,
    int maxLevel,
    int minLeafSize
)
:
    boundBoxTree(algo, maxLevel, minLeafSize)
{
    autosize(ctrl_, boxes.size(), maxLevel, minLeafSize);

    setBoxes(std::move(boxes));
    rebuild();
}


inline Foam::boundBoxTree::boundBoxTree
(
    boundBoxTree::algorithm algo,
    const UList<treeBoundBox>& boxes,
    int maxLevel,
    int minLeafSize
)
:
    boundBoxTree(algo, maxLevel, minLeafSize)
{
    autosize(ctrl_, boxes.size(), maxLevel, minLeafSize);

    setBoxes(boxes);
    rebuild();
}


inline Foam::boundBoxTree::boundBoxTree
(
    boundBoxTree::algorithm algo,
    const UList<boundBox>& boxes,
    int maxLevel,
    int minLeafSize
)
:
    boundBoxTree(algo, maxLevel, minLeafSize)
{
    autosize(ctrl_, boxes.size(), maxLevel, minLeafSize);

    setBoxes(boxes);
    rebuild();
}


inline Foam::boundBoxTree::boundBoxTree
(
    List<treeBoundBox>&& boxes,
    int maxLevel,
    int minLeafSize
)
:
    boundBoxTree
    (
        boundBoxTree::defaultAlgorithm(),
        std::move(boxes),
        maxLevel,
        minLeafSize
    )
{}


inline Foam::boundBoxTree::boundBoxTree
(
    const UList<treeBoundBox>& boxes,
    int maxLevel,
    int minLeafSize
)
:
    boundBoxTree
    (
        boundBoxTree::defaultAlgorithm(),
        boxes,
        maxLevel,
        minLeafSize
    )
{}


inline Foam::boundBoxTree::boundBoxTree
(
    const UList<boundBox>& boxes,
    int maxLevel,
    int minLeafSize
)
:
    boundBoxTree
    (
        boundBoxTree::defaultAlgorithm(),
        boxes,
        maxLevel,
        minLeafSize
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::treeBoundBox&
Foam::boundBoxTree::bb() const
{
    return
    (
        nodes_.empty()
      ? treeBoundBox::null()
      : nodes_[root()].bb()
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::treeBoundBox&
Foam::boundBoxTree::operator[](label leafi) const
{
    return
    (
        (leafi >= 0 && leafi < leafBoxes_.size())
      ? leafBoxes_[leafi]
      : treeBoundBox::null()
    );
}


// * * * * * * * * * * * * * * * Query Functions * * * * * * * * * * * * * //

inline Foam::label Foam::boundBoxTree::first_overlaps
(
    const boundBox& queryBb
) const
{
    if (!nodes_.empty() && queryBb.good())
    {
        overlapsBoxOp pred{queryBb};
        return find_first_of(pred);
    }

    return -1;
}


inline bool Foam::boundBoxTree::find_overlaps
(
    const boundBox& queryBb,
    DynamicList<label>& hits
) const
{
    all_overlaps(queryBb, hits);
    return !hits.empty();
}


inline const Foam::labelUList&
Foam::boundBoxTree::all_overlaps
(
    const boundBox& queryBb,
    DynamicList<label>& hits
) const
{
    hits.clear();

    if (!nodes_.empty() && queryBb.good())
    {
        overlapsBoxOp pred{queryBb};
        find_all(pred, hits);
    }

    return hits;
}


inline Foam::label Foam::boundBoxTree::first_overlaps
(
    const point& centre,
    scalar radiusSqr
) const
{
    if (!nodes_.empty())
    {
        overlapsSphereOp pred{centre, radiusSqr};
        return find_first_of(pred);
    }

    return -1;
}


inline bool Foam::boundBoxTree::find_overlaps
(
    const point& centre,
    scalar radiusSqr,
    DynamicList<label>& hits
) const
{
    all_overlaps(centre, radiusSqr, hits);
    return !hits.empty();
}


inline const Foam::labelUList&
Foam::boundBoxTree::all_overlaps
(
    const point& centre,
    scalar radiusSqr,
    DynamicList<label>& hits
) const
{
    hits.clear();

    if (!nodes_.empty())
    {
        overlapsSphereOp pred{centre, radiusSqr};
        find_all(pred, hits);
    }

    return hits;
}


inline Foam::label Foam::boundBoxTree::first_contains
(
    const point& p
) const
{
    if (!nodes_.empty())
    {
        containsPointOp pred{p};
        return find_first_of(pred);
    }

    return -1;
}


inline bool Foam::boundBoxTree::find_contains
(
    const point& p,
    DynamicList<label>& hits
) const
{
    hits.clear();

    if (!nodes_.empty())
    {
        containsPointOp pred{p};
        return find_all(pred, hits);
    }

    return false;
}


inline Foam::label Foam::boundBoxTree::first_contains
(
    const point& centre,
    scalar radiusSqr
) const
{
    if (!nodes_.empty())
    {
        containsSphereOp pred{centre, radiusSqr};
        return find_first_of(pred);
    }

    return -1;
}


inline Foam::label Foam::boundBoxTree::first_contains
(
    const point& p0,
    const point& p1
) const
{
    if (!nodes_.empty())
    {
        containsLineOp pred{p0, p1};
        return find_first_of(pred);
    }

    return -1;
}


inline bool Foam::boundBoxTree::find_contains
(
    const point& p0,
    const point& p1,
    DynamicList<label>& hits
) const
{
    all_contains(p0, p1, hits);
    return !hits.empty();
}


inline const Foam::labelUList&
Foam::boundBoxTree::all_contains
(
    const point& p0,
    const point& p1,
    DynamicList<label>& hits
) const
{
    hits.clear();

    if (!nodes_.empty())
    {
        containsLineOp pred{p0, p1};
        find_all(pred, hits);
    }

    return hits;
}


inline Foam::label Foam::boundBoxTree::first_contains
(
    const boundBox& other
) const
{
    if (!nodes_.empty() && other.good())
    {
        containsBoxOp pred{other};
        return find_first_of(pred);
    }

    return -1;
}


inline bool Foam::boundBoxTree::find_contains
(
    const boundBox& other,
    DynamicList<label>& hits
) const
{
    hits.clear();

    if (!nodes_.empty() && other.good())
    {
        containsBoxOp pred{other};
        return find_all(pred, hits);
    }

    return false;
}


inline Foam::label Foam::boundBoxTree::first_intersects
(
    const point& from,
    const point& to
) const
{
    intersectsLineOp pred{from, to};
    return find_first_of(pred);
}


inline Foam::label Foam::boundBoxTree::first_intersects
(
    const linePointRef& ln
) const
{
    return first_intersects(ln.start(), ln.end());
}


inline bool Foam::boundBoxTree::find_intersects
(
    const point& from,
    const point& to,
    DynamicList<label>& hits
) const
{
    all_intersects(from, to, hits);
    return !hits.empty();
}


inline bool Foam::boundBoxTree::find_intersects
(
    const linePointRef& ln,
    DynamicList<label>& hits
) const
{
    all_intersects(ln.start(), ln.end(), hits);
    return !hits.empty();
}


inline const Foam::labelUList&
Foam::boundBoxTree::all_intersects
(
    const point& from,
    const point& to,
    DynamicList<label>& hits
) const
{
    intersectsLineOp pred{from, to};
    find_all(pred, hits);
    return hits;
}


inline const Foam::labelUList&
Foam::boundBoxTree::all_intersects
(
    const linePointRef& ln,
    DynamicList<label>& hits
) const
{
    return all_intersects(ln.start(), ln.end(), hits);
}


// ************************************************************************* //
