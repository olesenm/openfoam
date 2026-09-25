/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "treeBoundBox.H"
#include "Random.H"
#include "Pair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Source of truth is boundBox::octFaces()/octEdges() (construct-on-first-use).
// Kept here as references only for backward compatibility.
const Foam::faceList& Foam::treeBoundBox::faces = Foam::boundBox::octFaces();
const Foam::edgeList& Foam::treeBoundBox::edges = Foam::boundBox::octEdges();

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Code position of point on bounding box faces
inline static direction
get_faceBits(const boundBox& bb, const point& p) noexcept
{
    direction octant = 0;

    if (p.x() == bb.min().x())
    {
        octant |= treeBoundBox::LEFTBIT;
    }
    else if (p.x() == bb.max().x())
    {
        octant |= treeBoundBox::RIGHTBIT;
    }

    if (p.y() == bb.min().y())
    {
        octant |= treeBoundBox::BOTTOMBIT;
    }
    else if (p.y() == bb.max().y())
    {
        octant |= treeBoundBox::TOPBIT;
    }

    if (p.z() == bb.min().z())
    {
        octant |= treeBoundBox::BACKBIT;
    }
    else if (p.z() == bb.max().z())
    {
        octant |= treeBoundBox::FRONTBIT;
    }

    return octant;
}


// Position of point relative to bounding box
inline static direction
get_posBits(const boundBox& bb, const point& p) noexcept
{
    direction octant = 0;

    if (p.x() < bb.min().x())
    {
        octant |= treeBoundBox::LEFTBIT;
    }
    else if (p.x() > bb.max().x())
    {
        octant |= treeBoundBox::RIGHTBIT;
    }

    if (p.y() < bb.min().y())
    {
        octant |= treeBoundBox::BOTTOMBIT;
    }
    else if (p.y() > bb.max().y())
    {
        octant |= treeBoundBox::TOPBIT;
    }

    if (p.z() < bb.min().z())
    {
        octant |= treeBoundBox::BACKBIT;
    }
    else if (p.z() > bb.max().z())
    {
        octant |= treeBoundBox::FRONTBIT;
    }

    return octant;
}


static bool intersects_impl
(
    const boundBox& bb,
    const point& overallStart,
    const vector& overallVec,
    const point& start,
    const point& end,
    point& pt,
    direction& ptOnFaces
)
{
    // Sutherlands algorithm:
    //   loop
    //     - start = intersection of line with one of the planes bounding
    //       the bounding box
    //     - stop if start inside bb (return true)
    //     - stop if start and end in same 'half' (e.g. both above bb)
    //       (return false)
    //
    // Uses posBits to efficiently determine 'half' in which start and end
    // point are.
    //
    // Note:
    //   - sets coordinate to exact position: e.g. pt.x() = min().x()
    //     since plane intersect routine might have truncation error.
    //     This makes sure that posBits tests 'inside'

    const direction endBits = get_posBits(bb, end);
    pt = start;

    // Allow maximum of 3 clips.
    for (direction i = 0; i < 4; ++i)
    {
        const direction ptBits = get_posBits(bb, pt);

        if (ptBits == 0)
        {
            // pt inside bb
            ptOnFaces = get_faceBits(bb, pt);
            return true;
        }

        if ((ptBits & endBits) != 0)
        {
            // pt and end in same block outside of bb
            ptOnFaces = get_faceBits(bb, pt);
            return false;
        }

        if (ptBits & treeBoundBox::LEFTBIT)
        {
            // Intersect with plane V=min, n=-1,0,0
            pt.x() = bb.min().x();

            if (Foam::mag(overallVec.x()) > VSMALL)
            {
                scalar s = (pt.x() - overallStart.x())/overallVec.x();
                pt.y() = overallStart.y() + overallVec.y()*s;
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
        }
        else if (ptBits & treeBoundBox::RIGHTBIT)
        {
            // Intersect with plane V=max, n=1,0,0
            pt.x() = bb.max().x();

            if (Foam::mag(overallVec.x()) > VSMALL)
            {
                scalar s = (pt.x() - overallStart.x())/overallVec.x();
                pt.y() = overallStart.y() + overallVec.y()*s;
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
        }
        else if (ptBits & treeBoundBox::BOTTOMBIT)
        {
            // Intersect with plane V=min, n=0,-1,0
            pt.y() = bb.min().y();

            if (Foam::mag(overallVec.y()) > VSMALL)
            {
                scalar s = (pt.y() - overallStart.y())/overallVec.y();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
        }
        else if (ptBits & treeBoundBox::TOPBIT)
        {
            // Intersect with plane V=max, n=0,1,0
            pt.y() = bb.max().y();

            if (Foam::mag(overallVec.y()) > VSMALL)
            {
                scalar s = (pt.y() - overallStart.y())/overallVec.y();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
        }
        else if (ptBits & treeBoundBox::BACKBIT)
        {
            // Intersect with plane V=min, n=0,0,-1
            pt.z() = bb.min().z();

            if (Foam::mag(overallVec.z()) > VSMALL)
            {
                scalar s = (pt.z() - overallStart.z())/overallVec.z();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.y() = overallStart.y() + overallVec.y()*s;
            }
        }
        else if (ptBits & treeBoundBox::FRONTBIT)
        {
            // Intersect with plane V=max, n=0,0,1
            pt.z() = bb.max().z();

            if (Foam::mag(overallVec.z()) > VSMALL)
            {
                scalar s = (pt.z() - overallStart.z())/overallVec.z();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.y() = overallStart.y() + overallVec.y()*s;
            }
        }
    }

    // Can end up here if the end point is on the edge of the boundBox
    return true;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::treeBoundBox::intersects
(
    const boundBox& bb,
    const point& start,
    const point& end
)
{
    point pt;
    direction ptBits;
    return intersects_impl(bb, start, end-start, start, end, pt, ptBits);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeBoundBox::treeBoundBox(const UList<point>& points)
:
    boundBox(points, false)  // syncPar = false
{
    #ifdef FULLDEBUG
    if (points.empty())
    {
        WarningInFunction
            << "No bounding box for zero-sized pointField" << nl;
    }
    #endif
}


Foam::treeBoundBox::treeBoundBox
(
    const UList<point>& points,
    const labelUList& indices
)
:
    boundBox(points, indices, false)  // syncPar = false
{
    #ifdef FULLDEBUG
    if (points.empty() || indices.empty())
    {
        WarningInFunction
            << "No bounding box for zero-sized pointField" << nl;
    }
    #endif
}


Foam::treeBoundBox Foam::treeBoundBox::extend
(
    Random& rndGen,
    scalar factor
) const
{
    treeBoundBox bb(*this);
    bb.inflate(rndGen, factor);
    return bb;
}


Foam::treeBoundBox Foam::treeBoundBox::extend
(
    Random& rndGen,
    scalar factor,
    scalar delta
) const
{
    treeBoundBox bb(*this);
    bb.inflate(rndGen, factor, delta);
    return bb;
}


Foam::treeBoundBox Foam::treeBoundBox::extend_random
(
    label seedValue,
    scalar factor
) const
{
    Random rndGen(seedValue);

    treeBoundBox bb(*this);
    bb.inflate(rndGen, factor);
    return bb;
}


Foam::treeBoundBox Foam::treeBoundBox::extend_random
(
    label seedValue,
    scalar factor,
    scalar delta
) const
{
    Random rndGen(seedValue);

    treeBoundBox bb(*this);
    bb.inflate(rndGen, factor, delta);
    return bb;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox Foam::treeBoundBox::subBbox
(
    const point& mid,
    const direction octant
) const
{
    if (octant > 7)
    {
        FatalErrorInFunction
            << "octant:" << int(octant) << " should be [0..7]"
            << abort(FatalError);
    }

    // Start the box with a single point (the mid-point) and push out the
    // min/max dimensions according to the octant.

    treeBoundBox bb(mid);

    if (octant & treeBoundBox::RIGHTHALF)
    {
        bb.max().x() = max().x();
    }
    else
    {
        bb.min().x() = min().x();
    }

    if (octant & treeBoundBox::TOPHALF)
    {
        bb.max().y() = max().y();
    }
    else
    {
        bb.min().y() = min().y();
    }

    if (octant & treeBoundBox::FRONTHALF)
    {
        bb.max().z() = max().z();
    }
    else
    {
        bb.min().z() = min().z();
    }

    return bb;
}


Foam::treeBoundBox Foam::treeBoundBox::subHalf
(
    const scalar mid,
    const direction whichFace
) const
{
    // Start with a copy of this bounding box and adjust limits accordingly
    // - corresponds to a clipping plane

    treeBoundBox bb(*this);

    switch (whichFace)
    {
        case LEFT   : bb.max().x() = mid; break;
        case RIGHT  : bb.min().x() = mid; break;

        case BOTTOM : bb.max().y() = mid; break;
        case TOP    : bb.min().y() = mid; break;

        case BACK   : bb.max().z() = mid; break;
        case FRONT  : bb.min().z() = mid; break;

        default:
        {
            FatalErrorInFunction
                << "face:" << int(whichFace) << " should be [0..5]"
                << abort(FatalError);
        }
    }

    return bb;
}


Foam::treeBoundBox Foam::treeBoundBox::subHalf
(
    const direction whichFace
) const
{
    direction cmpt =
    (
        (whichFace == faceId::LEFT || whichFace == faceId::RIGHT)
      ? vector::X
      : (whichFace == faceId::BOTTOM || whichFace == faceId::TOP)
      ? vector::Y
      : vector::Z
    );

    scalar mid = 0.5*(min()[cmpt] + max()[cmpt]);

    return subHalf(mid, whichFace);
}


Foam::treeBoundBox Foam::treeBoundBox::subBbox(const direction octant) const
{
    return subBbox(centre(), octant);
}


bool Foam::treeBoundBox::subOverlaps
(
    const direction octant,
    const boundBox& bb
) const
{
    // Slightly accelerated version of
    //     subBbox(octant).overlaps(bb)

    point subMin = centre();
    point subMax = subMin;

    if (octant & RIGHTHALF)
    {
        subMax.x() = max().x();
    }
    else
    {
        subMin.x() = min().x();
    }

    if (octant & TOPHALF)
    {
        subMax.y() = max().y();
    }
    else
    {
        subMin.y() = min().y();
    }

    if (octant & FRONTHALF)
    {
        subMax.z() = max().z();
    }
    else
    {
        subMin.z() = min().z();
    }

    // NB: ordering of corners *is* irrelevant
    return box_box_overlaps(subMin, subMax, bb.min(), bb.max());
}


bool Foam::treeBoundBox::intersects
(
    const point& overallStart,
    const vector& overallVec,
    const point& start,
    const point& end,
    point& pt,
    direction& ptOnFaces
) const
{
    return intersects_impl
    (
        *this,
        overallStart,
        overallVec,
        start,
        end,
        pt,
        ptOnFaces
    );
}


bool Foam::treeBoundBox::intersects
(
    const point& start,
    const point& end,
    point& pt
) const
{
    direction ptBits;
    return intersects(start, end-start, start, end, pt, ptBits);
}


bool Foam::treeBoundBox::contains(const vector& dir, const point& pt) const
{
    // Compare all components against min and max of bb

    for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
    {
        if (pt[cmpt] < min()[cmpt])
        {
            // Outside
            return false;
        }
        else if (pt[cmpt] == min()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] < 0)
            {
                return false;
            }
        }

        if (pt[cmpt] > max()[cmpt])
        {
            // Outside
            return false;
        }
        else if (pt[cmpt] == max()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] > 0)
            {
                return false;
            }
        }
    }

    // All components inside bb
    return true;
}


Foam::direction Foam::treeBoundBox::faceBits(const point& pt) const
{
    return get_faceBits(*this, pt);
}


Foam::direction Foam::treeBoundBox::posBits(const point& pt) const
{
    return get_posBits(*this, pt);
}


Foam::Pair<Foam::direction>
Foam::treeBoundBox::findExtremities(const point& p) const
{
    direction nearCorner = 0;
    direction farCorner = 0;

    for (direction cmpt = 0; cmpt < point::nComponents; ++cmpt)
    {
        // ie, octantBit enumeration (RIGHTHALF | TOPHALF | FRONTHALF)
        direction octbit = (1 << cmpt);

        if
        (
            Foam::mag(min()[cmpt] - p[cmpt])
          < Foam::mag(max()[cmpt] - p[cmpt])
        )
        {
            farCorner |= octbit;
        }
        else
        {
            nearCorner |= octbit;
        }
    }

    return { nearCorner, farCorner };
}


void Foam::treeBoundBox::calcExtremities
(
    const point& pt,
    point& near,
    point& far
) const
{
    for (direction cmpt = 0; cmpt < point::nComponents; ++cmpt)
    {
        if
        (
            Foam::mag(min()[cmpt] - pt[cmpt])
          < Foam::mag(max()[cmpt] - pt[cmpt])
        )
        {
            near[cmpt] = min()[cmpt];
            far[cmpt]  = max()[cmpt];
        }
        else
        {
            near[cmpt] = max()[cmpt];
            far[cmpt]  = min()[cmpt];
        }
    }
}


// ************************************************************************* //
