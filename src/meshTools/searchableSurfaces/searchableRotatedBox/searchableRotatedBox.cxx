/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "searchableRotatedBox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableRotatedBox, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableRotatedBox,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableRotatedBox,
        dict,
        rotatedBox
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableRotatedBox::searchableRotatedBox
(
    const IOobject& io,
    const vector& span,
    const coordSystem::cartesian& csys
)
:
    searchableSurface(io),
    box_
    (
        IOobject
        (
            io.name() + "_box",
            io.instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false      // never register
        ),
        treeBoundBox(Zero, span)
    ),
    transform_(csys.origin(), csys.e3(), csys.e1())
{
    points_ = transform_.globalPosition(box_.points());
}


Foam::searchableRotatedBox::searchableRotatedBox
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableRotatedBox
    (
        io,
        dict.get<vector>("span"),
        coordSystem::cartesian
        (
            dict.get<point>("origin"),
            dict.get<vector>("e3"),
            dict.get<vector>("e1")
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableRotatedBox::regions() const
{
    return box_.regions();
}


Foam::tmp<Foam::pointField> Foam::searchableRotatedBox::coordinates() const
{
    return transform_.globalPosition(box_.coordinates());
}


void Foam::searchableRotatedBox::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    box_.boundingSpheres(centres, radiusSqr);
    centres = transform_.globalPosition(centres);
}


Foam::tmp<Foam::pointField> Foam::searchableRotatedBox::points() const
{
    return points_;
}


bool Foam::searchableRotatedBox::overlaps(const boundBox& bb) const
{
    // (from treeDataPrimitivePatch.C)

    // 1. bounding box
    if (!bb.overlaps(bounds()))
    {
        return false;
    }

    // 2. Check if any corner points inside
    if (bb.contains_any(points_))
    {
        return true;
    }

    // 3. Difficult case: all points are outside but connecting edges might
    // go through cube.

    // 3a. my edges through bb faces
    for (const auto& e : boundBox::octEdges())
    {
        if (bb.intersects(points_[e.start()], points_[e.end()]))
        {
            return true;
        }
    }

    // 3b. bb edges through my faces

    // NB: must use octant-ordered corners since using
    // treeBoundBox edges/faces (octEdges, octFaces)

    const pointField bbCorners(bb.octCorners());

    for (const auto& f : boundBox::octFaces())
    {
        const point fc = f.centre(points_);

        for (const auto& e : boundBox::octEdges())
        {
            pointHit inter = f.intersection
            (
                bbCorners[e.start()],
                bbCorners[e.end()],
                fc,
                points_,
                intersection::HALF_RAY
            );

            if (inter.hit() && inter.distance() <= 1)
            {
                return true;
            }
        }
    }

    return false;
}


Foam::pointIndexHit Foam::searchableRotatedBox::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit boxNearest
    (
        box_.findNearest
        (
            transform_.localPosition(sample),
            nearestDistSqr
        )
    );

    boxNearest.point() = transform_.globalPosition(boxNearest.point());

    return boxNearest;
}


Foam::pointIndexHit Foam::searchableRotatedBox::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    NotImplemented;
    return pointIndexHit();
}


Foam::pointIndexHit Foam::searchableRotatedBox::findLine
(
    const point& start,
    const point& end
) const
{
    pointIndexHit boxHit
    (
        box_.findLine
        (
            transform_.localPosition(start),
            transform_.localPosition(end)
        )
    );

    boxHit.point() = transform_.globalPosition(boxHit.point());

    return boxHit;
}


Foam::pointIndexHit Foam::searchableRotatedBox::findLineAny
(
    const point& start,
    const point& end
) const
{
    return findLine(start, end);
}


void Foam::searchableRotatedBox::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.resize_nocopy(samples.size());

    scalar maxDistSqr = searchableSurface::defaultMaxDistSqr();

    forAll(samples, i)
    {
        if (i < nearestDistSqr.size()) { maxDistSqr = nearestDistSqr[i]; }

        info[i] = findNearest(samples[i], maxDistSqr);
    }
}


void Foam::searchableRotatedBox::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize_nocopy(start.size());

    forAll(start, i)
    {
        info[i] = findLine(start[i], end[i]);
    }
}


void Foam::searchableRotatedBox::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize_nocopy(start.size());

    forAll(start, i)
    {
        info[i] = findLineAny(start[i], end[i]);
    }
}


void Foam::searchableRotatedBox::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.resize_nocopy(start.size());

    // Work array
    DynamicList<pointIndexHit> hits;

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (SMALL is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress

    forAll(start, pointi)
    {
        info[pointi].clear();

        const auto& startp = start[pointi];
        const auto& endp = end[pointi];

        // See if any intersection between pt and end
        pointIndexHit inter = findLine(startp, endp);

        if (inter.hit())
        {
            // For progression
            const vector dirVec = (endp - startp);
            const scalar magSqrDirVec = dirVec.magSqr();
            const vector smallVec
            (
                ROOTSMALL*dirVec + vector::uniform(ROOTVSMALL)
            );

            // Track last hit
            label lastIndex = inter.index();
            point pt = inter.point() + smallVec;

            hits.clear();
            hits.push_back(inter);

            while (((pt-startp) & dirVec) <= magSqrDirVec)
            {
                // See if any intersection between pt and end
                inter = findLine(pt, endp);

                // Check for not hit or hit same face as before (can happen
                // if vector along surface of face)
                if (inter.not_hit() || (inter.index() == lastIndex))
                {
                    break;
                }

                // update
                lastIndex = inter.index();
                pt = inter.point() + smallVec;

                hits.push_back(inter);
            }

            // copy (keep DynamicList for following loops)
            info[pointi] = hits;
        }
    }
}


void Foam::searchableRotatedBox::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.resize_fill(info.size(), 0);
}


void Foam::searchableRotatedBox::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    // searchableBox does not use hitPoints so no need to transform
    box_.getNormal(info, normal);

    normal = transform_.globalVector(normal);
}


Foam::label Foam::searchableRotatedBox::measure
(
    List<vector>* centres,
    List<scalar>* radiusSqr,
    List<label>* regions,
    List<vector>* normals
) const
{
    label returnCount = box_.measure(centres, radiusSqr, regions, normals);

    // Transform centres
    if (centres)
    {
        for (auto& fc : *centres)
        {
            fc = transform_.globalPosition(fc);
        }
    }

    // Transform normals
    if (normals)
    {
        for (auto& n : *normals)
        {
            n = transform_.globalVector(n);
        }
    }

    return returnCount;
}


void Foam::searchableRotatedBox::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    box_.getVolumeType(transform_.localPosition(points), volType);
}


// ************************************************************************* //
