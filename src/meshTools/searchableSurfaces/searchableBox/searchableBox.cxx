/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "searchableBox.H"
#include "addToRunTimeSelectionTable.H"
#include <array>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableBox, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableBox,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableBox,
        dict,
        box
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read min/max or min/span
static void readBoxDim(const dictionary& dict, treeBoundBox& bb)
{
    IOobjectOption::readOption readOpt = IOobjectOption::MUST_READ;

    dict.readEntry("min", bb.min(), keyType::LITERAL, readOpt);

    if (dict.readIfPresent("span", bb.max(), keyType::LITERAL))
    {
        bb.max() += bb.min();
        readOpt = IOobjectOption::READ_IF_PRESENT;
    }

    dict.readEntry("max", bb.max(), keyType::LITERAL, readOpt);
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableBox::projectOntoCoordPlane
(
    const direction dir,
    const point& planePt,
    pointIndexHit& info
) const
{
    // Set point
    info.point()[dir] = planePt[dir];

    // Set face
    if (planePt[dir] == min()[dir])
    {
        info.setIndex(dir*2);
    }
    else if (planePt[dir] == max()[dir])
    {
        info.setIndex(dir*2+1);
    }
    else
    {
        FatalErrorInFunction
            << "Point on plane " << planePt
            << " is not on coordinate " << min()[dir]
            << " nor " << max()[dir] << nl
            << abort(FatalError);
    }
}


// Returns miss or hit with face (0..5) and region(always 0)
void Foam::searchableBox::findNearest_impl
(
    const point& bbMid,
    const point& sample,
    const scalar nearestDistSqr,
    pointIndexHit& info
) const
{
    // Point can be inside or outside. For every component direction can be
    // left of min, right of max or inbetween.
    // - outside points: project first one x plane (either min().x()
    // or max().x()), then onto y plane and finally z. You should be left
    // with intersection point
    // - inside point: find nearest side (compare to mid point). Project onto
    //   that.

    // The face is set to the last projected face.


    // Outside point projected onto cube. Assume faces 0..5.
    info.hitPoint(sample, -1);
    bool outside = false;

    // (for internal points) per direction what nearest cube side is
    point near;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        if (info.point()[dir] < min()[dir])
        {
            projectOntoCoordPlane(dir, min(), info);
            outside = true;
        }
        else if (info.point()[dir] > max()[dir])
        {
            projectOntoCoordPlane(dir, max(), info);
            outside = true;
        }
        else if (info.point()[dir] > bbMid[dir])
        {
            near[dir] = max()[dir];
        }
        else
        {
            near[dir] = min()[dir];
        }
    }


    // For outside points the info will be correct now. Handle inside points
    // using the three near distances. Project onto the nearest plane.
    if (!outside)
    {
        const vector dist(cmptMag(info.point() - near));

        direction projNorm(vector::Z);

        if (dist.x() < dist.y())
        {
            if (dist.x() < dist.z())
            {
                projNorm = vector::X;
            }
        }
        else
        {
            if (dist.y() < dist.z())
            {
                projNorm = vector::Y;
            }
        }

        projectOntoCoordPlane(projNorm, near, info);
    }


    // Check if outside. Optimisation: could do some checks on distance already
    // on components above
    if (info.point().distSqr(sample) > nearestDistSqr)
    {
        info.setMiss();
        info.setIndex(-1);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableBox::searchableBox
(
    const IOobject& io,
    const treeBoundBox& bb
)
:
    searchableSurface(io),
    treeBoundBox(bb)
{
    if (!treeBoundBox::good())
    {
        FatalErrorInFunction
            << "Illegal bounding box specification : "
            << static_cast<const treeBoundBox>(*this) << nl
            << exit(FatalError);
    }

    bounds() = static_cast<treeBoundBox>(*this);
}


Foam::searchableBox::searchableBox
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    treeBoundBox()
{
    readBoxDim(dict, *this);

    if (!treeBoundBox::good())
    {
        FatalErrorInFunction
            << "Illegal bounding box specification : "
            << static_cast<const treeBoundBox>(*this) << nl
            << exit(FatalError);
    }

    bounds() = static_cast<treeBoundBox>(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableBox::regions() const
{
    if (regions_.empty())
    {
        regions_.resize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


Foam::tmp<Foam::pointField> Foam::searchableBox::coordinates() const
{
    auto tctrs = tmp<pointField>::New(treeBoundBox::nFaces());
    treeBoundBox::faceCentres(tctrs.ref());

    return tctrs;
}


void Foam::searchableBox::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.resize_nocopy(treeBoundBox::nFaces());
    radiusSqr.resize_nocopy(treeBoundBox::nFaces());

    // Box faces are quads, points are axis-aligned,
    // so the face centre is trivial and the only need one face corner
    // to determine the bounding sphere radius

    for (label facei = 0; facei < treeBoundBox::nFaces(); ++facei)
    {
        label pointi = treeBoundBox::octFaces()[facei][0];
        point corner = treeBoundBox::octCorner(pointi);

        // Add a bit to make sure all points are tested inside
        centres[facei] = treeBoundBox::faceCentre(facei);
        radiusSqr[facei] = (centres[facei].distSqr(corner) + Foam::sqr(SMALL));
    }
}


Foam::tmp<Foam::pointField> Foam::searchableBox::points() const
{
    return treeBoundBox::points();
}


Foam::pointIndexHit Foam::searchableBox::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info;
    findNearest_impl(centre(), sample, nearestDistSqr, info);
    return info;

}


Foam::pointIndexHit Foam::searchableBox::findNearestOnEdge
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    const point bbMid(centre());

    // Outside point projected onto cube. Assume faces 0..5.
    pointIndexHit info(true, sample, -1);
    bool outside = false;

    // (for internal points) per direction what nearest cube side is
    point near;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        if (info.point()[dir] < min()[dir])
        {
            projectOntoCoordPlane(dir, min(), info);
            outside = true;
        }
        else if (info.point()[dir] > max()[dir])
        {
            projectOntoCoordPlane(dir, max(), info);
            outside = true;
        }
        else if (info.point()[dir] > bbMid[dir])
        {
            near[dir] = max()[dir];
        }
        else
        {
            near[dir] = min()[dir];
        }
    }


    // For outside points the info will be correct now. Handle inside points
    // using the three near distances. Project onto the nearest two planes.
    if (!outside)
    {
        // Get the per-component distance to nearest wall
        const vector dist(cmptMag(info.point() - near));

        // The sorted order (by magnitude)
        std::array<direction, vector::nComponents> indices{0, 1, 2};

        std::stable_sort
        (
            indices.begin(),
            indices.end(),
            [&dist](auto a, auto b){ return dist[a] < dist[b]; }
        );

        // Project onto nearest
        projectOntoCoordPlane(indices[0], near, info);
        // Project onto second nearest
        projectOntoCoordPlane(indices[1], near, info);
    }


    // Check if outside. Optimisation: could do some checks on distance already
    // on components above
    if (info.point().distSqr(sample) > nearestDistSqr)
    {
        info.setMiss();
        info.setIndex(-1);
    }

    return info;
}


Foam::pointIndexHit Foam::searchableBox::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    NotImplemented;
    return pointIndexHit();
}


void Foam::searchableBox::findLine_impl
(
    const point& start,
    const point& end,
    pointIndexHit& info
) const
{
    info = pointIndexHit(false, start, -1);

    bool foundInter = false;

    if (posBits(start) == 0)
    {
        if (posBits(end) == 0)
        {
            // Both start and end inside.
            foundInter = false;
        }
        else
        {
            // end is outside. Clip to bounding box.
            foundInter = intersects(end, start, info.point());
        }
    }
    else
    {
        // start is outside. Clip to bounding box.
        foundInter = intersects(start, end, info.point());
    }


    // Classify point
    if (foundInter)
    {
        info.setHit();

        for (direction dir = 0; dir < vector::nComponents; ++dir)
        {
            if (info.point()[dir] == min()[dir])
            {
                info.setIndex(2*dir);
                break;
            }
            else if (info.point()[dir] == max()[dir])
            {
                info.setIndex(2*dir+1);
                break;
            }
        }

        if (info.index() == -1)
        {
            FatalErrorInFunction
                << "point " << info.point()
                << " on segment " << start << end
                << " should be on face of " << *this
                << " but it isn't." << abort(FatalError);
        }
    }
}


Foam::pointIndexHit Foam::searchableBox::findLine
(
    const point& start,
    const point& end
) const
{
    pointIndexHit info;
    findLine_impl(start, end, info);
    return info;
}


Foam::pointIndexHit Foam::searchableBox::findLineAny
(
    const point& start,
    const point& end
) const
{
    pointIndexHit info;
    findLine_impl(start, end, info);
    return info;
}


void Foam::searchableBox::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.resize(samples.size());

    const point bbMid(centre());

    scalar maxDistSqr = searchableSurface::defaultMaxDistSqr();

    forAll(samples, i)
    {
        if (i < nearestDistSqr.size()) { maxDistSqr = nearestDistSqr[i]; }

        findNearest_impl(bbMid, samples[i], maxDistSqr, info[i]);
    }
}


void Foam::searchableBox::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize_nocopy(start.size());

    forAll(start, i)
    {
        findLine_impl(start[i], end[i], info[i]);
    }
}


void Foam::searchableBox::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize_nocopy(start.size());

    forAll(start, i)
    {
        findLine_impl(start[i], end[i], info[i]);
    }
}


void Foam::searchableBox::findLineAll
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
        pointIndexHit inter;
        findLine_impl(startp, endp, inter);

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
                findLine_impl(pt, endp, inter);

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


void Foam::searchableBox::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.resize_fill(info.size(), 0);
}


void Foam::searchableBox::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.resize_nocopy(info.size());

    forAll(info, i)
    {
        if
        (
            label facei = info[i].index();
            (info[i].hit() && (facei >= 0 && facei < treeBoundBox::nFaces()))
        )
        {
            normal[i] = treeBoundBox::faceNormals[facei];
        }
        else
        {
            // Set to what?
            normal[i] = Foam::zero{};
        }
    }
}


Foam::label Foam::searchableBox::measure
(
    List<vector>* centres,
    List<scalar>* boundingRadiusSqr,
    List<label>* regions,
    List<vector>* normals
) const
{
    const label output_size = treeBoundBox::nFaces();

    bool usedParams = false;

    // Extract face centres
    if (centres)
    {
        usedParams = true;
        centres->resize_nocopy(output_size); // MUST resize first
        treeBoundBox::faceCentres(*centres);
    }

    // Extract bounding sphere radius squared
    // Box faces are quads, points are axis-aligned,
    // so the face centre is trivial and the only need one face corner
    // to determine the bounding sphere radius

    if (boundingRadiusSqr)
    {
        usedParams = true;
        auto& radiusSqr = *boundingRadiusSqr;
        radiusSqr.resize_nocopy(output_size);

        for (label facei = 0; facei < treeBoundBox::nFaces(); ++facei)
        {
            label pointi = treeBoundBox::octFaces()[facei][0];
            point corner = treeBoundBox::octCorner(pointi);
            point fc = treeBoundBox::faceCentre(facei);

            // Add a bit to make sure all points are tested inside
            radiusSqr[facei] = (fc.distSqr(corner) + Foam::sqr(SMALL));
        }
    }

    // Extract regions
    if (regions)
    {
        usedParams = true;
        regions->resize_fill(output_size, 0);
    }

    // Extract normals
    if (normals)
    {
        usedParams = true;
        auto& normal = *normals;
        normal.resize_nocopy(treeBoundBox::nFaces());

        for (label facei = 0; facei < treeBoundBox::nFaces(); ++facei)
        {
            normal[facei] = treeBoundBox::faceNormals[facei];
        }
    }

    return (usedParams ? output_size : 0);
}


void Foam::searchableBox::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.resize_nocopy(points.size());

    forAll(points, pointi)
    {
        volType[pointi] =
        (
            treeBoundBox::contains(points[pointi])
          ? volumeType::INSIDE
          : volumeType::OUTSIDE
        );
    }
}


// ************************************************************************* //
