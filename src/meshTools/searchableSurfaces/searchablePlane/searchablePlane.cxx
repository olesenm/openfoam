/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "searchablePlane.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchablePlane, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchablePlane,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchablePlane,
        dict,
        plane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchablePlane::findLine_impl
(
    const point& start,
    const point& end,
    pointIndexHit& info
) const
{
    info.hitPoint(point::zero, 0);

    linePointRef l(start, end);

    const scalar t = plane::lineIntersect(l);

    if (t < 0 || t > 1)
    {
        info.setMiss();
        info.setIndex(-1);
    }
    else
    {
        info.point() = (start+t*l.vec());
    }
}


Foam::boundBox Foam::searchablePlane::calcBounds() const
{
    boundBox bb(boundBox::greatBox);

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (mag(normal()[dir]) - 1 < SMALL)
        {
            bb.min()[dir] = 0;
            bb.max()[dir] = 0;
            break;
        }
    }

    return bb;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchablePlane::searchablePlane
(
    const IOobject& io,
    const point& basePoint,
    const vector& normal
)
:
    searchableSurface(io),
    plane(basePoint, normal)
{
    bounds() = calcBounds();
}


Foam::searchablePlane::searchablePlane
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    plane(dict)
{
    bounds() = calcBounds();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchablePlane::regions() const
{
    if (regions_.empty())
    {
        regions_.resize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchablePlane::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.resize_fill(1, origin());
    radiusSqr.resize_fill(1, Foam::sqr(GREAT));
}


void Foam::searchablePlane::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.resize(samples.size());

    scalar maxDistSqr = searchableSurface::defaultMaxDistSqr();

    forAll(samples, i)
    {
        if (i < nearestDistSqr.size()) { maxDistSqr = nearestDistSqr[i]; }

        info[i].setPoint(plane::nearestPoint(samples[i]));

        if (info[i].point().distSqr(samples[i]) > maxDistSqr)
        {
            info[i].setIndex(-1);
            info[i].setMiss();
        }
        else
        {
            info[i].setIndex(0);
            info[i].setHit();
        }
    }
}


void Foam::searchablePlane::findLine
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


void Foam::searchablePlane::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    findLine(start, end, info);
}


void Foam::searchablePlane::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.resize_nocopy(start.size());

    forAll(start, i)
    {
        pointIndexHit near;
        findLine_impl(start[i], end[i], near);

        if (near.hit())
        {
            info[i].resize_fill(1, near);
        }
        else
        {
            info[i].clear();
        }
    }
}


void Foam::searchablePlane::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.resize_fill(info.size(), 0);
}


void Foam::searchablePlane::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& n
) const
{
    n.resize_fill(info.size(), plane::normal());
}


Foam::label Foam::searchablePlane::measure
(
    List<vector>* centres,
    List<scalar>* radiusSqr,
    List<label>* regions,
    List<vector>* normals
) const
{
    const label output_size = 1;

    bool usedParams = false;

    if (centres)
    {
        usedParams = true;
        centres->resize_fill(output_size, origin());
    }

    if (radiusSqr)
    {
        usedParams = true;
        radiusSqr->resize_fill(output_size, Foam::sqr(GREAT));
    }

    if (regions)
    {
        usedParams = true;
        regions->resize_fill(output_size, 0);
    }

    if (normals)
    {
        usedParams = true;
        normals->resize_fill(output_size, plane::normal());
    }

    return (usedParams ? output_size : 0);
}


void Foam::searchablePlane::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorInFunction
        << "Volume type not supported for plane."
        << exit(FatalError);
}


// ************************************************************************* //
