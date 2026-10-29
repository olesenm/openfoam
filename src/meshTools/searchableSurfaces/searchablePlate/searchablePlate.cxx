/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "searchablePlate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchablePlate, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchablePlate,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchablePlate,
        dict,
        plate
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::searchablePlate::calcNormal(const point& span)
{
    direction normalDir = vector::nComponents;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        if (span[dir] < 0)
        {
            // Negative entry. Flag and exit.
            normalDir = vector::nComponents;
            break;
        }
        else if (span[dir] < VSMALL)
        {
            if (normalDir != vector::nComponents)
            {
                // Multiple zero entries. Flag and exit.
                normalDir = vector::nComponents;
                break;
            }

            normalDir = dir;
        }
    }

    if (normalDir == vector::nComponents)
    {
        FatalErrorInFunction
            << "Span should have two positive and one zero entry: "
            << span << nl
            << exit(FatalError);
    }

    return normalDir;
}


// Returns miss or hit with face (always 0)
void Foam::searchablePlate::findNearest_impl
(
    const point& sample,
    const scalar nearestDistSqr,
    pointIndexHit& info
) const
{
    // For every component direction can be
    // left of min, right of max or inbetween.
    // - outside points: project first one x plane (either min().x()
    // or max().x()), then onto y plane and finally z. You should be left
    // with intersection point
    // - inside point: find nearest side (compare to mid point). Project onto
    //   that.

    // Project point on plane.
    info.hitPoint(sample, 0);
    info.point()[normalDir_] = origin_[normalDir_];

    // Clip to edges if outside
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        if (dir != normalDir_)
        {
            if (info.point()[dir] < origin_[dir])
            {
                info.point()[dir] = origin_[dir];
            }
            else if (info.point()[dir] > origin_[dir]+span_[dir])
            {
                info.point()[dir] = origin_[dir]+span_[dir];
            }
        }
    }

    // Check if outside. Optimisation: could do some checks on distance already
    // on components above
    if (nearestDistSqr < info.point().distSqr(sample))
    {
        info.setMiss();
        info.setIndex(-1);
    }
}


void Foam::searchablePlate::findLine_impl
(
    const point& start,
    const point& end,
    pointIndexHit& info
) const
{
    info.hitPoint(point::zero, 0);

    const vector dir(end-start);

    if (mag(dir[normalDir_]) < VSMALL)
    {
        info.setMiss();
        info.setIndex(-1);
        return;
    }

    const scalar t = (origin_[normalDir_]-start[normalDir_]) / dir[normalDir_];

    if (t < 0 || t > 1)
    {
        info.setMiss();
        info.setIndex(-1);
        return;
    }


    // hit (probably)
    info.point() = start+t*dir;
    info.point()[normalDir_] = origin_[normalDir_];

    // Clip to edges
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        if (dir != normalDir_)
        {
            if (info.point()[dir] < origin_[dir])
            {
                info.setMiss();
                info.setIndex(-1);
                break;
            }
            else if (info.point()[dir] > origin_[dir]+span_[dir])
            {
                info.setMiss();
                info.setIndex(-1);
                break;
            }
        }
    }

    // Debug
    if (info.hit())
    {
        treeBoundBox bb(origin_, origin_+span_);
        bb.min()[normalDir_] -= 1e-6;
        bb.max()[normalDir_] += 1e-6;

        if (!bb.contains(info.point()))
        {
            FatalErrorInFunction
                << "bb:" << bb << endl
                << "origin_:" << origin_ << endl
                << "span_:" << span_ << endl
                << "normalDir_:" << normalDir_ << endl
                << "hitPoint:" << info.point()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchablePlate::searchablePlate
(
    const IOobject& io,
    const point& origin,
    const vector& span
)
:
    searchableSurface(io),
    origin_(origin),
    span_(span),
    normalDir_(calcNormal(span_))
{
    DebugInFunction
        << " origin:" << origin_
        << " origin+span:" << origin_+span_
        << " normal:" << vector::componentNames[normalDir_] << nl;

    bounds() = boundBox(origin_, origin_ + span_);
}


Foam::searchablePlate::searchablePlate
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchablePlate(io, dict.get<point>("origin"), dict.get<vector>("span"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchablePlate::regions() const
{
    if (regions_.empty())
    {
        regions_.resize(1);
        regions_.first() = "region0";
    }
    return regions_;
}


Foam::tmp<Foam::pointField> Foam::searchablePlate::coordinates() const
{
    return tmp<pointField>::New(1, (origin_ + 0.5*span_));
}


void Foam::searchablePlate::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    // Add a bit to make sure all points are tested inside
    centres.resize_fill(1, (origin_ + 0.5*span_));
    radiusSqr.resize_fill(1, (Foam::magSqr(0.5*span_) + Foam::sqr(SMALL)));
}


Foam::tmp<Foam::pointField> Foam::searchablePlate::points() const
{
    auto tpts = tmp<pointField>::New(4, origin_);
    auto& pts = tpts.ref();

    pts[2] += span_;

    if (span_.x() < span_.y() && span_.x() < span_.z())
    {
        pts[1].y() += span_.y();
        pts[3].z() += span_.z();
    }
    else if (span_.y() < span_.z())
    {
        pts[1].x() += span_.x();
        pts[3].z() += span_.z();
    }
    else
    {
        pts[1].x() += span_.x();
        pts[3].y() += span_.y();
    }

    return tpts;
}


bool Foam::searchablePlate::overlaps(const boundBox& bb) const
{
    return bb.overlaps(bounds());
}


void Foam::searchablePlate::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.resize_nocopy(samples.size());

    forAll(samples, i)
    {
        findNearest_impl(samples[i], nearestDistSqr[i], info[i]);
    }
}


void Foam::searchablePlate::findLine
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


void Foam::searchablePlate::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    findLine(start, end, info);
}


void Foam::searchablePlate::findLineAll
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


void Foam::searchablePlate::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.resize_fill(info.size(), 0);
}


void Foam::searchablePlate::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.resize_fill(info.size(), Foam::zero{});

    for (auto& n : normal)
    {
        n[normalDir_] = 1.0;
    }
}


void Foam::searchablePlate::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorInFunction
        << "Volume type not supported for plate."
        << exit(FatalError);
}


// ************************************************************************* //
