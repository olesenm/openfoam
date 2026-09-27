/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "boundBox.H"

// * * * * * * * * * * * * * Geometrical Information * * * * * * * * * * * * //

// Box corners as per hex cellmodel

template<Foam::direction I>
inline Foam::point Foam::boundBox::hexCorner() const
{
    /**/ if constexpr (0 == I)  // == octCorner<0>()
    {
        return min_;
    }
    else if constexpr (1 == I)  // == octCorner<1>()
    {
        return point(max_.x(), min_.y(), min_.z());
    }
    else if constexpr (2 == I)  // == octCorner<3>()
    {
        return point(max_.x(), max_.y(), min_.z());
    }
    else if constexpr (3 == I)  // == octCorner<2>()
    {
        return point(min_.x(), max_.y(), min_.z());
    }
    else if constexpr (4 == I)  // == octCorner<4>()
    {
        return point(min_.x(), min_.y(), max_.z());
    }
    else if constexpr (5 == I)  // == octCorner<5>()
    {
        return point(max_.x(), min_.y(), max_.z());
    }
    else if constexpr (6 == I)  // == octCorner<7>()
    {
        return max_;
    }
    else if constexpr (7 == I)  // == octCorner<6>()
    {
        return point(min_.x(), max_.y(), max_.z());
    }
    else
    {
        static_assert(I < 8, "Corner index [0..7]");
        return point::zero;
    }
}


// Box corners by octant encoding

template<Foam::direction Octant>
inline Foam::point Foam::boundBox::octCorner() const
{
    static_assert(Octant < 8, "Corner index [0..7]");
    return point
    (
        (Octant & 1) ? max_.x() : min_.x(),
        (Octant & 2) ? max_.y() : min_.y(),
        (Octant & 4) ? max_.z() : min_.z()
    );
}


inline Foam::point Foam::boundBox::octCorner(direction octant) const
{
    return point
    (
        (octant & 1) ? max_.x() : min_.x(),
        (octant & 2) ? max_.y() : min_.y(),
        (octant & 4) ? max_.z() : min_.z()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::boundBox::boundBox()
:
    min_(invertedBox.min()),
    max_(invertedBox.max())
{}


inline Foam::boundBox::boundBox(Foam::zero_one)
:
    min_(point::zero),
    max_(point::one)
{}


inline Foam::boundBox::boundBox(const point& p)
:
    min_(p),
    max_(p)
{}


inline Foam::boundBox::boundBox(const point& min, const point& max)
:
    min_(min),
    max_(max)
{}


inline Foam::boundBox::boundBox(const Pair<point>& bb)
:
    min_(bb.first()),
    max_(bb.second())
{}


inline Foam::boundBox::boundBox(const boundBox& a, const boundBox& b)
:
    min_(::Foam::min(a.min(), b.min())),
    max_(::Foam::max(a.max(), b.max()))
{}


inline Foam::boundBox::boundBox(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::boundBox::empty() const noexcept
{
    // Is empty/invalid if any component has (max < min), ie, !(min <= max)
    return
    (
        (max_.x() < min_.x())
     || (max_.y() < min_.y())
     || (max_.z() < min_.z())
    );
}


inline bool Foam::boundBox::good() const noexcept
{
    return !empty();
}


inline const Foam::point& Foam::boundBox::min() const noexcept
{
    return min_;
}


inline const Foam::point& Foam::boundBox::max() const noexcept
{
    return max_;
}


inline Foam::point& Foam::boundBox::min() noexcept
{
    return min_;
}


inline Foam::point& Foam::boundBox::max() noexcept
{
    return max_;
}


inline Foam::point Foam::boundBox::centre() const
{
    return point
    (
        (0.5*min_.x() + 0.5*max_.x()),
        (0.5*min_.y() + 0.5*max_.y()),
        (0.5*min_.z() + 0.5*max_.z())
    );
}


inline Foam::scalar Foam::boundBox::centre(direction cmpt) const
{
    return (0.5*min_[cmpt] + 0.5*max_[cmpt]);
}


inline Foam::point Foam::boundBox::lerp(scalar t) const
{
    return point
    (
        Foam::lerp(min_.x(), max_.x(), t),
        Foam::lerp(min_.y(), max_.y(), t),
        Foam::lerp(min_.z(), max_.z(), t)
    );
}


inline Foam::point Foam::boundBox::lerp(const vector& t) const
{
    return point
    (
        Foam::lerp(min_.x(), max_.x(), t.x()),
        Foam::lerp(min_.y(), max_.y(), t.y()),
        Foam::lerp(min_.z(), max_.z(), t.z())
    );
}


inline Foam::vector Foam::boundBox::span() const
{
    return (max_ - min_);
}


inline Foam::scalar Foam::boundBox::span(direction cmpt) const
{
    return (max_[cmpt] - min_[cmpt]);
}


inline Foam::scalar Foam::boundBox::mag() const
{
    return min_.dist(max_);
}


inline Foam::scalar Foam::boundBox::magSqr() const
{
    return min_.distSqr(max_);
}


inline Foam::scalar Foam::boundBox::volume() const
{
    return cmptProduct(span());
}


inline Foam::scalar Foam::boundBox::minDim() const
{
    return cmptMin(span());
}


inline Foam::scalar Foam::boundBox::maxDim() const
{
    return cmptMax(span());
}


inline Foam::scalar Foam::boundBox::avgDim() const
{
    return cmptAv(span());
}


inline Foam::direction Foam::boundBox::minDir() const
{
    direction cmpt = 0;

    scalar best = ROOTVGREAT;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar dist = (max_[dir] - min_[dir]);
        if (dist < best && dist > 0)
        {
            best = dist;
            cmpt = dir;
        }
    }

    return cmpt;
}


inline Foam::direction Foam::boundBox::maxDir() const
{
    direction cmpt = 0;

    scalar best = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar dist = (max_[dir] - min_[dir]);
        if (dist > best)
        {
            best = dist;
            cmpt = dir;
        }
    }

    return cmpt;
}


inline Foam::Pair<Foam::scalar> Foam::boundBox::distSqr(const point& p) const
{
    scalar sumMinSqr = 0;
    scalar sumMaxSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        scalar a = Foam::magSqr(min_[dir] - p[dir]);
        scalar b = Foam::magSqr(max_[dir] - p[dir]);

        sumMinSqr += std::min(a, b);
        sumMaxSqr += std::max(a, b);
    }

    return { sumMinSqr, sumMaxSqr };
}


inline Foam::scalar Foam::boundBox::minDistSqr(const point& p) const
{
    scalar sumMinSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        scalar a = Foam::magSqr(min_[dir] - p[dir]);
        scalar b = Foam::magSqr(max_[dir] - p[dir]);

        sumMinSqr += std::min(a, b);
    }

    return sumMinSqr;
}


inline Foam::scalar Foam::boundBox::maxDistSqr(const point& p) const
{
    scalar sumMaxSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        scalar a = Foam::magSqr(min_[dir] - p[dir]);
        scalar b = Foam::magSqr(max_[dir] - p[dir]);

        sumMaxSqr += std::max(a, b);
    }

    return sumMaxSqr;
}


Foam::scalar Foam::boundBox::minDist(const point& p) const
{
    return std::sqrt(minDistSqr(p));
}


Foam::scalar Foam::boundBox::maxDist(const point& p) const
{
    return std::sqrt(maxDistSqr(p));
}


inline int Foam::boundBox::nDim() const
{
    int ngood = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar dist = (max_[dir] - min_[dir]);
        if (dist < 0)
        {
            return -1;
        }
        else if (dist > 0)
        {
            ++ngood;
        }
    }

    return ngood;
}


inline void Foam::boundBox::reset()
{
    min_ = invertedBox.min();
    max_ = invertedBox.max();
}


inline void Foam::boundBox::reset(Foam::zero_one)
{
    min_ = point::zero;
    max_ = point::one;
}


inline void Foam::boundBox::reset(const point& pt)
{
    min_ = pt;
    max_ = pt;
}


inline void Foam::boundBox::reset(const point& min, const point& max)
{
    min_ = min;
    max_ = max;
}


inline void Foam::boundBox::add(const boundBox& bb)
{
    min_ = ::Foam::min(min_, bb.min_);
    max_ = ::Foam::max(max_, bb.max_);
}


inline void Foam::boundBox::add(const Pair<point>& bb)
{
    min_ = ::Foam::min(min_, bb.first());
    max_ = ::Foam::max(max_, bb.second());
}


inline void Foam::boundBox::add(const point& pt)
{
    min_ = ::Foam::min(min_, pt);
    max_ = ::Foam::max(max_, pt);
}


inline void Foam::boundBox::add(const point& pt0, const point& pt1)
{
    add(pt0);
    add(pt1);
}


inline void Foam::boundBox::add(const UList<point>& points)
{
    for (const point& p : points)
    {
        add(p);
    }
}


inline void Foam::boundBox::add(const tmp<pointField>& tpoints)
{
    add(tpoints());
    tpoints.clear();
}


inline void Foam::boundBox::grow(scalar delta)
{
    min_.x() -= delta; min_.y() -= delta; min_.z() -= delta;
    max_.x() += delta; max_.y() += delta; max_.z() += delta;
}


inline void Foam::boundBox::grow(const vector& delta)
{
    min_ -= delta;
    max_ += delta;
}


inline void Foam::boundBox::inflate(scalar factor)
{
    grow(factor*mag());
}


inline void Foam::boundBox::inflate_components(scalar epsilon)
{
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar middle = (0.5*min_[dir] + 0.5*max_[dir]);
        const scalar half_span = (1 + epsilon)*(max_[dir] - middle);

        if (half_span < 0)
        {
            // Clamp away tiny negative half-spans (rounding noise) to zero.
            // Ignore non-trivial inverted spans completely.

            if (Foam::mag(half_span) < ROOTVSMALL)
            {
                min_[dir] = middle;
                max_[dir] = middle;
            }
        }
        else
        {
            min_[dir] = (middle - half_span);
            max_[dir] = (middle + half_span);
        }
    }
}


inline bool Foam::boundBox::box_box_overlaps
(
    const point& minA, const point& maxA,   // boxA
    const point& minB, const point& maxB    // boxB
)
{
    return
    (
        minA.x() <= maxB.x() && minB.x() <= maxA.x()
     && minA.y() <= maxB.y() && minB.y() <= maxA.y()
     && minA.z() <= maxB.z() && minB.z() <= maxA.z()
    );
}


inline bool Foam::boundBox::box_sphere_overlaps
(
    const point& corner0,
    const point& corner1,
    const point& centre,
    const scalar radiusSqr
)
{
    // Find out where centre is in relation to bb.
    // Find nearest point on bb.
    scalar distSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar d0 = corner0[dir] - centre[dir];
        const scalar d1 = corner1[dir] - centre[dir];

        if ((d0 > 0) != (d1 > 0))
        {
            // centre inside both extrema. This component does not add any
            // distance.
        }
        else
        {
            distSqr += std::min(Foam::magSqr(d0), Foam::magSqr(d1));

            if (distSqr > radiusSqr)
            {
                return false;
            }
        }
    }

    return true;
}


inline bool Foam::boundBox::box_sphere_contains
(
    const point& min_corner,    // box corner (min)
    const point& max_corner,    // box corner (max)
    const point& centre,
    const scalar radiusSqr      // Input (squared radius, avoids sqrt)
)
{
    // The box contains the sphere iff, on every axis, the margin between
    // the sphere centre and both box faces is >= radius.
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const auto margin = std::min
        (
            (centre[dir] - min_corner[dir]),
            (max_corner[dir] - centre[dir])
        );

        // Need a sufficient (and non-negative) margin on each component
        if (margin < 0 || Foam::sqr(margin) < radiusSqr)
        {
            return false;
        }
    }

    return true;
}


inline bool Foam::boundBox::intersects(const linePointRef& ln) const
{
    return intersects(ln.start(), ln.end());
}


inline bool Foam::boundBox::overlaps(const boundBox& bb) const
{
    return box_box_overlaps(min_, max_, bb.min(), bb.max());
}


inline bool Foam::boundBox::overlaps
(
    const point& centre,
    scalar radiusSqr
) const
{
    return box_sphere_overlaps(min_, max_, centre, radiusSqr);
}


inline bool Foam::boundBox::contains(const point& pt) const noexcept
{
    return
    (
        min_.x() <= pt.x() && pt.x() <= max_.x()
     && min_.y() <= pt.y() && pt.y() <= max_.y()
     && min_.z() <= pt.z() && pt.z() <= max_.z()
    );
}


inline bool Foam::boundBox::containsInside(const point& pt) const noexcept
{
    return
    (
        min_.x() < pt.x() && pt.x() < max_.x()
     && min_.y() < pt.y() && pt.y() < max_.y()
     && min_.z() < pt.z() && pt.z() < max_.z()
    );
}


inline bool Foam::boundBox::contains(const boundBox& bb) const noexcept
{
    return contains_all(bb.min(), bb.max());
}


inline bool Foam::boundBox::contains
(
    const point& centre,
    scalar radiusSqr
) const
{
    return box_sphere_contains(min_, max_, centre, radiusSqr);
}


inline bool Foam::boundBox::contains_all
(
    const point& p0,
    const point& p1
) const
{
    return (contains(p0) && contains(p1));
}


inline bool Foam::boundBox::contains_any
(
    const point& p0,
    const point& p1
) const
{
    return (contains(p0) || contains(p1));
}


inline bool Foam::boundBox::contains_all
(
    const point& p0,
    const point& p1,
    const point& p2
) const
{
    return (contains(p0) && contains(p1) && contains(p2));
}


inline bool Foam::boundBox::contains_any
(
    const point& p0,
    const point& p1,
    const point& p2
) const
{
    return (contains(p0) || contains(p1) || contains(p2));
}


inline bool Foam::boundBox::contains_all(const linePointRef& ln) const
{
    return contains_all(ln.start(), ln.end());
}


inline bool Foam::boundBox::contains_any(const linePointRef& ln) const
{
    return contains_any(ln.start(), ln.end());
}


inline Foam::boundBox
Foam::boundBox::returnReduce(int communicator) const
{
    boundBox work(*this);
    work.reduce(communicator);
    return work;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
