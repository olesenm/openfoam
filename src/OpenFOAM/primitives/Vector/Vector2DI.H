/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline constexpr Foam::Vector2D<Cmpt>::Vector2D(Foam::zero)
:
    Vector2D::vsType(Foam::zero{})
{}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>::Vector2D
(
    const VectorSpace<Vector2D<Cmpt>, Cmpt, 2>& vs
)
:
    Vector2D::vsType(vs)
{}


template<class Cmpt>
inline constexpr Foam::Vector2D<Cmpt>::Vector2D
(
    const Cmpt& vx,
    const Cmpt& vy
)
:
    // aggregate initialization
    Vector2D::vsType(Foam::direct_construct_t{}, vx, vy)
{}


template<class Cmpt>
inline constexpr Foam::Vector2D<Cmpt>::Vector2D
(
    const Vector<Cmpt>& vs
)
:
    // aggregate initialization
    Vector2D::vsType(Foam::direct_construct_t{}, vs.x(), vs.y())
{}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>::Vector2D(Istream& is)
:
    Vector2D::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline constexpr Foam::scalar Foam::Vector2D<Cmpt>::magSqr() const
{
    return
    (
        Foam::magSqr(this->x())
      + Foam::magSqr(this->y())
    );
}


template<class Cmpt>
inline Foam::scalar Foam::Vector2D<Cmpt>::mag() const
{
    return std::sqrt(this->magSqr());
}


template<class Cmpt>
inline constexpr Foam::scalar
Foam::Vector2D<Cmpt>::distSqr(const Vector2D<Cmpt>& v2) const
{
    return
    (
        Foam::magSqr(v2.x() - this->x())
      + Foam::magSqr(v2.y() - this->y())
    );
}


template<class Cmpt>
inline Foam::scalar Foam::Vector2D<Cmpt>::dist(const Vector2D<Cmpt>& v2) const
{
    return std::sqrt(this->distSqr(v2));
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>& Foam::Vector2D<Cmpt>::normalise(const scalar tol)
{
    const scalar s = this->mag();

    if (s < tol)
    {
        *this = Zero;
    }
    else
    {
        *this /= s;
    }

    return *this;
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>&
Foam::Vector2D<Cmpt>::removeCollinear(const Vector2D<Cmpt>& unitVec)
{
    *this -= (*this & unitVec) * unitVec;
    return *this;
}


// * * * * * * * * * * * * * Comparison Operations * * * * * * * * * * * * * //

template<class Cmpt>
inline bool
Foam::Vector2D<Cmpt>::less_xy(const Vector2D<Cmpt>& a, const Vector2D<Cmpt>& b)
{
    return
    (
        (a.x() < b.x())
     || (!(b.x() < a.x()) && (a.y() < b.y()))
    );
}


template<class Cmpt>
inline bool
Foam::Vector2D<Cmpt>::less_yx(const Vector2D<Cmpt>& a, const Vector2D<Cmpt>& b)
{
    return
    (
        (a.y() < b.y())
     || (!(b.y() < a.y()) && (a.x() < b.x()))
    );
}


template<class Cmpt>
inline constexpr Cmpt
Foam::Vector2D<Cmpt>::inner(const Vector2D<Cmpt>& b) const
{
    return (x()*b.x() + y()*b.y());
}


template<class Cmpt>
inline constexpr Cmpt
Foam::Vector2D<Cmpt>::perp(const Vector2D<Cmpt>& b) const
{
    return x()*b.y() - y()*b.x();
}


template<class Cmpt>
inline bool Foam::Vector2D<Cmpt>::isClose
(
    const Vector2D<Cmpt>& b,
    const scalar tol
) const
{
    return (Foam::mag(x() - b.x()) < tol && Foam::mag(y() - b.y()) < tol);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Linear interpolation of vectors a and b by factor t
template<class Cmpt>
inline constexpr Vector2D<Cmpt> lerp
(
    const Vector2D<Cmpt>& a,
    const Vector2D<Cmpt>& b,
    const scalar t
)
{
    const scalar onet = (1-t);

    return Vector2D<Cmpt>
    (
        onet*a.x() + t*b.x(),
        onet*a.y() + t*b.y()
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline constexpr Cmpt
operator&(const Vector2D<Cmpt>& v1, const Vector2D<Cmpt>& v2)
{
    return v1.inner(v2);
}


template<class Cmpt>
inline constexpr Vector2D<Cmpt>
operator*(const Cmpt& s, const Vector2D<Cmpt>& v)
{
    return Vector2D<Cmpt>(s*v.x(), s*v.y());
}


template<class Cmpt>
inline constexpr Vector2D<Cmpt>
operator*(const Vector2D<Cmpt>& v, const Cmpt& s)
{
    return s*v;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
