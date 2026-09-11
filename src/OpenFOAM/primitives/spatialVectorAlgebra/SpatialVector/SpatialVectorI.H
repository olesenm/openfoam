/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
inline constexpr Foam::SpatialVector<Cmpt>::SpatialVector(Foam::zero)
:
    SpatialVector::vsType(Foam::zero{})
{}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector
(
    const typename SpatialVector::vsType& vs
)
:
    SpatialVector::vsType(vs)
{}


template<class Cmpt>
inline constexpr Foam::SpatialVector<Cmpt>::SpatialVector
(
    const Vector<Cmpt>& w,
    const Vector<Cmpt>& l
)
:
    // aggregate initialization
    SpatialVector::vsType
    (
        Foam::direct_construct_t{},
        w.x(), w.y(), w.z(),
        l.x(), l.y(), l.z()
    )
{}


template<class Cmpt>
inline constexpr Foam::SpatialVector<Cmpt>::SpatialVector
(
    const Cmpt& v0,
    const Cmpt& v1,
    const Cmpt& v2,
    const Cmpt& v3,
    const Cmpt& v4,
    const Cmpt& v5
)
:
    // aggregate initialization
    SpatialVector::vsType
    (
        Foam::direct_construct_t{},
        v0, v1, v2,
        v3, v4, v5
    )
{}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector(Istream& is)
:
    SpatialVector::vsType(is)
{}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::dual::dual(const SpatialVector& v)
:
    v_(v)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline constexpr Foam::Vector<Cmpt> Foam::SpatialVector<Cmpt>::w() const
{
    return Vector<Cmpt>(this->v_[WX], this->v_[WY], this->v_[WZ]);
}

template<class Cmpt>
inline constexpr Foam::Vector<Cmpt> Foam::SpatialVector<Cmpt>::l() const
{
    return Vector<Cmpt>(this->v_[LX], this->v_[LY], this->v_[LZ]);
}


template<class Cmpt>
const Foam::SpatialVector<Cmpt>& Foam::SpatialVector<Cmpt>::dual::v() const
{
    return v_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline typename Foam::SpatialVector<Cmpt>::dual
Foam::SpatialVector<Cmpt>::operator*() const
{
    return dual(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Return the cross-product between two spatial vectors
template<class Cmpt>
inline SpatialVector<Cmpt> operator^
(
    const SpatialVector<Cmpt>& u,
    const SpatialVector<Cmpt>& v
)
{
    return SpatialVector<Cmpt>
    (
       -u.wz()*v.wy() + u.wy()*v.wz(),
        u.wz()*v.wx() - u.wx()*v.wz(),
       -u.wy()*v.wx() + u.wx()*v.wy(),
       -u.lz()*v.wy() + u.ly()*v.wz() - u.wz()*v.ly() + u.wy()*v.lz(),
        u.lz()*v.wx() - u.lx()*v.wz() + u.wz()*v.lx() - u.wx()*v.lz(),
       -u.ly()*v.wx() + u.lx()*v.wy() - u.wy()*v.lx() + u.wx()*v.ly()
    );
}


//- Return the dual cross-product between two spatial vectors
template<class Cmpt>
inline SpatialVector<Cmpt> operator^
(
    const SpatialVector<Cmpt>& v,
    const typename SpatialVector<Cmpt>::dual& df
)
{
    const SpatialVector<Cmpt>& f = df.v();

    return SpatialVector<Cmpt>
    (
       -v.wz()*f.wy() + v.wy()*f.wz() - v.lz()*f.ly() + v.ly()*f.lz(),
        v.wz()*f.wx() - v.wx()*f.wz() + v.lz()*f.lx() - v.lx()*f.lz(),
       -v.wy()*f.wx() + v.wx()*f.wy() - v.ly()*f.lx() + v.lx()*f.ly(),
       -v.wz()*f.ly() + v.wy()*f.lz(),
        v.wz()*f.lx() - v.wx()*f.lz(),
       -v.wy()*f.lx() + v.wx()*f.ly()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
