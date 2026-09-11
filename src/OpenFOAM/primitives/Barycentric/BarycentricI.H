/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    ANB WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline constexpr Foam::Barycentric<Cmpt>::Barycentric(Foam::zero)
:
    Barycentric::vsType(Foam::zero{})
{}


template<class Cmpt>
inline constexpr Foam::Barycentric<Cmpt>::Barycentric
(
    const Cmpt& va,
    const Cmpt& vb,
    const Cmpt& vc,
    const Cmpt& vd
)
:
    // aggregate initialization
    Barycentric::vsType(Foam::direct_construct_t{}, va, vb, vc, vd)
{}


template<class Cmpt>
inline Foam::Barycentric<Cmpt>::Barycentric
(
    const Cmpt& va,
    const Cmpt& vb,
    const Cmpt& vc
)
:
    Barycentric<Cmpt>(va, vb, vc, (pTraits<Cmpt>::one - va - vb - vc))
{}


// * * * * * * * * * * * * * * * Member Operations * * * * * * * * * * * * * //

template<class Cmpt>
inline Cmpt Foam::Barycentric<Cmpt>::inner(const Barycentric<Cmpt>& b2) const
{
    const Barycentric<Cmpt>& b1 = *this;

    return (b1.a()*b2.a() + b1.b()*b2.b() + b1.c()*b2.c() + b1.d()*b2.d());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Cmpt operator&(const Barycentric<Cmpt>& b1, const Barycentric<Cmpt>& b2)
{
    return b1.inner(b2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
