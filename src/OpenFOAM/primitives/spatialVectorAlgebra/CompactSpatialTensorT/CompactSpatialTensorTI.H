/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
inline constexpr Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT
(
    Foam::zero
)
:
    CompactSpatialTensorT::msType(Foam::zero{})
{}


template<class Cmpt>
inline Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT
(
    const typename CompactSpatialTensorT::msType& ms
)
:
    CompactSpatialTensorT::msType(ms)
{}


template<class Cmpt>
inline constexpr Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT
(
    const Cmpt& t00, const Cmpt& t01, const Cmpt& t02,
    const Cmpt& t10, const Cmpt& t11, const Cmpt& t12,
    const Cmpt& t20, const Cmpt& t21, const Cmpt& t22,
    const Cmpt& t30, const Cmpt& t31, const Cmpt& t32,
    const Cmpt& t40, const Cmpt& t41, const Cmpt& t42,
    const Cmpt& t50, const Cmpt& t51, const Cmpt& t52
)
:
    // aggregate initialization
    CompactSpatialTensorT::msType
    (
        Foam::direct_construct_t{},
        t00, t01, t02,
        t10, t11, t12,
        t20, t21, t22,
        t30, t31, t32,
        t40, t41, t42,
        t50, t51, t52
    )
{}


template<class Cmpt>
inline Foam::CompactSpatialTensorT<Cmpt>::CompactSpatialTensorT(Istream& is)
:
    CompactSpatialTensorT::msType(is)
{}


// ************************************************************************* //
