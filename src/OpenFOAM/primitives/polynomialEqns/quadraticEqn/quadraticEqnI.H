/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

inline constexpr Foam::quadraticEqn::quadraticEqn(Foam::zero)
:
    quadraticEqn::vsType(Foam::zero{})
{}


inline constexpr Foam::quadraticEqn::quadraticEqn
(
    scalar a,
    scalar b,
    scalar c
)
:
    // aggregate initialization
    quadraticEqn::vsType(Foam::direct_construct_t{}, a, b, c)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::quadraticEqn::value(const scalar x) const
{
    return x*(x*a() + b()) + c();
}


inline Foam::scalar Foam::quadraticEqn::derivative(const scalar x) const
{
    return x*2*a() + b();
}


inline Foam::scalar Foam::quadraticEqn::error(const scalar x) const
{
    return
        SMALL*mag(x)*(mag(x*a()) + mag(b()))
      + SMALL*(mag(x*(x*a() + b())) + mag(c()));
}


// ************************************************************************* //
