/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "CrankNicolson.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{
    defineTypeNameAndDebug(CrankNicolson, 0);
    addToRunTimeSelectionTable(sixDoFSolver, CrankNicolson, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::CrankNicolson::CrankNicolson
(
    const dictionary& dict,
    sixDoFRigidBodyMotion& body
)
:
    sixDoFSolver(dict, body),
    aoc_(dict.getOrDefault<scalar>("aoc", 0.5)),
    voc_(dict.getOrDefault<scalar>("voc", 0.5))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::CrankNicolson::~CrankNicolson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFSolvers::CrankNicolson::solve
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    // Update the constraints to the object
    updateConstraints();

    // Implicit acceleration solve. For the default aoc=voc=0.5 this is
    // equivalent to Newmark with gamma=0.5 and beta=0.25.
    const scalar beta = voc_*aoc_;

    newmarkTranslationalAcceleration
    (
        fGlobal,
        aoc_,
        beta,
        deltaT,
        firstIter
    );

    // Update torque from external moment and restraint moments
    updateTorque(tauGlobal, firstIter);

    // Correct linear velocity
    v() = tConstraints()
      & (v0() + aDamp()*deltaT*(aoc_*a() + (1 - aoc_)*a0()));

    // Correct position
    centreOfRotation() =
        centreOfRotation0() + deltaT*(voc_*v() + (1 - voc_)*v0());

    // Correct angular momentum
    pi() = rConstraints()
      & (pi0() + aDamp()*deltaT*(aoc_*tau() + (1 - aoc_)*tau0()));

    // Correct orientation
    Tuple2<tensor, vector> Qpi =
        rotate(Q0(), (voc_*pi() + (1 - voc_)*pi0()), deltaT);
    Q() = Qpi.first();
}

// ************************************************************************* //
