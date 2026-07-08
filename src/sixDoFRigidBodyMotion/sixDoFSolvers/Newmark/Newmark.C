/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "Newmark.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{
    defineTypeNameAndDebug(Newmark, 0);
    addToRunTimeSelectionTable(sixDoFSolver, Newmark, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::Newmark::Newmark
(
    const dictionary& dict,
    sixDoFRigidBodyMotion& body
)
:
    sixDoFSolver(dict, body),
    gamma_(dict.getOrDefault<scalar>("gamma", 0.5)),
    beta_
    (
        max
        (
            0.25*sqr(gamma_ + 0.5),
            dict.getOrDefault<scalar>("beta", 0.25)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::Newmark::~Newmark()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFSolvers::Newmark::solve
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

    // Implicit Newmark-beta Version 1 translational acceleration
    newmarkTranslationalAcceleration
    (
        fGlobal,
        gamma_,
        beta_,
        deltaT,
        firstIter
    );

    // Update torque from external moment and restraint moments
    updateTorque(tauGlobal, firstIter);

    // Correct linear velocity
    v() =
        tConstraints()
      & (
            v0()
          + aDamp()*deltaT*((1 - gamma_)*a0() + gamma_*a())
        );

    // Correct position
    centreOfRotation() =
        centreOfRotation0()
      + (
            tConstraints()
          & (
                deltaT*v0()
              + aDamp()*0.5*sqr(deltaT)*((1 - 2*beta_)*a0() + 2*beta_*a())
            )
        );

    // Correct angular momentum
    pi() =
        rConstraints()
      & (
            pi0()
          + aDamp()*deltaT*((1 - gamma_)*tau0() + gamma_*tau())
        );

    // Correct orientation
    vector piDeltaT =
        rConstraints()
      & (
            deltaT*pi0()
          + aDamp()*0.5*sqr(deltaT)*((1 - 2*beta_)*tau0() + 2*beta_*tau())
        );
    Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
    Q() = Qpi.first();
}

// ************************************************************************* //
