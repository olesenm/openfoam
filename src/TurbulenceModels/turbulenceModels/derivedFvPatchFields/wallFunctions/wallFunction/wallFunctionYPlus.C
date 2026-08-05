/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "wallFunctionYPlus.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "volFields.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::wallFunctionYPlus
(
    const turbulenceModel& turbModel,
    const label patchi,
    const bool useWallFunction
)
{
    tmp<volScalarField> tnut = turbModel.nut();
    const volScalarField::Boundary& nutBf = tnut.cref().boundaryField();

    const auto* nutWallPatch =
        isA<nutWallFunctionFvPatchScalarField>(nutBf[patchi]);

    if (useWallFunction && nutWallPatch)
    {
        return nutWallPatch->yPlus();
    }

    const fvPatch& patch = turbModel.mesh().boundary()[patchi];

    if (!isA<wallFvPatch>(patch))
    {
        return tmp<scalarField>::New(patch.size(), Zero);
    }

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuEff = turbModel.nuEff(patchi);
    const scalarField& nuEff = tnuEff.cref();

    const tmp<scalarField> tnu = turbModel.nu(patchi);
    const scalarField& nu = tnu.cref();

    const fvPatchVectorField& U = turbModel.U().boundaryField()[patchi];

    return y*sqrt(nuEff*mag(U.snGrad()))/nu;
}


// ************************************************************************* //
