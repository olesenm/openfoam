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

Class
    Foam::wallFunctionYPlus

Description
    Calculate the yPlus value for patches of turbulence-enabled simulations.

    For wall patches using a nut wall function, the yPlus calculation of
    the wall-function is used. Otherwise, the classic yPlus expression is
    calculated.

SourceFiles
    wallFunctionYPlus.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallFunctionYPlus_H
#define Foam_wallFunctionYPlus_H

#include "labelFwd.H"
#include "primitiveFieldsFwd.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class turbulenceModel;

//- Return the yPlus value for a turbulence-model-enabled patch
tmp<scalarField> wallFunctionYPlus
(
    const turbulenceModel&,
    const label patchi,
    const bool useWallFunction
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
