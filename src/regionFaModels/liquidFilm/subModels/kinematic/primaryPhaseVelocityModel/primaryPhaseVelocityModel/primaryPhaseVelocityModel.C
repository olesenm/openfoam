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

#include "primaryPhaseVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(primaryPhaseVelocityModel, 0);
defineRunTimeSelectionTable(primaryPhaseVelocityModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

primaryPhaseVelocityModel::primaryPhaseVelocityModel
(
    const word& modelType,
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    film_(film),
    mesh_(film.regionMesh())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& primaryPhaseVelocityModel::coeffsDict
(
    const word& modelType,
    const dictionary& dict
)
{
    return dict.subDict(modelType + "Coeffs", keyType::LITERAL);
}


Foam::tmp<Foam::areaVectorField>
primaryPhaseVelocityModel::newUpField() const
{
    return areaVectorField::New
    (
        "Up",
        IOobjectOption::NO_REGISTER,
        mesh_,
        dimensionedVector(dimVelocity, Zero)
    );
}


void primaryPhaseVelocityModel::projectTangent(vectorField& Up) const
{
    const vectorField& nHat = mesh_.faceAreaNormals().internalField();

    Up -= nHat*(Up & nHat);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
