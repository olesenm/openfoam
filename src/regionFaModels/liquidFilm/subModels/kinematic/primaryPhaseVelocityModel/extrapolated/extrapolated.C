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

#include "extrapolated.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{
    defineTypeNameAndDebug(extrapolated, 0);
    addToRunTimeSelectionTable
    (
        primaryPhaseVelocityModel,
        extrapolated,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::areaSurfaceFilmModels::extrapolated::extrapolated
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    primaryPhaseVelocityModel(typeName, film, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::areaVectorField>
Foam::regionModels::areaSurfaceFilmModels::extrapolated::Up() const
{
    const auto& Uprimary =
        film_.primaryMesh().lookupObject<volVectorField>(film_.UName());

    auto tUp = newUpField();
    auto& Up = tUp.ref();

    // ‌Set ‌up ‌mapping ‌values ‌per ‌patch
    const labelList& patches = mesh_.whichPolyPatches();

    PtrMap<vectorField> patchValues(2*patches.size());

    // U ‌tangential‌on primary
    for (const label patchi : patches)
    {
        const fvPatchVectorField& Uw = Uprimary.boundaryField()[patchi];

        patchValues.set(patchi, -Uw.snGrad());
    }

    // ‌Map ‌U ‌tang ‌to ‌surface
    film_.vsm().mapToSurface(patchValues, Up.primitiveFieldRef());

    Up.primitiveFieldRef() *= film_.h().primitiveField();

    // ‌U ‌tangent ‌on ‌surface
    projectTangent(Up.primitiveFieldRef());

    return tUp;
}


// ************************************************************************* //
