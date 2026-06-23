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

#include "prescribed.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{
    defineTypeNameAndDebug(prescribed, 0);
    addToRunTimeSelectionTable
    (
        primaryPhaseVelocityModel,
        prescribed,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::areaSurfaceFilmModels::prescribed::prescribed
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    primaryPhaseVelocityModel(typeName, film, dict),
    value_(Zero),
    fieldName_(word::null)
{
    const dictionary& coeffs = coeffsDict(typeName, dict);

    if (!coeffs.readIfPresent("field", fieldName_))
    {
        coeffs.readEntry("value", value_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::areaVectorField>
Foam::regionModels::areaSurfaceFilmModels::prescribed::Up() const
{
    auto tUp = newUpField();
    auto& Up = tUp.ref();

    if (!fieldName_.empty())
    {
        const auto& Ufield =
            film_.primaryMesh().lookupObject<volVectorField>(fieldName_);

        film_.vsm().mapInternalToSurface(Ufield, Up.primitiveFieldRef());
    }
    else
    {
        Up.primitiveFieldRef() = value_;
    }

    projectTangent(Up.primitiveFieldRef());

    return tUp;
}


// ************************************************************************* //
