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

#include "yPlus.H"
#include "turbulenceModel.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFunctionYPlus.H"
#include "wallFunctionRoughness.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{
    defineTypeNameAndDebug(yPlus, 0);
    addToRunTimeSelectionTable
    (
        primaryPhaseVelocityModel,
        yPlus,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::areaSurfaceFilmModels::yPlus::yPlus
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    primaryPhaseVelocityModel(typeName, film, dict),
    wallCoeffs_(coeffsDict(typeName, dict)),
    yPlusTarget_(coeffsDict(typeName, dict).get<scalar>("yPlus")),
    Ks_(coeffsDict(typeName, dict).getOrDefault<scalar>("Ks", 0)),
    Cs_(coeffsDict(typeName, dict).getOrDefault<scalar>("Cs", 0))
{
    if (yPlusTarget_ < VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Target yPlus = " << yPlusTarget_ << " must be positive."
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::areaVectorField>
Foam::regionModels::areaSurfaceFilmModels::yPlus::Up() const
{
    const fvMesh& mesh = film_.primaryMesh();

    const auto* turbModelPtr = mesh.cfindObject<turbulenceModel>
    (
        turbulenceModel::propertiesName
    );

    if (!turbModelPtr)
    {
        FatalErrorInFunction
            << "yPlus primary-phase velocity model requires a turbulence model"
            << " on the primary region." << nl
            << exit(FatalError);
    }

    const auto& turbModel = *turbModelPtr;

    const auto& Uprimary = mesh.lookupObject<volVectorField>(film_.UName());

    auto tUp = newUpField();
    auto& Up = tUp.ref();

    const List<labelPair>& patchFaces = mesh_.whichPatchFaces();
    const vectorField& nHat = mesh_.faceAreaNormals().internalField();

    const label szPatches = mesh.boundary().size();
    bitSet calculated(szPatches, false);
    List<scalarField> yPlusByPatch(szPatches);
    List<scalarField> yByPatch(szPatches);
    List<scalarField> nuByPatch(szPatches);

    forAll(patchFaces, i)
    {
        const label patchi = patchFaces[i].first();
        const label facei = patchFaces[i].second();

        if (patchi < 0) continue;

        const fvPatch& pp = mesh.boundary()[patchi];
        const label celli = pp.faceCells()[facei];

        const vector Uc(Uprimary[celli]);
        const vector& n = nHat[i];

        const vector Upar(Uc - n*(Uc & n));
        const scalar magUpar = mag(Upar);

        if (magUpar < ROOTVSMALL)
        {
            Up[i] = Zero;
            continue;
        }

        const vector UparHat(Upar/magUpar);

        if (!calculated[patchi])
        {
            constexpr bool useWallFunction = true;

            tmp<scalarField> tyPlus =
                wallFunctionYPlus(turbModel, patchi, useWallFunction);

            yPlusByPatch[patchi] = tyPlus.cref();
            yByPatch[patchi] = turbModel.y()[patchi];

            tmp<scalarField> tnu = turbModel.nu(patchi);
            nuByPatch[patchi] = tnu.cref();

            calculated[patchi] = true;
        }

        const scalar yPlus = yPlusByPatch[patchi][facei];
        const scalar y = yByPatch[patchi][facei];
        const scalar nuw = nuByPatch[patchi][facei];

        if (yPlus < ROOTVSMALL)
        {
            Up[i] = Upar;
            continue;
        }

        if (y < ROOTVSMALL || nuw < ROOTVSMALL)
        {
            FatalErrorInFunction
                << "Invalid wall distance or viscosity on primary patch "
                << mesh.boundary()[patchi].name()
                << " at face " << facei << nl
                << "    y = " << y << nl
                << "    nu = " << nuw << nl
                << exit(FatalError);
        }

        const scalar uStar = yPlus*nuw/y;

        if (uStar < ROOTVSMALL)
        {
            Up[i] = Upar;
            continue;
        }

        scalar Edash = wallCoeffs_.E();
        const scalar KsPlus = uStar*Ks_/nuw;

        if (KsPlus > 2.25)
        {
            Edash /= wallFunctionRoughness(KsPlus, Cs_);
        }

        const scalar magUp =
            (uStar/wallCoeffs_.kappa())
           *Foam::log(Foam::max(Edash*yPlusTarget_, ROOTVSMALL));

        Up[i] = magUp*UparHat;
    }

    projectTangent(Up.primitiveFieldRef());

    return tUp;
}


// ************************************************************************* //
