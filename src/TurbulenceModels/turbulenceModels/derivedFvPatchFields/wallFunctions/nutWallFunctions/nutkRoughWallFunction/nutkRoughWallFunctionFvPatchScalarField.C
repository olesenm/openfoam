/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "nutkRoughWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFunctionRoughness.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::nutkRoughWallFunctionFvPatchScalarField::fnRough
(
    const scalar KsPlus,
    const scalar Cs
) const
{
    return wallFunctionRoughness(KsPlus, Cs);
}


Foam::tmp<Foam::scalarField> Foam::nutkRoughWallFunctionFvPatchScalarField::
calcNut() const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();
    const scalar E = wallCoeffs_.E();

    auto tnutw = tmp<scalarField>::New(*this);
    auto& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        const label celli = patch().faceCells()[facei];

        const scalar uStar = Cmu25*sqrt(k[celli]);
        const scalar yPlus = uStar*y[facei]/nuw[facei];
        const scalar KsPlus = uStar*Ks_[facei]/nuw[facei];

        scalar Edash = E;
        if (2.25 < KsPlus)
        {
            Edash /= fnRough(KsPlus, Cs_[facei]);
        }

        const scalar limitingNutw = max(nutw[facei], nuw[facei]);

        // To avoid oscillations limit the change in the wall viscosity
        // which is particularly important if it temporarily becomes zero
        nutw[facei] =
            max
            (
                min
                (
                    nuw[facei]
                   *(yPlus*kappa/log(max(Edash*yPlus, 1+1e-4)) - 1),
                    2*limitingNutw
                ), 0.5*limitingNutw
            );

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[facei]
                << endl;
        }
    }

    return tnutw;
}


void Foam::nutkRoughWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    Ks_(p.size(), Zero),
    Cs_(p.size(), Zero)
{}


Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper)
{}


Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const this_bctype& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nutkRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    this->parent_bctype::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void Foam::nutkRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    this->parent_bctype::rmap(ptf, addr);

    const auto& nrwfpsf = refCast<const this_bctype>(ptf);

    Ks_.rmap(nrwfpsf.Ks_, addr);
    Cs_.rmap(nrwfpsf.Cs_, addr);
}


void Foam::nutkRoughWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    nutWallFunctionFvPatchScalarField::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutkRoughWallFunctionFvPatchScalarField
    );
}

// ************************************************************************* //
