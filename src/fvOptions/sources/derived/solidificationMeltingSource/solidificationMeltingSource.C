/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "solidificationMeltingSource.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "gravityMeshObject.H"
#include "zeroGradientFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(solidificationMeltingSource, 0);
    addToRunTimeSelectionTable(option, solidificationMeltingSource, dictionary);
}
}

const Foam::Enum
<
    Foam::fv::solidificationMeltingSource::thermoMode
>
Foam::fv::solidificationMeltingSource::thermoModeTypeNames_
({
    { thermoMode::mdThermo, "thermo" },
    { thermoMode::mdLookup, "lookup" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::solidificationMeltingSource::readPhaseChangeCoeffs
(
    const dictionary& dict
)
{
    alphaSolidTptr_ = Function1<scalar>::NewIfPresent
    (
        "alphaSolidT",
        dict,
        &mesh_
    );

    const bool hasAlphaSolidT = alphaSolidTptr_.good();
    const bool hasTsol = dict.found("Tsol");
    const bool hasTliq = dict.found("Tliq");
    const bool hasTmelt = dict.found("Tmelt");

    if (hasTliq && !hasTsol)
    {
        FatalIOErrorInFunction(dict)
            << "Entry Tliq requires Tsol"
            << exit(FatalIOError);
    }

    if (hasTsol)
    {
        Tsol_ = dict.get<scalar>("Tsol");
        Tliq_ = dict.getOrDefault<scalar>("Tliq", Tsol_);

        if (Tliq_ < Tsol_)
        {
            FatalIOErrorInFunction(dict)
                << "Liquidus temperature cannot be below solidus temperature"
                << nl
                << "Tsol = " << Tsol_ << nl
                << "Tliq = " << Tliq_
                << exit(FatalIOError);
        }
    }

    if (hasAlphaSolidT)
    {
        phaseChangeMode_ = mdAlphaSolidT;

        if (hasTmelt)
        {
            Tmelt_ = dict.get<scalar>("Tmelt");
            hasTemperatureRef_ = true;
        }
        else if (hasTsol)
        {
            Tmelt_ = Tsol_;
            hasTemperatureRef_ = true;
        }
        else
        {
            hasTemperatureRef_ = false;
        }
    }
    else if (hasTsol)
    {
        phaseChangeMode_ = mdTsolTliq;
        Tmelt_ = Tsol_;
        hasTemperatureRef_ = true;
    }
    else
    {
        phaseChangeMode_ = mdTmelt;
        Tmelt_ = dict.get<scalar>("Tmelt");
        Tsol_ = Tmelt_;
        Tliq_ = Tmelt_;
        hasTemperatureRef_ = true;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::fv::solidificationMeltingSource::Cp() const
{
    switch (mode_)
    {
        case mdThermo:
        {
            const auto& thermo =
                mesh_.lookupObject<basicThermo>(basicThermo::dictName);

            return thermo.Cp();
            break;
        }
        case mdLookup:
        {
            if (CpName_ == "CpRef")
            {
                const scalar CpRef = coeffs_.get<scalar>("CpRef");

                return volScalarField::New
                (
                    IOobject::scopedName(name_, "Cp"),
                    mesh_,
                    dimensionedScalar
                    (
                        "Cp",
                        dimEnergy/dimMass/dimTemperature,
                        CpRef
                    ),
                    fvPatchFieldBase::extrapolatedCalculatedType()
                );
            }
            else
            {
                return mesh_.lookupObject<volScalarField>(CpName_);
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled thermo mode: " << thermoModeTypeNames_[mode_]
                << abort(FatalError);
        }
    }

    return nullptr;
}


void Foam::fv::solidificationMeltingSource::update(const volScalarField& Cp)
{
    if (curTimeIndex_ == mesh_.time().timeIndex())
    {
        return;
    }

    if (debug)
    {
        Info<< type() << ": " << name_ << " - updating phase indicator" << endl;
    }

    if (mesh_.topoChanging())
    {
        deltaT_.resize(cells_.size());
    }

    // update old time alpha1 field
    alpha1_.oldTime();

    const auto& T = mesh_.lookupObject<volScalarField>(TName_);

    forAll(cells_, i)
    {
        label celli = cells_[i];

        scalar Tc = T[celli];
        scalar Cpc = Cp[celli];

        switch (phaseChangeMode_)
        {
            case mdAlphaSolidT:
            {
                const scalar alpha1Target =
                    scalar(1)
                  - Foam::clamp(alphaSolidTptr_->value(Tc), zero_one{});

                const scalar alpha1New =
                    relax_*alpha1Target
                  + (scalar(1) - relax_)*alpha1_[celli];

                alpha1_[celli] = Foam::clamp(alpha1New, zero_one{});
                deltaT_[i] = hasTemperatureRef_ ? (Tc - Tmelt_) : scalar(0);
                break;
            }
            case mdTsolTliq:
            {
                if (Tliq_ == Tsol_)
                {
                    const scalar alpha1New =
                        alpha1_[celli] + relax_*Cpc*(Tc - Tsol_)/L_;

                    alpha1_[celli] = Foam::clamp(alpha1New, zero_one{});
                    deltaT_[i] = Tc - Tsol_;
                }
                else
                {
                    const scalar alpha1Target =
                        Foam::clamp((Tc - Tsol_)/(Tliq_ - Tsol_), zero_one{});

                    const scalar alpha1New =
                        relax_*alpha1Target
                      + (scalar(1) - relax_)*alpha1_[celli];

                    alpha1_[celli] = Foam::clamp(alpha1New, zero_one{});
                    deltaT_[i] =
                        Tc - (Tsol_ + (Tliq_ - Tsol_)*alpha1_[celli]);
                }
                break;
            }
            case mdTmelt:
            {
                const scalar alpha1New =
                    alpha1_[celli] + relax_*Cpc*(Tc - Tmelt_)/L_;

                alpha1_[celli] = Foam::clamp(alpha1New, zero_one{});
                deltaT_[i] = Tc - Tmelt_;
                break;
            }
            default:  // should not happen
                break;
        }
    }

    alpha1_.correctBoundaryConditions();

    curTimeIndex_ = mesh_.time().timeIndex();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidificationMeltingSource::solidificationMeltingSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    mode_(thermoModeTypeNames_.get("thermoMode", coeffs_)),
    phaseChangeMode_(mdTmelt),
    alphaSolidTptr_(nullptr),
    alpha1_
    (
        IOobject
        (
            IOobject::scopedName(name_, "alpha1"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    deltaT_(cells_.size(), 0),
    TName_(coeffs_.getOrDefault<word>("T", "T")),
    CpName_(coeffs_.getOrDefault<word>("Cp", "Cp")),
    UName_(coeffs_.getOrDefault<word>("U", "U")),
    phiName_(coeffs_.getOrDefault<word>("phi", "phi")),
    Tmelt_(0),
    Tsol_(0),
    Tliq_(0),
    L_(coeffs_.get<scalar>("L")),
    relax_(coeffs_.getOrDefault<scalar>("relax", 0.9)),
    rhoRef_(coeffs_.get<scalar>("rhoRef")),
    Cu_(coeffs_.getOrDefault<scalar>("Cu", 100000)),
    q_(coeffs_.getOrDefault<scalar>("q", 0.001)),
    beta_(coeffs_.get<scalar>("beta")),
    curTimeIndex_(-1),
    hasTemperatureRef_(false)
{
    readPhaseChangeCoeffs(coeffs_);

    fieldNames_.resize(2);
    fieldNames_[0] = UName_;

    switch (mode_)
    {
        case mdThermo:
        {
            const auto& thermo =
                mesh_.lookupObject<basicThermo>(basicThermo::dictName);

            fieldNames_[1] = thermo.he().name();
            break;
        }
        case mdLookup:
        {
            fieldNames_[1] = TName_;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled thermo mode: " << thermoModeTypeNames_[mode_]
                << abort(FatalError);
        }
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::solidificationMeltingSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    apply(geometricOneField(), eqn);
}


void Foam::fv::solidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    apply(rho, eqn);
}


void Foam::fv::solidificationMeltingSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField Cp(this->Cp());

    update(Cp);

    const vector& g = meshObjects::gravity::New(mesh_.time()).value();

    scalarField& Sp = eqn.diag();
    vectorField& Su = eqn.source();
    const scalarField& V = mesh_.V();

    forAll(cells_, i)
    {
        label celli = cells_[i];

        scalar Vc = V[celli];
        scalar alpha1c = alpha1_[celli];

        scalar S = -Cu_*sqr(1.0 - alpha1c)/(pow3(alpha1c) + q_);
        vector Sb(rhoRef_*g*beta_*deltaT_[i]);

        Sp[celli] += Vc*S;
        Su[celli] += Vc*Sb;
    }
}


void Foam::fv::solidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    // Momentum source uses a Boussinesq approximation - redirect
    addSup(eqn, fieldi);
}


bool Foam::fv::solidificationMeltingSource::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        readPhaseChangeCoeffs(coeffs_);

        coeffs_.readEntry("L", L_);

        coeffs_.readIfPresent("relax", relax_);

        thermoModeTypeNames_.readEntry("thermoMode", coeffs_, mode_);

        coeffs_.readEntry("rhoRef", rhoRef_);
        coeffs_.readIfPresent("T", TName_);
        coeffs_.readIfPresent("U", UName_);

        coeffs_.readIfPresent("Cu", Cu_);
        coeffs_.readIfPresent("q", q_);

        coeffs_.readEntry("beta", beta_);

        return true;
    }

    return false;
}


// ************************************************************************* //
