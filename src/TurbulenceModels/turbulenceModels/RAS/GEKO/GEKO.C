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

#include "GEKO.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::validate_input() const
{
    bool hadErrors = false;

    const auto reportError = [&hadErrors, this]() -> Ostream&
    {
        if (!hadErrors)
        {
            FatalIOErrorInFunction(this->coeffDict_);
            hadErrors = true;
        }

        return FatalIOError.stream();
    };

    const auto validatePositive =
        [&reportError](const char* name, scalar value)
    {
        if (value < SMALL)
        {
            reportError()
                << "Invalid range for " << name << ": " << value
                << " - valid range is (0,inf].\n";
        }
    };

    validatePositive("Cmu", Cmu_.value());
    validatePositive("Comega2", Comega2_.value());
    validatePositive("sigmakTilde", sigmakTilde_.value());
    validatePositive("sigmaOmegaTilde", sigmaOmegaTilde_.value());
    validatePositive("CRealize", CRealize_.value());
    validatePositive("kappa", kappa_.value());
    validatePositive("Aplus", Aplus_.value());

    forAll(CSEP_, i)
    {
        if (CSEP_[i] < SMALL)
        {
            reportError()
                << "Invalid value for CSEP field: " << CSEP_[i]
                << ", index: " << i
                << ". Valid range is (SMALL,inf].\n";
        }
    }

    if (hadErrors)
    {
        FatalIOError << exit(FatalIOError);
    }
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correctNut()
{
    const volScalarField S
    (
        sqrt(2.0*magSqr(symm(fvc::grad(this->U_))))
    );

    correctNut(S);
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correctNut(const volScalarField& S)
{
    // Realizability limiter (MM:Eq. A16)
    this->nut_ = k_/max(omega_, S/CRealize_);

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
GEKO<BasicTurbulenceModel>::GEKO
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),


    wallDistanceFree_
    (
        Switch::getOrAddToDict
        (
            "wallDistanceFree",
            this->coeffDict_,
            false
        )
    ),
    productionLimiter_
    (
        Switch::getOrAddToDict
        (
            "productionLimiter",
            this->coeffDict_,
            true
        )
    ),
    katoLaunder_
    (
        Switch::getOrAddToDict
        (
            "katoLaunder",
            this->coeffDict_,
            false
        )
    ),
    dilatationCorrection_
    (
        Switch::getOrAddToDict
        (
            "dilatationCorrection",
            this->coeffDict_,
            false
        )
    ),
    writeCalibrationFields_
    (
        Switch::getOrAddToDict
        (
            "writeCalibrationFields",
            this->coeffDict_,
            false
        )
    ),
    machineLearning_
    (
        Switch::getOrAddToDict
        (
            "machineLearning",
            this->coeffDict_,
            false
        )
    ),

    Cmu_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Comega2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Comega2",
            this->coeffDict_,
            0.083
        )
    ),
    sigmakTilde_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmakTilde",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaOmegaTilde_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaOmegaTilde",
            this->coeffDict_,
            1.17
        )
    ),
    CRealize_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CRealize",
            this->coeffDict_,
            scalar(1.0)/std::sqrt(scalar(3.0))  // (MM:Eq. 17)
        )
    ),
    kappa_
    (
        dimensionedScalar::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CFbTurb_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CFbTurb",
            this->coeffDict_,
            2.0
        )
    ),
    CFbLam_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CFbLam",
            this->coeffDict_,
            1.0  // increase to 25 in case of transition simulations (MM:p. 4)
        )
    ),
    Aplus_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Aplus",
            this->coeffDict_,
            15.0
        )
    ),
    Cc1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cc1",
            this->coeffDict_,
            1.7
        )
    ),
    Cc2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cc2",
            this->coeffDict_,
            1.4
        )
    ),
    CPKlim_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CPKlim",
            this->coeffDict_,
            10.0
        )
    ),


    Csep_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Csep",
            this->coeffDict_,
            1.75
        )
    ),
    Cnw_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cnw",
            this->coeffDict_,
            0.5
        )
    ),
    Cmix_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmix",
            this->coeffDict_,
            0.0
        )
    ),
    Cjet_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cjet",
            this->coeffDict_,
            1.0
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobjectOption::MUST_READ,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobjectOption::MUST_READ,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        ),
        this->mesh_
    ),


    CSEP_
    (
        IOobject
        (
            IOobject::groupName("CSEP", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        ),
        this->mesh_,
        Csep_
    ),
    CNW_
    (
        IOobject
        (
            IOobject::groupName("CNW", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        ),
        this->mesh_,
        Cnw_
    ),
    CMIX_
    (
        IOobject
        (
            IOobject::groupName("CMIX", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        ),
        this->mesh_,
        Cmix_
    ),
    CJET_
    (
        IOobject
        (
            IOobject::groupName("CJET", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        ),
        this->mesh_,
        Cjet_
    ),


    FSEP_
    (
        this->mesh_.newIOobject
        (
            IOobject::groupName("FSEP", alphaRhoPhi.group())
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    FNW_
    (
        this->mesh_.newIOobject
        (
            IOobject::groupName("FNW", alphaRhoPhi.group())
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    FMIX_
    (
        this->mesh_.newIOobject
        (
            IOobject::groupName("FMIX", alphaRhoPhi.group())
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    FJET_
    (
        this->mesh_.newIOobject
        (
            IOobject::groupName("FJET", alphaRhoPhi.group())
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    FBlend_
    (
        this->mesh_.newIOobject
        (
            IOobject::groupName("FBlend", alphaRhoPhi.group())
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)  // wall-distance free variant
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    validate_input();

    if (machineLearning_)
    {
        CkPtr_ = autoPtr<volScalarField>::New
        (
            IOobject
            (
                IOobject::groupName("Ck", alphaRhoPhi.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobjectOption::READ_IF_PRESENT,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            this->mesh_,
            Foam::zero{},
            dimless
        );

        ComegaPtr_ = autoPtr<volScalarField>::New
        (
            IOobject
            (
                IOobject::groupName("Comega", alphaRhoPhi.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobjectOption::READ_IF_PRESENT,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            this->mesh_,
            Foam::zero{},
            dimless
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> GEKO<BasicTurbulenceModel>::R() const
{
    // The gradient is Hessian form (∂U_j/∂x_i), but symm(...) is invariant
    tmp<volTensorField> tgradU = fvc::grad(this->U_);

    tmp<volScalarField> tk = this->k();

    // (MM:Eq. A7)
    // R = (2/3)k I - nu_t * (2S_ij - (2/3)divU δ_ij)
    // The dilatation correction (-(2/3)divU term) is optional
    return volSymmTensorField::New
    (
        IOobject::groupName("R", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        ((2.0/3.0)*I)*tk()
      - this->nut_
       *(
            dilatationCorrection_
          ? devTwoSymm(tgradU())
          : (twoSymm(tgradU()))
        ),
        tk().boundaryField().types()
    );
}


template<class BasicTurbulenceModel>
bool GEKO<BasicTurbulenceModel>::read()
{
    if (!eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        return false;
    }

    wallDistanceFree_.readIfPresent("wallDistanceFree", this->coeffDict());
    productionLimiter_.readIfPresent("productionLimiter", this->coeffDict());
    katoLaunder_.readIfPresent("katoLaunder", this->coeffDict());
    dilatationCorrection_.readIfPresent
    (
        "dilatationCorrection",
        this->coeffDict()
    );
    writeCalibrationFields_.readIfPresent
    (
        "writeCalibrationFields",
        this->coeffDict()
    );
    // Do not allow changing machineLearning switch at runtime since the fields
    // are only conditionally read and solved

    Cmu_.readIfPresent(this->coeffDict());
    Comega2_.readIfPresent(this->coeffDict());
    sigmakTilde_.readIfPresent(this->coeffDict());
    sigmaOmegaTilde_.readIfPresent(this->coeffDict());
    CRealize_.readIfPresent(this->coeffDict());
    kappa_.readIfPresent(this->coeffDict());
    CFbTurb_.readIfPresent(this->coeffDict());
    CFbLam_.readIfPresent(this->coeffDict());
    Aplus_.readIfPresent(this->coeffDict());
    Cc1_.readIfPresent(this->coeffDict());
    Cc2_.readIfPresent(this->coeffDict());
    CPKlim_.readIfPresent(this->coeffDict());

    validate_input();

    Csep_.readIfPresent(this->coeffDict());
    Cnw_.readIfPresent(this->coeffDict());
    Cmix_.readIfPresent(this->coeffDict());
    Cjet_.readIfPresent(this->coeffDict());

    return true;
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Calculate the velocity gradient tensor in Hessian form (delta_i u_j)
    // Tranpose of the classical Jacobian form (delta_j u_i)
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU.cref();

    // Strain-rate tensor magnitude (MM:Eq. 4)
    const volScalarField S(sqrt(2.0*magSqr(symm(gradU))));

    // Rotation tensor magnitude (MM:Eq. A14)
    const volScalarField OmegaMag(sqrt(2.0*magSqr(skew(gradU))));


    // Calculate the production terms

    // Turbulent kinetic energy production rate per mass (MM:Eq. 15)
    tmp<volScalarField> tG;
    if (katoLaunder_)
    {
        tG = tmp<volScalarField>::New(this->GName(), nut*S*OmegaMag);
    }
    else
    {
        tG = tmp<volScalarField>::New(this->GName(), nut*sqr(S));
    }
    volScalarField& G = tG.ref();

    // Omega production rate per-mass (MM:Eq. A4):
    const volScalarField Pomega_m
    (
        G
       /max
       (
            nut + scalar(0.001)*this->nu(),
            dimensionedScalar(nut.dimensions(), SMALL)
       )
    );

    if (productionLimiter_)
    {
        // Limited production per-mass (MM:Eq. A15)
        G.clamp_max(CPKlim_*Cmu_*k_*omega_);
    }


    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();


    // Calculate the blending function field (MM:Eq. A8)

    if (!wallDistanceFree_)
    {
        const volScalarField& y = wallDist::New(this->mesh_).y();

        tmp<volScalarField> tkTilde = max(k_, CFbLam_*this->nu()*omega_);

        tmp<volScalarField> tLT = sqrt(tkTilde)/(Cmu_*omega_);

        tmp<volScalarField> txBlend = CFbTurb_*tLT/y;

        FBlend_ = Foam::tanh(pow4(txBlend));
        FBlend_.correctBoundaryConditions();
    }


    // Calculate the blending-function coefficient fields (MM:Eqs. A17-A18)
    // sigmaK or sigmaOmega cannot be zero
    const volScalarField sigmaK
    (
        max
        (
            sigmakTilde_*(scalar(1) + (CSEP_ - scalar(1))*scalar(0.25)*FBlend_),
            dimensionedScalar(dimless, SMALL)
        )
    );

    const volScalarField sigmaOmega
    (
        max
        (
            sigmaOmegaTilde_*(scalar(1) + (CSEP_ - scalar(1))*FBlend_),
            dimensionedScalar(dimless, SMALL)
        )
    );

    const scalar sqrtCmu = Foam::sqrt(Cmu_.value());
    const volScalarField Comega1
    (
        (scalar(1)/Cmu_)
       *(
            Comega2_
          - sqr(kappa_)*sqrtCmu/sigmaOmega
        )
    );


    // Calculate the flow-separation calibration field, FSEP (MM:Eqs. A9-10)

    const volScalarField yTildePlus
    (
        (scalar(1)/kappa_)*(k_/omega_)/this->nu()
    );

    const volScalarField fD
    (
        scalar(1)/(scalar(1) + sqr(yTildePlus/Aplus_))
    );

    const volScalarField Cpsi
    (
        scalar(0.2454)*Foam::pow(CSEP_, scalar(-0.803))
    );

    const volScalarField psi(Cpsi*(CSEP_ - scalar(1))*FBlend_);

    FSEP_ = scalar(1) + fD*psi;


    // Calculate the jet-flow calibration field, FJET (MM:Eqs. A14)

    const volScalarField Soo
    (
        min
        (
            min(S, OmegaMag),
            omega_
        )
       /(scalar(0.3)*omega_)
    );

    const volScalarField xJet
    (
        scalar(4.0)
       *(
            (Soo - scalar(0.9))/scalar(0.15)
          - scalar(0.5)
        )
    );

    FJET_ = scalar(0.5)*(scalar(1) + Foam::tanh(xJet));


    // Calculate the mixing function field, FMIX (MM:Eq. A13)

    FMIX_ =
        FSEP_
      + (
            CMIX_
          + scalar(0.13)*CJET_*(FJET_ - scalar(1))
        )
       *(scalar(1) - FBlend_);


    // Calculate the non-equilibrium calibration field, FNW (MM:Eqs. A11)

    // Cross-diffusion term (MM:Eq. A6)
    const volScalarField CD
    (
        rho*(scalar(2)/sigmaOmega)
       *(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    // Destruction term (MM:Eq. A12)
    const volScalarField Destr(Comega2_*FMIX_*rho*sqr(omega_));

    const volScalarField FLim(Destr/max(Destr, scalar(2)*CD));

    forAll(FNW_, celli)
    {
        if (CD[celli] > scalar(0))
        {
            FNW_[celli] = scalar(1);
        }
        else
        {
            FNW_[celli] =
                (
                    CNW_[celli]
                  - (Cc1_.value() + Cc2_.value()*CNW_[celli])*fD[celli]
                )*FLim[celli];
        }
    }


    // Calculate the specific dissipation rate field (MM:Eq. A2)

    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(sigmaOmega), omega_)
     ==
        alpha()*rho()*Comega1()*Pomega_m()
      - fvm::Sp(alpha()*rho()*Comega2_*FMIX_()*omega_(), omega_)
      + alpha()*FNW_()*CD()
      + fvOptions(alpha, rho, omega_)
    );

    if (machineLearning_)
    {
        const auto& Comega = *ComegaPtr_;

        volScalarField::Internal mlomega
        (
            alpha()*rho()*Comega()*sqr(omega_())
        );
        omegaEqn.ref() -= mlomega;
    }

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Calculate the turbulent kinetic energy field (MM:Eq. A1)

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(sigmaK), k_)
     ==
        alpha()*rho()*G()
      - fvm::Sp(alpha()*rho()*Cmu_*omega_(), k_)
      + fvOptions(alpha, rho, k_)
    );

    if (machineLearning_)
    {
        const auto& Ck = *CkPtr_;

        volScalarField::Internal mlk(alpha()*rho()*Ck()*k_()*omega_());
        kEqn.ref() -= mlk;
    }

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    kEqn.ref().boundaryManipulate(k_.boundaryFieldRef());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    // Update nut with realizability limiter (MM:Eq. A16)
    correctNut(S);

    // Write the calibration and machine-learning fields if enabled
    if (this->mesh_.time().writeTime())
    {
        if (writeCalibrationFields_)
        {
            // The calibration fields are only accessed via internal-field-only
            // in the fvMatrix assembly. Their boundary values are never
            // consumed; therefore, MPI synchronisation is avoided except field
            // writes.
            FSEP_.correctBoundaryConditions();
            FNW_.correctBoundaryConditions();
            FMIX_.correctBoundaryConditions();
            FJET_.correctBoundaryConditions();

            FSEP_.write();
            FNW_.write();
            FMIX_.write();
            FJET_.write();
            FBlend_.write();
        }

        if (CkPtr_)
        {
            CkPtr_->write();
        }
        if (ComegaPtr_)
        {
            ComegaPtr_->write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
