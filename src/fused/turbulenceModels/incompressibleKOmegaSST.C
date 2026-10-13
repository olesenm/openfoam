/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvMatrixExpression.H"
#include "incompressibleKOmegaSST.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
#include "fvCFD.H"

#include "GeometricFieldExpression.H"
#include "addToRunTimeSelectionTable.H"

#include <ratio>
#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{


//- Wrap of constant as a list expression
typedef Expression::UniformListWrap<scalar> constant;
//- Wrap of List as an expression
typedef Expression::ListConstRefWrap<scalar> expr;
//- Evaluator of an expression
typedef Expression::ListRefWrap<scalar> evaluator;


//- Wrap of constant as an expression
typedef Expression::UniformGeometricFieldWrap<volScalarField> volConstant;
//- Wrap of volScalarField as an expression
typedef Expression::GeometricFieldConstRefWrap<volScalarField> volExpr;
//- Evaluator of an expression
typedef Expression::GeometricFieldRefWrap<volScalarField> volEvaluator;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(incompressibleKOmegaSST, 0);
addToRunTimeSelectionTable(RASModel, incompressibleKOmegaSST, dictionary);

static double time_total = 0;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void incompressibleKOmegaSST::correctNut()
{
    // Correct the turbulence viscosity
    //correctNut(2*magSqr(symm(fvc::grad(this->U_))));

    // Field for sizing constants. Only sizes used, not values. Should be of
    // same type as constant.
    const volScalarField& dummy = this->omega_;

    // Correct the turbulence viscosity
    auto tgradU = fvc::grad(this->U_);
    const volScalarField magSqrSymmGrad(Foam::magSqr(Foam::symm(tgradU)));
    correctNut(volConstant(dummy, 2.0)*volExpr(magSqrSymmGrad));
}


Foam::tmp<Foam::volScalarField> incompressibleKOmegaSST::S2
(
    const volTensorField& gradU
) const
{
    return 2*magSqr(symm(gradU));
}


tmp<volScalarField::Internal> incompressibleKOmegaSST::GbyNu0
(
    const volTensorField& gradU,
    const volScalarField& /* S2 not used */
) const
{
    return tmp<volScalarField::Internal>::New
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        gradU() && devTwoSymm(gradU())
    );
}


tmp<fvScalarMatrix> incompressibleKOmegaSST::kSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        k_,
        dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
    );
}


tmp<fvScalarMatrix> incompressibleKOmegaSST::omegaSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        omega_,
        dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
    );
}


tmp<fvScalarMatrix> incompressibleKOmegaSST::Qsas
(
    const volScalarField::Internal& S2
    //const volScalarField::Internal& gamma,
    //const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>::New
    (
        omega_,
        dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

incompressibleKOmegaSST::incompressibleKOmegaSST
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<incompressible::RASModel>
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

    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void incompressibleKOmegaSST::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}


bool incompressibleKOmegaSST::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}


void incompressibleKOmegaSST::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    const std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();

    // Sizing
    const label nCells = this->nut_.size();

    // Create some expression constants
    const constant twoThird(nCells, (2.0/3.0));
    const constant one(nCells, 1.0);
    const constant two(nCells, 2.0);
    const constant minOne(nCells, -1.0);

    const constant alpha(nCells, 1.0);
    const constant rho(nCells, 1.0);

    const constant alphaK1(nCells, alphaK1_);
    const constant alphaK2(nCells, alphaK2_);

    const constant alphaOmega1(nCells, alphaOmega1_);
    const constant alphaOmega2(nCells, alphaOmega2_);

    const constant gamma1(nCells, gamma1_);
    const constant gamma2(nCells, gamma2_);

    const constant beta1(nCells, beta1_);
    const constant beta2(nCells, beta2_);

    const constant betaStar(nCells, betaStar_);

    const constant a1(nCells, a1_);
    const constant b1(nCells, b1_);
    const constant c1(nCells, c1_);

    const constant kInf(nCells, kInf_);
    const constant omegaInf(nCells, omegaInf_);

    const volConstant volTwo(this->nut_, 2.0);
    const volConstant volAlphaOmega2(this->nut_, alphaOmega2_);


    // Local references
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<incompressible::RASModel>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));
    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    // Construct G with correct dimensions and use expression templates
    // to set values. TBD: constructor for DimensionedField with ET?
    volScalarField::Internal G(this->GName(), GbyNu0);
    static_cast<scalarList&>(G) = expr(nut.internalField())*expr(GbyNu0);

    // - boundary condition changes a cell value
    // - normally this would be triggered through correctBoundaryConditions
    // - which would do
    //      - fvPatchField::evaluate() which calls
    //      - fvPatchField::updateCoeffs()
    // - however any processor boundary conditions already start sending
    //   at initEvaluate so would send over the old value.
    // - avoid this by explicitly calling updateCoeffs early and then
    //   only doing the boundary conditions that rely on initEvaluate
    //   (currently only coupled ones)

    //- 1. Explicitly swap values on coupled boundary conditions
    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    // omegaWallFunctions change the cell value! Make sure to push these to
    // coupled neighbours. Note that we want to avoid the re-updateCoeffs
    // of the wallFunctions so make sure to bypass the evaluate on
    // those patches and only do the coupled ones.
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    ////- 2. Make sure the boundary condition calls updateCoeffs from
    ////     initEvaluate
    ////     (so before any swap is done - requires all coupled bcs to be
    ////      after wall bcs. Unfortunately this conflicts with cyclicACMI)
    //omega_.correctBoundaryConditions();

    auto tkOmega(fvc::grad(k_) & fvc::grad(omega_));

    //const volScalarField CDkOmega
    //(
    //    (2*alphaOmega2_)*(tkOmega)/omega_
    //);
    const auto eCDkOmega
    (
        (volTwo*volAlphaOmega2)*volExpr(tkOmega())/volExpr(omega_)
    );

    // Evaluate into F1
    const volScalarField F1("F1", mesh(), this->F1(eCDkOmega));
    //const volScalarField F23(this->F23());
    const auto eF2(this->F2().internalField());

    // Create temporary storage with calculated bc
    volScalarField work(scalar(1.0)*F1);
    // Fix dimensions. TDB : add to expressions
    work.dimensions() *= sqr(dimLength)/dimTime;

    {
        //const volScalarField::Internal gamma(this->gamma(F1));
        // Get expression for gamma
        const auto egamma(this->gamma(volExpr(F1)));
        const auto ebeta(this->beta(volExpr(F1)));

        //GbyNu0 = GbyNu(GbyNu0, F23(), S2());
        // Note: operate only on List, not DimensionedField level
        scalarList& GbyNu0i = GbyNu0;
        if (F3_)
        {
            // F23 is F2*F3
            const auto eF23 = eF2*F3().internalField();
            GbyNu0i = GbyNu(GbyNu0i.expr(), eF23, expr(S2()));
        }
        else
        {
            // F23 is same as F2
            GbyNu0i = GbyNu(GbyNu0i.expr(), eF2, expr(S2()));
        }

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian
            (
                volEvaluator(work).evaluate(DomegaEff(F1.expr())),
                omega_
            )
         ==
        //     alpha()*rho()*gamma*GbyNu0
        //   - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
        //   - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
        //   - fvm::SuSp
        //     (
        //         alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
        //         omega_
        //     )
        //   + alpha()*rho()*beta*sqr(omegaInf_)
            Qsas(S2())      //, gamma, beta)
          + omegaSource()
          + fvOptions(omega_)
        );
        {
            // Source terms as expression templates
            auto& m = omegaEqn.ref();

            m =
            (
                m.expr()
            ==
                Expression::Su(alpha*rho*egamma.internalField()*GbyNu0, m)
              - Expression::SuSp
                (
                    twoThird*alpha*rho*egamma.internalField()*divU,
                    m
                )
              - Expression::Sp(alpha*rho*ebeta.internalField()*omega_(), m)
              - Expression::SuSp
                (
                    alpha*rho*(expr(F1()) - one)*eCDkOmega.internalField()/omega_(),
                    m
                )
              + Expression::Su(alpha*rho*ebeta.internalField()*sqr(omegaInf), m)
            );
        }

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }


    {
        // Turbulent kinetic energy equation
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(k_)
          + fvm::div(alphaRhoPhi, k_)
          - fvm::laplacian
            (
                volEvaluator(work).evaluate(DkEff(volExpr(F1))),
                k_
            )
         ==
        //     alpha()*rho()*Pk(G)
        //   - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
        //   - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
        //   + alpha()*rho()*betaStar_*omegaInf_*kInf_
            kSource()
          + fvOptions(k_)
        );

        {
            // Source terms as expression templates
            fvScalarMatrix& m = kEqn.ref();

            m =
            (
                m.expr()
             == Expression::Su(alpha*rho*Pk(G.expr()), m)
              - Expression::SuSp(twoThird*alpha*rho*expr(divU), m)
              - Expression::Sp(alpha*rho*betaStar*expr(omega_), m)
              + Expression::Su(alpha*rho*betaStar*omegaInf*kInf, m)
            );
        }


        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(k_);
        bound(k_, this->kMin_);
    }

    correctNut(volExpr(S2));

    const std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time_span = t2 - t1;
    time_total += time_span.count();
    Info<< "incompressibleKOmegaSST: spent:" << time_span.count()
        << " total:" << time_total << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
