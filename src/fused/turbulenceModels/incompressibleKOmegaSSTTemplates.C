/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 M.Janssens
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

#include "fvCFD.H"
#include "GeometricFieldExpression.H"
#include "fvOptions.H"
#include "wallDist.H"

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


//- Wrap of constant as an expression
typedef Expression::UniformGeometricFieldWrap<volTensorField> vtConstant;
//- Wrap of volTensorField as an expression
typedef Expression::GeometricFieldConstRefWrap<volTensorField> vtExpr;
//- Evaluator of an expression
typedef Expression::GeometricFieldRefWrap<volTensorField> vtEvaluator;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CDType>
auto incompressibleKOmegaSST::F1(const CDType& CDkOmega) const -> decltype(auto)
{
    // Field for sizing constants
    const volScalarField& dummy = this->omega_;

    // Constants in volScalarField form
    const volConstant one(dummy, scalar(1));
    const volConstant hundred(dummy, scalar(100));
    const volConstant fiveHundred(dummy, scalar(500));
    const volConstant betaStar(dummy, betaStar_);
    const volConstant twoOverBetaStar(dummy, scalar(2)/betaStar_);
    const volConstant rho(dummy, 1.0);

    const volScalarField& y = wallDist::New(this->mesh_).y();

    // Wrap some fields
    const auto wk = volExpr(k_);
    const auto wy = volExpr(y);
    const auto womega = volExpr(omega_);
    const auto wmu = volExpr(this->mu());

    auto CDkOmegaPlus = max
    (
        CDkOmega,
        volConstant(dummy, dimensionedScalar(dimless/sqr(dimTime), 1.0e-10))
    );

    auto arg1 = min
    (
        min
        (
            max
            (
                (one/betaStar)*sqrt(wk)/(womega*wy),
                fiveHundred*(wmu/rho)/(sqr(wy)*womega)
            ),
            volConstant(dummy, 4*alphaOmega2_)*wk/(CDkOmegaPlus*sqr(wy))
        ),
        volConstant(dummy, 10)
    );

    return tanh(pow4(arg1));
}


auto incompressibleKOmegaSST::F2() const -> decltype(auto)
{
    // Field for sizing constants. Only sizes used, not values. Should be of
    // same type as constant and persist until after evaluation.
    const volScalarField& dummy = this->omega_;

    // Constants in volScalarField form
    const volConstant hundred(dummy, scalar(100));
    const volConstant fiveHundred(dummy, scalar(500));
    const volConstant twoOverBetaStar(dummy, scalar(2)/betaStar_);

    const volConstant rho(dummy, 1.0);

    const volScalarField& y = wallDist::New(this->mesh_).y();

    // Wrap some fields
    const auto wk = volExpr(k_);
    const auto wy = volExpr(y);
    const auto womega = volExpr(omega_);
    const auto wmu = volExpr(this->mu()());

    auto arg2 = min
    (
        max
        (
            twoOverBetaStar * sqrt(wk) / (womega*wy),
            fiveHundred
          * (wmu / rho)
          / (sqr(wy)*womega)
        ),
        hundred
    );

    return tanh(sqr(arg2));
}


auto incompressibleKOmegaSST::F3() const -> decltype(auto)
{
    // Field for sizing constants. Only sizes used, not values. Should be of
    // same type as constant and persist until after evaluation.
    const volScalarField& dummy = this->omega_;

    // Constants in volScalarField form
    const volConstant c150(dummy, scalar(150));
    const volConstant c10(dummy, scalar(10));
    const volConstant rho(dummy, 1.0);

    const volScalarField& y = wallDist::New(this->mesh_).y();

    // Wrap some fields
    const auto wy = volExpr(y);
    const auto womega = volExpr(omega_);
    const auto wmu = volExpr(this->mu()());

    auto arg3 = min
    (
        c150*(wmu/rho)/(womega*sqr(wy)),
        c10
    );

    return volConstant(dummy, 1) - tanh(pow4(arg3));
}


template<class S2Type>
void incompressibleKOmegaSST::correctNut(const S2Type& S2)
{
    // Field for sizing constants
    const volScalarField& dummy = this->omega_;

    const volConstant a1(dummy, a1_);
    const volConstant b1(dummy, b1_);

    // Wrap some fields
    const auto wk = volExpr(k_);
    const auto womega = volExpr(omega_);
    auto wnut = volEvaluator(this->nut_);

    // Correct the turbulence viscosity
    if (F3_)
    {
        wnut = a1*wk/max(a1*womega, b1*F2()*F3()*sqrt(S2));
    }
    else
    {
        wnut = a1*wk/max(a1*womega, b1*F2()*sqrt(S2));
    }
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class E1>
auto incompressibleKOmegaSST::Pk(const E1& G) const -> decltype(auto)
{
    // Field for sizing constants. Only sizes used, not values. Should be of
    // same type as constant and persist until after evaluation.
    const label nCells = this->omega_.size();

    // Constants in scalarList form
    const constant c1(nCells, c1_);
    const constant betaStar(nCells, betaStar_);

    // Wrap some fields
    const auto wk = expr(k_.internalField());
    const auto womega = expr(omega_.internalField());

    return min(G, (c1*betaStar)*wk*womega);
}


template<class E1, class E2, class E3>
auto incompressibleKOmegaSST::GbyNu
(
    const E1& GbyNu0,
    const E2& F2,
    const E3& S2
) const -> decltype(auto)
{
    // Field for sizing constants. Only sizes used, not values. Should be of
    // same type as constant and persist until after evaluation.
    const label nCells = this->omega_.size();

    // Constants in scalarList form
    const constant a1(nCells, a1_);
    const constant b1(nCells, b1_);
    const constant c1(nCells, c1_);
    const constant betaStar(nCells, betaStar_);

    // Wrap some fields
    const auto womega = expr(omega_.internalField());

    return min
    (
        GbyNu0,
        (c1/a1)*betaStar*womega*max(a1*womega, b1*F2*sqrt(S2))
    );
}


//- Return the blended field
template<class E1>
auto incompressibleKOmegaSST::blend
(
    const E1& F1,
    const scalar& p1,
    const scalar& p2
) const -> decltype(auto)
{
    // Field for sizing constants
    const volScalarField& dummy = this->omega_;
    const volConstant psi1(dummy, p1);
    const volConstant psi2(dummy, p2);

    return F1*(psi1 - psi2) + psi2;
}


template<class E1>
auto incompressibleKOmegaSST::alphaK(const E1& F1) const -> decltype(auto)
{
    return blend(F1, alphaK1_.value(), alphaK2_.value());
}


template<class E1>
auto incompressibleKOmegaSST::alphaOmega(const E1& F1) const -> decltype(auto)
{
    return blend(F1, alphaOmega1_.value(), alphaOmega2_.value());
}

template<class E1>
auto incompressibleKOmegaSST::beta(const E1& F1) const -> decltype(auto)
{
    return blend(F1, beta1_.value(), beta2_.value());
}


template<class E1>
auto incompressibleKOmegaSST::gamma(const E1& F1) const -> decltype(auto)
{
    return blend(F1, gamma1_.value(), gamma2_.value());
}


template<class E1>
auto incompressibleKOmegaSST::DkEff(const E1& F1) const -> decltype(auto)
{
    return alphaK(F1)*volExpr(this->nut_) + volExpr(this->nu());
}


template<class E1>
auto incompressibleKOmegaSST::DomegaEff(const E1& F1) const -> decltype(auto)
{
    return alphaOmega(F1)*volExpr(this->nut_) + volExpr(this->nu());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
