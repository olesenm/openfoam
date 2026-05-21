/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::RASModels::incompressibleKOmegaSST

Group
    grpRASTurbulence

Description
    Implementation of the k-omega-SST turbulence model for
    incompressible and compressible flows.

    Light wrapper around base class.

    Turbulence model described in:
    \verbatim
        Menter, F. R. & Esch, T. (2001).
        Elements of Industrial Heat Transfer Prediction.
        16th Brazilian Congress of Mechanical Engineering (COBEM).
    \endverbatim

    with updated coefficients from
    \verbatim
        Menter, F. R., Kuntz, M., and Langtry, R. (2003).
        Ten Years of Industrial Experience with the SST Turbulence Model.
        Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano,
        & M. Tummers, Begell House, Inc., 625 - 632.
    \endverbatim

    but with the consistent production terms from the 2001 paper as form in the
    2003 paper is a typo, see
    \verbatim
        http://turbmodels.larc.nasa.gov/sst.html
    \endverbatim

    and the addition of the optional F3 term for rough walls from
    \verbatim
        Hellsten, A. (1998).
        "Some Improvements in Menter’s k-omega-SST turbulence model"
        29th AIAA Fluid Dynamics Conference, AIAA-98-2554.
    \endverbatim

    Note that this implementation is written in terms of alpha diffusion
    coefficients rather than the more traditional sigma (alpha = 1/sigma) so
    that the blending can be applied to all coefficuients in a consistent
    manner.  The paper suggests that sigma is blended but this would not be
    consistent with the blending of the k-epsilon and k-omega models.

    Also note that the error in the last term of equation (2) relating to
    sigma has been corrected.

    Wall-functions are applied in this implementation by using equations (14)
    to specify the near-wall omega as appropriate.

    The blending functions (15) and (16) are not currently used because of the
    uncertainty in their origin, range of applicability and that if y+ becomes
    sufficiently small blending u_tau in this manner clearly becomes nonsense.

    The default model coefficients are
    \verbatim
        incompressibleKOmegaSSTCoeffs
        {
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;
        }
    \endverbatim

SourceFiles
    incompressibleKOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleKOmegaSST_H
#define incompressibleKOmegaSST_H

//#include "incompressibleKOmegaSSTBase.H"
#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                   Class incompressibleKOmegaSST Declaration
\*---------------------------------------------------------------------------*/

class incompressibleKOmegaSST
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private Member Functions

        //- No copy construct
        incompressibleKOmegaSST(const incompressibleKOmegaSST&) = delete;

        //- No copy assignment
        void operator=(const incompressibleKOmegaSST&) = delete;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar b1_;
            dimensionedScalar c1_;

            //- Flag to include the F3 term
            Switch F3_;


        // Fields

            //- Turbulent kinetic energy field [m^2/s^2]
            volScalarField k_;

            //- Specific dissipation rate field [1/s]
            volScalarField omega_;


        // Decay control

            //- Flag to include the decay control
            Switch decayControl_;
            dimensionedScalar kInf_;
            dimensionedScalar omegaInf_;


    // Protected Member Functions

        //- Set decay control with kInf and omegaInf
        void setDecayControl(const dictionary& dict);


        // Expression templates

            template<class CDType>
            auto F1(const CDType&) const -> decltype(auto);

            auto F2() const -> decltype(auto);

            auto F3() const -> decltype(auto);

            template<class S2Type>
            void correctNut(const S2Type& S2);

            //- Return k production rate
            template<class E1>
            auto Pk(const E1& G) const -> decltype(auto);

            //- Return G/nu
            template<class E1, class E2, class E3>
            auto GbyNu(const E1&, const E2&, const E3&) const
                -> decltype(auto);

            //- Return the blended field
            template<class E1>
            auto blend
            (
                const E1& F1,
                const scalar& p1,
                const scalar& p2
            ) const -> decltype(auto);

            template<class E1>
            auto alphaK(const E1& F1) const -> decltype(auto);

            template<class E1>
            auto alphaOmega(const E1& F1) const -> decltype(auto);

            template<class E1>
            auto beta(const E1& F1) const -> decltype(auto);

            template<class E1>
            auto gamma(const E1& F1) const -> decltype(auto);


            //- Return the effective diffusivity for k
            template<class E1>
            auto DkEff(const E1& F1) const -> decltype(auto);

            //- Return the effective diffusivity for omega
            template<class E1>
            auto DomegaEff(const E1& F1) const -> decltype(auto);


        virtual void correctNut();

        //- Return square of strain rate
        virtual tmp<volScalarField> S2
        (
            const volTensorField& gradU
        ) const;

        ////- Return k production rate
        //virtual tmp<volScalarField::Internal> Pk
        //(
        //    const volScalarField::Internal& G
        //) const;

        ////- Return epsilon/k which for standard RAS is betaStar*omega
        //virtual tmp<volScalarField::Internal> epsilonByk
        //(
        //    const volScalarField& F1,
        //    const volTensorField& gradU
        //) const;

        //- Return (G/nu)_0
        virtual tmp<volScalarField::Internal> GbyNu0
        (
            const volTensorField& gradU,
            const volScalarField& S2
        ) const;

        ////- Return G/nu
        //virtual tmp<volScalarField::Internal> GbyNu
        //(
        //    const volScalarField::Internal& GbyNu0,
        //    const volScalarField::Internal& F2,
        //    const volScalarField::Internal& S2
        //) const;

        virtual tmp<fvScalarMatrix> kSource() const;

        virtual tmp<fvScalarMatrix> omegaSource() const;

        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField::Internal& S2
            //const volScalarField::Internal& gamma,
            //const volScalarField::Internal& beta
        ) const;


public:

    //- Runtime type information
    TypeName("incompressibleKOmegaSST");


    // Constructors

        //- Construct from components
        incompressibleKOmegaSST
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~incompressibleKOmegaSST() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        ////- Return the effective diffusivity for k
        //tmp<volScalarField> DkEff(const volScalarField& F1) const
        //{
        //    return tmp<volScalarField>::New
        //    (
        //        "DkEff",
        //        alphaK(F1)*this->nut_ + this->nu()
        //    );
        //}

        ////- Return the effective diffusivity for omega
        //tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        //{
        //    return tmp<volScalarField>::New
        //    (
        //        "DomegaEff",
        //        alphaOmega(F1)*this->nut_ + this->nu()
        //    );
        //}

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "incompressibleKOmegaSSTTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
