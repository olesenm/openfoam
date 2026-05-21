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

Class
    Foam::RASModels::GEKO

Group
    grpRASTurbulence

Description
    \c GEKO implements the Menter & Matyushenko (2025) generalized k-omega
    two-equation turbulence model with the support:

    - Wall-distance-free functionality.
    - Production limiter and Kato-Launder stagnation correction.
    - Dilatation correction in the Reynolds stress tensor.
    - Machine-learning fields to support data-driven turbulence modelling.

    Reference:
    \verbatim
        Standard model (tag:MM):
            Menter, F. R., & Matyushenko, A. (2025).
            Generalized k−omega (GEKO) Two-Equation Turbulence Model.
            AIAA Journal, 63(11), 4590-4606.
    \endverbatim

Usage
    Minimal example by using \c constant/turbulenceProperties - see the
    declaration content for the meaning of the entries:
    \verbatim
    RAS
    {
        // Mandatory entries
        RASModel        GEKO;

        // Optional entries

        // Model switches
        wallDistanceFree        <bool>;
        productionLimiter       <bool>;
        katoLaunder             <bool>;
        dilatationCorrection    <bool>;
        writeCalibrationFields  <bool>;
        machineLearning         <bool>;

        // Main model coefficients for GEKO model
        Cmu             <scalar>;
        Comega2         <scalar>;
        sigmakTilde     <scalar>;
        sigmaOmegaTilde <scalar>;
        CRealize        <scalar>;

        // Auxiliary model coefficients for GEKO model
        kappa           <scalar>;
        CFbTurb         <scalar>;
        CFbLam          <scalar>;
        Aplus           <scalar>;
        Cc1             <scalar>;
        Cc2             <scalar>;
        CPKlim          <scalar>;

        // GEKO model coefficients
        Csep            <scalar>;
        Cnw             <scalar>;
        Cmix            <scalar>;
        Cjet            <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    Mandatory input fields:
    \verbatim
      Name    | Type            | Description               | Units
      k       | volScalarField  | Turbulent kinetic energy  | m^2/s^2
      omega   | volScalarField  | Specific dissipation rate | 1/s
    \endverbatim

    Optional input fields:
    \verbatim
      Name    | Type            | Description               | Units
      CSEP    | volScalarField  | Calibration for flow separation  | -
      CNW     | volScalarField  | Calibration for near-wall flows | -
      CMIX    | volScalarField  | Calibration for free-shear flows | -
      CJET    | volScalarField  | Calibration for jet flows | -
    \endverbatim

    The inherited entries are elaborated in:
      - \link eddyViscosity.H \endlink
      - \link RASModel.H \endlink

SourceFiles
    GEKO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_RASModels_GEKO_H
#define Foam_RASModels_GEKO_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class GEKO Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class GEKO
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Data

    // Model switches

        //- Flag to enable wall-distance-free model variant
        Switch wallDistanceFree_;      // (MM:p. 14; if true: F_blend = 1)

        //- Flag to disable production limiter for turbulent kinetic energy term
        Switch productionLimiter_;     // (MM:Eq. 16; default true)

        //- Flag to disable Kato-Launder correction for production term in
        //- stagnation regions
        Switch katoLaunder_;           // (MM:Eq. 15)

        //- Flag to enable dilatation correction term in the Reynolds stress
        //- tensor
        Switch dilatationCorrection_;  // (MM:Eq. A7)

        //- Flag to enable writing of calibration fields
        Switch writeCalibrationFields_;

        //- Flag to enable machine-learning fields
        //  (MM:p. 13; if true: Ck and Comega fields are read and solved)
        Switch machineLearning_;


    // Model coefficients

        // (MM:Table 3; Main model coefficients for GEKO model)
        dimensionedScalar Cmu_;
        dimensionedScalar Comega2_;
        dimensionedScalar sigmakTilde_;
        dimensionedScalar sigmaOmegaTilde_;
        dimensionedScalar CRealize_;

        // (MM:Table 4; Auxiliary model coefficients for GEKO model)
        dimensionedScalar kappa_;
        dimensionedScalar CFbTurb_;
        dimensionedScalar CFbLam_;
        dimensionedScalar Aplus_;
        dimensionedScalar Cc1_;
        dimensionedScalar Cc2_;

        // (MM:Eq. 16)
        dimensionedScalar CPKlim_;


    // GEKO model coefficients

        //- Coefficient for the CSEP field
        dimensionedScalar Csep_;

        //- Coefficient for the CNW field
        dimensionedScalar Cnw_;

        //- Coefficient for the CMIX field
        dimensionedScalar Cmix_;

        //- Coefficient for the CJET field
        dimensionedScalar Cjet_;


    // Operand Fields

        //- Turbulent kinetic energy field [m^2/s^2]
        volScalarField k_;

        //- Specific dissipation rate field [1/s]
        volScalarField omega_;


    // GEKO model fields

        //- Coefficient field to optimise flow separation from smooth surfaces
        //- such as boundary layers etc.
        volScalarField CSEP_;

        //- Coefficient field to optimise flows in non-equilibrium near-wall
        //- regions for such as heat transfer, skin friction etc.
        volScalarField CNW_;

        //- Coefficient field to optimise spreading rates for free-shear flows
        volScalarField CMIX_;

        //- Coefficient field to optimise round jet flows
        volScalarField CJET_;


    // Auxiliary fields

        //- Calibration field for flow separation
        volScalarField FSEP_;

        //- Calibration field for near-wall non-equilibrium effects
        volScalarField FNW_;

        //- Calibration field for mixing in free shear flows
        volScalarField FMIX_;

        //- Calibration field for jet flows
        volScalarField FJET_;

        //- Blending function field (MM:Eq. A8)
        volScalarField FBlend_;


    // Optional machine-learning fields (MM:p. 13)

        //- Machine-learning field to optimise turbulent kinetic energy
        autoPtr<volScalarField> CkPtr_;

        //- Machine-learning field to optimise specific dissipation rate
        autoPtr<volScalarField> ComegaPtr_;


    // Private Member Functions

        //- Validate various input fields and coefficients
        void validate_input() const;

        //- Function to update nut without providing the strain-rate magnitude
        virtual void correctNut();

        //- Function to update nut with realizability limiter
        virtual void correctNut(const volScalarField& S);


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("GEKO");


    // Generated Methods

        //- No copy construct
        GEKO(const GEKO&) = delete;

        //- No copy assignment
        void operator=(const GEKO&) = delete;


    // Constructors

        //- Construct from components
        GEKO
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~GEKO() = default;


    // Member Functions

    // Access

        //- Return the turbulence kinetic energy field [m^2/s^2]
        virtual tmp<volScalarField> k() const { return k_; }

        //- Return the specific dissipation rate field [1/s]
        virtual tmp<volScalarField> omega() const { return omega_; }


    // Evaluation

        //- Return the Reynolds stress tensor; optional dilatation corr term
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& sigmaK) const
        {
            // (MM:Eq. A1; note that sigmaK is not sigmakTilde)
            return tmp<volScalarField>::New
            (
                "DkEff",
                (this->nut_/sigmaK + this->nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& sigmaOmega) const
        {
            // (MM:Eq. A1; note that sigmaOmega is not sigmaOmegaTilde)
            return tmp<volScalarField>::New
            (
                "DomegaEff",
                (this->nut_/sigmaOmega + this->nu())
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // I-O

        //- Re-read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GEKO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
