/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

#include "comfort.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(comfort, 0);
    addToRunTimeSelectionTable(functionObject, comfort, dictionary);
}
}


// Temperature bounds based on EN ISO 7730 (10 - 40 degC)
static const Foam::scalarMinMax Tbounds(283.15, 313.15);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::functionObjects::comfort::magU() const
{
    tmp<volScalarField> tmagU = mag(lookupObject<volVectorField>("U"));
    volScalarField& magU = tmagU.ref();

    // Flag to use the averaged velocity field in the domain.
    // Consistent with EN ISO 7730 but does not make physical sense
    if (meanVelocity_)
    {
        magU = magU.weightedAverage(mesh_.V());
    }

    return tmagU;
}


void Foam::functionObjects::comfort::initWallInfo() const
{
    if (wallInfoInit_) return;

    const fvBoundaryMesh& patches = mesh_.boundary();
    label wallCount = 0;
    forAll(patches, patchi)
    {
        if (isType<wallFvPatch>(patches[patchi]))
        {
            ++wallCount;
        }
    }

    wallArea_ = 0;
    wallPatchIDs_.setSize(wallCount);

    scalar localWallArea = 0;
    label wallIndex = 0;
    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];
        if (isType<wallFvPatch>(p))
        {
            wallPatchIDs_[wallIndex++] = patchi;
            localWallArea += sum(p.magSf());
        }
    }

    wallArea_ = localWallArea;
    reduce(wallArea_, sumOp<scalar>());

    wallInfoInit_ = true;
}


Foam::dimensionedScalar Foam::functionObjects::comfort::Trad() const
{
    if (TradSet_)
    {
        if (!Tbounds.contains(Trad_.value()))
        {
            WarningInFunction
                << "The calculated mean wall radiation temperature is out of\n"
                << "the bounds specified in EN ISO 7730:2005\n"
                << "Valid range is 10 degC < T < 40 degC\n"
                << "The actual value is: " << (Trad_.value() - 273.15) << nl
                << endl;
        }
        return Trad_;
    }


    // The mean radiation is calculated by the mean wall temperatures
    // which are summed and divided by the area | only walls are taken into
    // account. This approach might be correct for a squared room but will
    // definitely be inconsistent for complex room geometries. The norm does
    // not provide any information about the calculation of this quantity.
    const volScalarField::Boundary& TBf =
        lookupObject<volScalarField>("T").boundaryField();


    scalar TareaIntegral = 0;

    forAll(wallPatchIDs_, wi)
    {
        const label patchi = wallPatchIDs_[wi];
        const fvPatchScalarField& Tp = TBf[patchi];
        const fvPatch& p = Tp.patch();
        TareaIntegral += weightedSum(p.magSf(), Tp);
    }

    reduce(TareaIntegral, sumOp<scalar>());

    dimensionedScalar TradMean
    (
        dimTemperature,
        TareaIntegral/(wallArea_ + VSMALL)
    );


    if (!Tbounds.contains(TradMean.value()))
    {
        WarningInFunction
            << "The calculated mean wall radiation temperature is out of\n"
            << "the bounds specified in EN ISO 7730:2005\n"
            << "Valid range is 10 degC < T < 40 degC\n"
            << "The actual value is: " << (TradMean.value() - 273.15) << nl
            << endl;
    }

    return TradMean;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::comfort::pSat() const
{
    static const dimensionedScalar kPaToPa(dimPressure, 1000);
    static const dimensionedScalar A(dimless, 16.6563);
    static const dimensionedScalar B(dimTemperature, 4030.183);
    static const dimensionedScalar C(dimTemperature, -38.15);

    tmp<volScalarField> tpSat;

    // Calculate the saturation pressure if no user input is given
    if (pSat_.value() == 0)
    {
        const auto& T = lookupObject<volScalarField>("T");

        // Equation based on ISO 7730:2006
        tpSat = kPaToPa*exp(A - B/(T + C));
    }
    else
    {
        tpSat = volScalarField::New("pSat", mesh_, pSat_);
    }

    return tpSat;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::comfort::Tcloth
(
    volScalarField& hc,
    const dimensionedScalar& metabolicRateSI,
    const dimensionedScalar& extWorkSI,
    const volScalarField& T,
    const dimensionedScalar& Trad
)
{
    static const dimensionedScalar factor1(dimTemperature, 308.85);
    static const dimensionedScalar factor2
    (
        dimTemperature/metabolicRateSI.dimensions(),
        0.028
    );
    static const dimensionedScalar factor3
    (
        dimMass/pow3(dimTime)/pow4(dimTemperature),
        3.96e-8
    );
    static const dimensionedScalar coeffHCf
    (
        hc.dimensions()/sqrt(dimVelocity), 12.1
    );
    static const dimensionedScalar coeffHCn
    (
        hc.dimensions()/pow025(dimTemperature), 2.38
    );
    const dimensionedScalar MminusW(metabolicRateSI - extWorkSI);
    const dimensionedScalar Trad4(pow4(Trad));

    // Heat transfer coefficient based on forced convection [W/m^2/K]
    const volScalarField hcForced(coeffHCf*sqrt(magU()));

    // Tcl [K] (surface cloth temperature)
    auto tTcl = volScalarField::New("Tcl", mesh_, dimTemperature);
    auto& Tcl = tTcl.ref();
    Tcl = T;  // initial guess
    Tcl.storePrevIter();

    // Heat transfer coefficient based on natural convection
    auto thcNatural = volScalarField::New("hcNatural", mesh_, hc.dimensions());
    auto& hcNatural = thcNatural.ref();

    label iter = 0;
    // Iterative solving of equation (2)
    do
    {
        Tcl = 0.5*(Tcl + Tcl.prevIter());
        Tcl.storePrevIter();

        // Heat transfer coefficient based on natural convection
        hcNatural = coeffHCn*pow025(mag(Tcl - T));

        // Set heat transfer coefficient based on equation (3)
        hc =
            pos(hcForced - hcNatural)*hcForced
          + neg0(hcForced - hcNatural)*hcNatural;

        // Calculate surface temperature based on equation (2)
        Tcl =
            factor1
          - factor2*MminusW
          - Icl_*factor3*fcl_*(pow4(Tcl) - Trad4)
          - Icl_*fcl_*hc*(Tcl - T);

        // Make sure that Tcl is in some physical limit (same range as we used
        // for the radiative estimation - based on ISO EN 7730:2005)
        Tcl.clamp_range(Tbounds);

    } while (!converged(Tcl) && ++iter < maxClothIter_);

    if (iter == maxClothIter_)
    {
        WarningInFunction
            << "The surface cloth temperature 'Tcl' did not converge within "
            << iter << " iterations" << nl;
    }

    return tTcl;
}


bool Foam::functionObjects::comfort::converged
(
    const volScalarField& phi
) const
{
    const auto& curr = phi.primitiveField();
    const auto& prev = phi.prevIter().primitiveField();
    forAll(curr, i)
    {
        if (mag(curr[i] - prev[i]) >= tolerance_)
        {
            return false;
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::comfort::comfort
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    clothing_("clothing", dimless, 0),
    metabolicRate_("metabolicRate", dimMass/pow3(dimTime), 0.8),
    extWork_("extWork", dimMass/pow3(dimTime), 0),
    Trad_("Trad", dimTemperature, 0),
    relHumidity_("relHumidity", dimless, 0.5),
    pSat_("pSat", dimPressure, 0),
    Icl_("Icl", pow3(dimTime)*dimTemperature/dimMass, 0),
    fcl_("fcl", dimless, 0),
    tolerance_(1e-4),
    wallArea_(0),
    maxClothIter_(100),
    wallInfoInit_(false),
    TradSet_(false),
    meanVelocity_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::comfort::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    clothing_.readIfPresent(dict);
    metabolicRate_.readIfPresent(dict);
    extWork_.readIfPresent(dict);
    pSat_.readIfPresent(dict);
    tolerance_ = dict.getOrDefault("tolerance", 1e-4);
    maxClothIter_ = dict.getOrDefault("maxClothIter", 100);
    meanVelocity_ = dict.getOrDefault("meanVelocity", false);

    // Read relative humidity if provided and convert from % to fraction
    if (dict.found(relHumidity_.name()))
    {
        relHumidity_.read(dict);
        relHumidity_ /= scalar(100);
    }

    // Read radiation temperature if provided
    if (dict.found(Trad_.name()))
    {
        TradSet_ = true;
        Trad_.read(dict);
    }

    Icl_ = dimensionedScalar(Icl_.dimensions(), 0.155)*clothing_;

    fcl_.value() =
        Icl_.value() <= 0.078
      ? 1.0 + 1.290*Icl_.value()
      : 1.05 + 0.645*Icl_.value();

    initWallInfo();

    return true;
}


bool Foam::functionObjects::comfort::execute()
{
    // Assign and build fields
    const dimensionedScalar Trad(this->Trad());
    const volScalarField pSatRH(this->pSat()*relHumidity_);

    const dimensionedScalar metabolicRateSI(58.15*metabolicRate_);
    const dimensionedScalar extWorkSI(58.15*extWork_);
    const dimensionedScalar metaDiff(metabolicRateSI - extWorkSI);

    const auto& T = lookupObject<volScalarField>("T");

    // Heat transfer coefficient [W/m^2/K]
    // This field is updated in Tcloth()
    volScalarField hc
    (
        IOobject
        (
            "hc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime)/dimTemperature, Zero)
    );

    // Calculate the surface temperature of the cloth by an iterative
    // process using equation (2) from DIN EN ISO 7730 [degC]
    const volScalarField TclothFld
    (
        this->Tcloth
        (
            hc,
            metabolicRateSI,
            extWorkSI,
            T,
            Trad
        )
    );

    // Calculate the PMV quantity
    static const dimensionedScalar factor1(pow3(dimTime)/dimMass, 0.303);
    static const dimensionedScalar factor2
    (
        dimless/metabolicRateSI.dimensions(),
        -0.036
    );
    static const dimensionedScalar factor3(factor1.dimensions(), 0.028);
    static const dimensionedScalar factor4(dimLength/dimTime, 3.05e-3);
    static const dimensionedScalar factor5(dimPressure, 5733);
    static const dimensionedScalar factor6(dimTime/dimLength, 6.99);
    static const dimensionedScalar factor8(metabolicRateSI.dimensions(), 58.15);
    static const dimensionedScalar factor9(dimless/dimPressure, 1.7e-5);
    static const dimensionedScalar factor10(dimPressure, 5867);
    static const dimensionedScalar factor11(dimless/dimTemperature, 0.0014);
    static const dimensionedScalar factor12(dimTemperature, 307.15);
    static const dimensionedScalar factor13
    (
        dimMass/pow3(dimTime)/pow4(dimTemperature),
        3.96e-8
    );

    const scalar factor7
    (
        // Special treatment of Term4
        // if metaRate - extWork < factor8, set to zero
        metaDiff.value() < factor8.value() ? 0 : 0.42
    );

    const dimensionedScalar Trad4(pow4(Trad));

    Info<< "Calculating the predicted mean vote (PMV)" << endl;

    // Equation (1)
    tmp<volScalarField> PMV =
        (
            // Term1: Thermal sensation transfer coefficient
            (factor1*exp(factor2*metabolicRateSI) + factor3)
           *(
                metaDiff

                // Term2: Heat loss difference through skin
              - factor4
               *(
                    factor5
                  - factor6*metaDiff
                  - pSatRH
                )

                // Term3: Heat loss through sweating
              - factor7*(metabolicRateSI - extWorkSI - factor8)

                // Term4: Heat loss through latent respiration
              - factor9*metabolicRateSI*(factor10 - pSatRH)

                // Term5: Heat loss through dry respiration
              - factor11*metabolicRateSI*(factor12 - T)

                // Term6: Heat loss through radiation
              - factor13*fcl_*(pow4(TclothFld) - Trad4)

                // Term7: Heat loss through convection
              - fcl_*hc*(TclothFld - T)
            )
        );

    Info<< "Calculating the predicted percentage of dissatisfaction (PPD)"
        << endl;

    // Equation (5)
    tmp<volScalarField> PPD =
        100 - 95*exp(-0.03353*pow4(PMV()) - 0.21790*sqr(PMV()));

    Info<< "Calculating the draught rating (DR)\n";

    static const dimensionedScalar Umin(dimVelocity, 0.05);
    static const dimensionedScalar Umax(dimVelocity, 0.5);
    static const dimensionedScalar pre(dimless, 0.37);
    static const dimensionedScalar C1(dimVelocity, 3.14);

    // Limit the velocity field to the values given in EN ISO 7733
    volScalarField Umag(mag(lookupObject<volVectorField>("U")));
    Umag.clamp_range(Umin, Umax);

    // Calculate the turbulent intensity if turbulent kinetic energy field k
    // exists
    volScalarField TI
    (
        IOobject
        (
            "TI",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    if (foundObject<volScalarField>("k"))
    {
        const auto& k = lookupObject<volScalarField>("k");
        TI = Foam::sqrt(2.0/3.0*k)/Umag;
    }

    // For unit correctness
    static const dimensionedScalar correctUnit
    (
        dimensionSet(0, -1.62, 1.62, -1, 0, 0, 0),
        1
    );

    // Equation (6)
    tmp<volScalarField> DR =
        correctUnit*(factor12 - T)*pow(Umag - Umin, 0.62)*(pre*Umag*TI + C1);

    // Calculate the operative temperature
    tmp<volScalarField> Top = 0.5*(T + Trad);

    // Need modifiable field names:
    word fieldNamePMV = "PMV";
    word fieldNamePPD = "PPD";
    word fieldNameDR = "DR";
    word fieldNameTop = "Top";

    return
    (
        store(fieldNamePMV, PMV)
     && store(fieldNamePPD, PPD)
     && store(fieldNameDR, DR)
     && store(fieldNameTop, Top)
    );
}


bool Foam::functionObjects::comfort::write()
{
    return
    (
        writeObject("PMV")
     && writeObject("PPD")
     && writeObject("DR")
     && writeObject("Top")
    );
}


// ************************************************************************* //
