/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "kinematicThinFilmIBM.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kinematicThinFilmIBM, 0);
addToRunTimeSelectionTable(liquidFilmBase, kinematicThinFilmIBM, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicThinFilmIBM::kinematicThinFilmIBM
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    kinematicThinFilm(modelType, mesh, dict),
    immersedBoundary_(primaryMesh().time(), regionMesh())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void kinematicThinFilmIBM::preEvolveRegion()
{
    kinematicThinFilm::preEvolveRegion();

    // Update the IBM models
    immersedBoundary_.correct(Uf_);
}


void kinematicThinFilmIBM::evolveRegion()
{
    DebugInFunction << endl;

    auto& Ustar = immersedBoundary_.Ustar();

    phi2s_ = fac::interpolate(h_)*phif_;

    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        // After the first one, outer corrections use updated h_, so mass
        // source in momentum eq should be omitted in remaining outer iters.
        const bool includeMassSource = (includeMassSourceOuter_ && oCorr == 1);

        areaVectorField rhs(-fac::grad(pf_*h_)/rho_ + pf_*fac::grad(h_)/rho_);

        // IBM force
        areaVectorField f
        (
            (UEqn(Ustar, includeMassSource) - rhs) & Ustar
        );

        // Ensure force is only applied to IBM faces
        immersedBoundary_.zeroFilter(f);

        pf_.storePrevIter();

        faVectorMatrix UsEqn(UEqn(Uf_, includeMassSource) == f);

        UsEqn.relax();

        faOptions().constrain(UsEqn);

        if (momentumPredictor_)
        {
            solve(UsEqn == rhs);
        }

        for (int corr=1; corr<=nCorr_; corr++)
        {
            areaScalarField UsA(UsEqn.A());

            Uf_ = UsEqn.H()/UsA;
            Uf_.correctBoundaryConditions();
            faOptions().correct(Uf_);

            phif_ =
                (fac::interpolate(Uf_) & regionMesh().Le())
                - fac::interpolate(1.0/(rho_*UsA))
                * fac::lnGrad(pf_*h_)*regionMesh().magLe()
                + fac::interpolate(pf_/(rho_*UsA))
                * fac::lnGrad(h_)*regionMesh().magLe();

            for (int nFilm=1; nFilm<=nFilmCorr_; nFilm++)
            {
                faScalarMatrix hEqn
                (
                    fam::ddt(h_)
                  + fam::div(phif_, h_)
                 ==
                    faOptions()(rho_, h_, dimVelocity)
                  + rhoSp_
                );

                hEqn.relax();
                faOptions().constrain(hEqn);
                hEqn.solve();
                faOptions().correct(h_);

                if (nFilm == nFilmCorr_)
                {
                    phi2s_ = hEqn.flux();
                }
            }

            // Bound h_
            h_ = max(h_, h0_);

            pf_ = rho_*gn_*h_ - sigma_*fac::laplacian(h_) + pnSp_ + ppf_;
            pf_.correctBoundaryConditions();
            pf_.relax();

            Uf_ -= (1.0/(rho_*UsA))*fac::grad(pf_*h_)
                - (pf_/(rho_*UsA))*fac::grad(h_);
            Uf_.correctBoundaryConditions();
            faOptions().correct(Uf_);
        }
    }

    Info<< "Film h min/max   = " << gMinMax(h_) << nl
        << "Film mag(U) min/max   = " << gMinMaxMag(Uf_) << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
