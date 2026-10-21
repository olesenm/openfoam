/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "PBiCG.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PBiCG, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<PBiCG>
        addPBiCGAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PBiCG::PBiCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const scalarField& interfaceBouCoeffsData,
    const scalarField& interfaceIntCoeffsData,
    const labelUList& coeffsOffsets,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffsData,
        interfaceIntCoeffsData,
        coeffsOffsets,
        interfaces,
        solverControls
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::PBiCG::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    const label nCells = psi.size();
    const labelRange cellRange(0, nCells);

    solveScalarField pA(nCells);

    solveScalarField wA(nCells);

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffsData_, coeffsOffsets_, interfaces_, cmpt);

    // --- Calculate initial residual field
    ConstPrecisionAdaptor<solveScalar, scalar> tsource(source);
    solveScalarField rA(tsource() - wA);
    //solveScalar* __restrict__ rAPtr = rA.begin();

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        true
    );

    // --- Calculate normalisation factor
    const solveScalar normFactor = this->normFactor(psi, tsource(), wA, pA);

    if ((log_ >= 2) || (lduMatrix::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        solveScalarField pT(nCells, 0);
        //solveScalar* __restrict__ pTPtr = pT.begin();

        solveScalarField wT(nCells);
        //gsolveScalar* __restrict__ wTPtr = wT.begin();

        // --- Calculate T.psi
        matrix_.Tmul(wT, psi, interfaceIntCoeffsData_, coeffsOffsets_, interfaces_, cmpt);

        // --- Calculate initial transpose residual field
        solveScalarField rT(tsource() - wT);
        //solveScalar* __restrict__ rTPtr = rT.begin();

        // --- Initial value not used
        solveScalar wArT = 0;

        // --- Select and construct the preconditioner
        if (!preconPtr_)
        {
            preconPtr_ = lduMatrix::preconditioner::New
            (
                *this,
                controlDict_
            );
        }

        // --- Solver iteration
        do
        {
            // --- Store previous wArT
            const solveScalar wArTold = wArT;

            // --- Precondition residuals
            preconPtr_->precondition(wA, rA, cmpt);
            preconPtr_->preconditionT(wT, rT, cmpt);

            // --- Update search directions:
            wArT = gSumProd(wA, rT, matrix().mesh().comm());

            if (solverPerf.nIterations() == 0)
            {
                apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        cellRange.begin(),
                        cellRange.end(),
                        [&](const label cell)
                        {
                            pA[cell] = wA[cell];
                            pT[cell] = wT[cell];
                        }
                    );
                });
            }
            else
            {
                const solveScalar beta = wArT/wArTold;

                apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        cellRange.begin(),
                        cellRange.end(),
                        [&](const label cell)
                        {
                            pA[cell] = wA[cell] + beta*pA[cell];
                            pT[cell] = wT[cell] + beta*pT[cell];
                        }
                    );
                });
            }


            // --- Update preconditioned residuals
            matrix_.Amul(wA, pA, interfaceBouCoeffsData_, coeffsOffsets_, interfaces_, cmpt);
            matrix_.Tmul(wT, pT, interfaceIntCoeffsData_, coeffsOffsets_, interfaces_, cmpt);

            const solveScalar wApT = gSumProd(wA, pT, matrix().mesh().comm());

            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApT)/normFactor))
            {
                break;
            }


            // --- Update solution and residual:

            const solveScalar alpha = wArT/wApT;

            apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    [&](const label cell)
                    {
                        psi[cell] += alpha*pA[cell];
                        rA[cell] -= alpha*wA[cell];
                        rT[cell] -= alpha*wT[cell];
                    }
                );
            });

            solverPerf.finalResidual() =
                gSumMag(rA, matrix().mesh().comm())
               /normFactor;
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    // Recommend PBiCGStab if PBiCG fails to converge
    const label upperMaxIters = max(maxIter_, lduMatrix::defaultMaxIter);

    if (solverPerf.nIterations() > upperMaxIters)
    {
        FatalErrorInFunction
            << "PBiCG has failed to converge within the maximum iterations: "
            << upperMaxIters << nl
            << "    Please try the more robust PBiCGStab solver."
            << exit(FatalError);
    }

    if (preconPtr_)
    {
        preconPtr_->setFinished(solverPerf);
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        false
    );

    return solverPerf;
}


// ************************************************************************* //
