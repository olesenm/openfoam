/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "PBiCGStab.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PBiCGStab, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PBiCGStab>
        addPBiCGStabSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<PBiCGStab>
        addPBiCGStabAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PBiCGStab::PBiCGStab
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::PBiCGStab::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    const label nCells = psi.size();
    const labelRange cellRange(0, nCells);

    solveScalarField pA(nCells);
    solveScalarField yA(nCells);

    // --- Calculate A.psi
    matrix_.Amul(yA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    solveScalarField rA(source.size());

    apply_policy(static_cast<std::size_t>(source.size()), [&](auto policy)
    {
        std::for_each
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            [&](const label cell)
            {
                rA[cell] = source[cell] - yA[cell];
            }
        );
    });

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        true
    );

    // --- Calculate normalisation factor
    const solveScalar normFactor = this->normFactor(psi, source, yA, pA);

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
        solveScalarField AyA(nCells);

        solveScalarField sA(nCells);

        solveScalarField zA(nCells);

        solveScalarField tA(nCells);

        // --- Store initial residual
        const solveScalarField rA0(rA);

        // --- Initial values not used
        solveScalar rA0rA = 0;
        solveScalar alpha = 0;
        solveScalar omega = 0;

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
            // --- Store previous rA0rA
            const solveScalar rA0rAold = rA0rA;

            rA0rA = gSumProd(rA0, rA, matrix().mesh().comm());

            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(rA0rA)))
            {
                break;
            }

            // --- Update pA
            if (solverPerf.nIterations() == 0)
            {
                apply_policy(rA.size(), [&](auto policy)
                {
                    std::copy
                    (
                        policy,
                        rA.begin(),
                        rA.end(),
                        pA.begin()
                    );
                });
            }
            else
            {
                // --- Test for singularity
                if (solverPerf.checkSingularity(mag(omega)))
                {
                    break;
                }

                const solveScalar beta = (rA0rA/rA0rAold)*(alpha/omega);

                apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        cellRange.begin(),
                        cellRange.end(),
                        [&](const label cell)
                        {
                            pA[cell] =
                                rA[cell] + beta*(pA[cell] - omega*AyA[cell]);
                        }
                    );
                });
            }

            // --- Precondition pA
            preconPtr_->precondition(yA, pA, cmpt);

            // --- Calculate AyA
            matrix_.Amul(AyA, yA, interfaceBouCoeffs_, interfaces_, cmpt);

            const solveScalar rA0AyA =
                gSumProd(rA0, AyA, matrix().mesh().comm());

            alpha = rA0rA/rA0AyA;

            // --- Calculate sA
            apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    [&](const label cell)
                    {
                        sA[cell] = rA[cell] - alpha*AyA[cell];
                    }
                );
            });

            // --- Test sA for convergence
            solverPerf.finalResidual() =
                gSumMag(sA, matrix().mesh().comm())/normFactor;

            if
            (
                solverPerf.nIterations() >= minIter_
             && solverPerf.checkConvergence(tolerance_, relTol_, log_)
            )
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
                            psi[cell] += alpha*yA[cell];
                        }
                    );
                }); 

                solverPerf.nIterations()++;

                return solverPerf;
            }

            // --- Precondition sA
            preconPtr_->precondition(zA, sA, cmpt);

            // --- Calculate tA
            matrix_.Amul(tA, zA, interfaceBouCoeffs_, interfaces_, cmpt);

            const solveScalar tAtA = gSumSqr(tA, matrix().mesh().comm());

            // --- Calculate omega from tA and sA
            //     (cheaper than using zA with preconditioned tA)
            omega = gSumProd(tA, sA, matrix().mesh().comm())/tAtA;

            // --- Update solution and residual
            apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    [&](const label cell)
                    {
                        psi[cell] += alpha*yA[cell] + omega*zA[cell];
                        rA[cell] = sA[cell] - omega*tA[cell];
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


Foam::solverPerformance Foam::PBiCGStab::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}


// ************************************************************************* //
