/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "GAMGSolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::interpolate
(
    solveScalarField& psi,
    solveScalarField& Apsi,
    const lduMatrix& m,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ psiPtr = psi.begin();

    const label* const __restrict__ uPtr = m.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = m.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ diagPtr = m.diag().begin();
    const scalar* const __restrict__ upperPtr = m.upper().begin();
    const scalar* const __restrict__ lowerPtr = m.lower().begin();

    Apsi = 0;
    solveScalar* __restrict__ ApsiPtr = Apsi.begin();

    const label startRequest = UPstream::nRequests();

    m.initMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    const label nFaces = m.upper().size();
    #ifdef FOAM_OFFLOAD
    auto rfun = [](const scalar x)
    {
        return x;
    };

    matrixOperation
    (
        Apsi.begin(),
        Apsi,
        m.lduAddr(),
        matrixCoeffsMultiplyFunctor<solveScalar, scalar, decltype(rfun)>
        (
            psi,
            m.upper(),
            m.lduAddr().upperAddr(),
            rfun
        ),
        matrixCoeffsMultiplyFunctor<solveScalar, scalar, decltype(rfun)>
        (
            psi,
            m.lower(),
            m.lduAddr().lowerAddr(),
            rfun
        )
    );
    #else
    for (label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
    }
    #endif

    m.updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt,
        startRequest
    );

    const label nCells = m.diag().size();

    const labelRange cellRange(nCells);

    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            [&](const label celli)
            {
                psiPtr[celli] = -ApsiPtr[celli]/(diagPtr[celli]);
            }
        );
    });
}


void Foam::GAMGSolver::interpolate
(
    solveScalarField& psi,
    solveScalarField& Apsi,
    const lduMatrix& m,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const labelList& restrictAddressing,
    #ifdef FOAM_OFFLOAD
    const labelList& restrictSortAddressing,
    const labelList& restrictTargetAddressing,
    const labelList& restrictTargetStartAddressing,
    #endif
    const solveScalarField& psiC,
    const direction cmpt
) const
{
    interpolate
    (
        psi,
        Apsi,
        m,
        interfaceBouCoeffs,
        interfaces,
        cmpt
    );

    const label nCells = m.diag().size();
    solveScalar* __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = m.diag().begin();
    const solveScalar* const __restrict__ psiCPtr = psiC.begin();


    const label nCCells = psiC.size();
    solveScalarField corrC(nCCells, 0);
    solveScalar* __restrict__ corrCPtr = corrC.begin();

    solveScalarField diagC(nCCells, 0);
    solveScalar* __restrict__ diagCPtr = diagC.begin();

    #ifdef FOAM_OFFLOAD
    const label nTargets = restrictTargetStartAddressing.size() - 1;

    const labelRange fineCellRange(nTargets);

    apply_policy(nTargets, [&](auto policy)
    {
        std::for_each
        (
            policy,
            fineCellRange.begin(),
            fineCellRange.end(),
            [&](const label i)
            {
                const label ccelli = restrictTargetAddressing[i];

                solveScalar corr = Zero;
                solveScalar diag = Zero;

                const label start = restrictTargetStartAddressing[i];
                const label end = restrictTargetStartAddressing[i+1];

                for (label j = start; j < end; ++j)
                {
                    const label celli = restrictSortAddressing[j];

                    corr += diagPtr[celli]*psiPtr[celli];
                    diag += diagPtr[celli];
                }

                corrCPtr[ccelli] = corr;
                diagCPtr[ccelli] = diag;
            }
        );
    });


    const labelRange coarseCellRange(nCCells);

    apply_policy(nCCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            coarseCellRange.begin(),
            coarseCellRange.end(),
            [&](const label ccelli)
            {
                corrCPtr[ccelli] =
                    psiCPtr[ccelli] - corrCPtr[ccelli]/diagCPtr[ccelli];
            }
        );
    });


    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            fineCellRange.begin(),
            fineCellRange.end(),
            [&](const label celli)
            {
                const label ccelli = restrictTargetAddressing[celli];

                psiPtr[celli] += corrCPtr[ccelli];
            }
        );
    });
    #else
    for (label celli=0; celli<nCells; celli++)
    {
        corrCPtr[restrictAddressing[celli]] += diagPtr[celli]*psiPtr[celli];
        diagCPtr[restrictAddressing[celli]] += diagPtr[celli];
    }

    for (label ccelli=0; ccelli<nCCells; ccelli++)
    {
        corrCPtr[ccelli] = psiCPtr[ccelli] - corrCPtr[ccelli]/diagCPtr[ccelli];
    }

    for (label celli=0; celli<nCells; celli++)
    {
        psiPtr[celli] += corrCPtr[restrictAddressing[celli]];
    }
    #endif
}


// ************************************************************************* //
