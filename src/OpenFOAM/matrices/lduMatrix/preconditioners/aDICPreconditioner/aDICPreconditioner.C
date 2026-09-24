/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2025 Cineca
    Copyright (C) 2026 UKRI-STFC.
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

#include "aDICPreconditioner.H"
#include "tmp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aDICPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<aDICPreconditioner>
        addaDICPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aDICPreconditioner::aDICPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD_(sol.matrix().diag().size())
{
    const label nCells = rD_.size();
    const scalar* const __restrict__ diagPtr = sol.matrix().diag().cbegin();
    solveScalar* __restrict__ rDPtr = rD_.begin();

    apply_policy(nCells, [&](auto policy)
    {
        std::copy(policy, diagPtr, diagPtr + nCells, rDPtr);
    });

    calcReciprocalD(rD_, sol.matrix());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::aDICPreconditioner::calcReciprocalD
(
    solveScalarField& rD,
    const lduMatrix& matrix
)
{
    // The (untouched) original diagonal - read-only for the whole pass.
    const solveScalar* const __restrict__ diagPtr = rD.cbegin();

    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().cbegin();
    const scalar* const __restrict__ upperPtr = matrix.upper().cbegin();

    // CSR-style addressing of the faces for which a cell is the neighbour,
    // used to gather face contributions per-cell instead of scattering
    // them with an atomic add (not supported on GPUs).
    const label* const __restrict__ losortStartPtr =
        matrix.lduAddr().losortStartAddr().cbegin();
    const label* const __restrict__ losortPtr =
        matrix.lduAddr().losortAddr().cbegin();

    const label nCells = rD.size();

    // Output buffer is distinct from diagPtr: cells are processed in
    // arbitrary order and may read each other's original diagonal value
    // via lPtr, so the result cannot be written in-place until every
    // cell has finished reading.
    tmp<solveScalarField> tsumD(new solveScalarField(nCells));
    solveScalarField& sumD = tsumD.ref();
    solveScalar* __restrict__ sumDPtr = sumD.begin();

    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            labelRange(nCells).begin(),
            labelRange(nCells).end(),
            [=](const label celli)
            {
                solveScalar sum = diagPtr[celli];

                const label fStart = losortStartPtr[celli];
                const label fEnd = losortStartPtr[celli + 1];

                for (label i = fStart; i < fEnd; ++i)
                {
                    const label facei = losortPtr[i];
                    sum -= upperPtr[facei]*upperPtr[facei]/diagPtr[lPtr[facei]];
                }

                sumDPtr[celli] = sum;
            }
        );
    });

    // Calculate the reciprocal of the preconditioned diagonal
    solveScalar* __restrict__ rDPtr = rD.begin();

    apply_policy(nCells, [&](auto policy)
    {
        std::transform
        (
            policy,
            sumDPtr, sumDPtr + nCells,
            rDPtr,
            [](const solveScalar sum) { return solveScalar(1.0)/sum; }
        );
    });
}


void Foam::aDICPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction
) const
{
    solveScalar* __restrict__ wAPtr = wA.begin();
    const solveScalar* __restrict__ rAPtr = rA.cbegin();
    const solveScalar* __restrict__ rDPtr = rD_.cbegin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().cbegin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().cbegin();

    const scalar* const __restrict__ upperPtr =
        solver_.matrix().upper().cbegin();

    // CSR-style addressing: faces for which a cell is the owner (native
    // face order) and faces for which a cell is the neighbour (via the
    // losort permutation). Used to gather face contributions per-cell
    // instead of scattering them with an atomic add.
    const label* const __restrict__ ownStartPtr =
        solver_.matrix().lduAddr().ownerStartAddr().cbegin();
    const label* const __restrict__ losortStartPtr =
        solver_.matrix().lduAddr().losortStartAddr().cbegin();
    const label* const __restrict__ losortPtr =
        solver_.matrix().lduAddr().losortAddr().cbegin();

    const label nCells = wA.size();

    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            labelRange(nCells).begin(),
            labelRange(nCells).end(),
            [=](const label celli)
            {
                wAPtr[celli] = rDPtr[celli]*rAPtr[celli];
            }
        );
    });

    tmp<solveScalarField> wATmp = tmp<solveScalarField>::New(wA);
    solveScalarField& wAtmp = wATmp.ref();
    solveScalar* __restrict__ wAtmpPtr = wAtmp.begin();

    // Sweep 1: gather contributions onto each cell from the faces for
    // which it is the neighbour (reads only the owner-side wA values,
    // which are untouched during this sweep).
    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            labelRange(nCells).begin(),
            labelRange(nCells).end(),
            [=](const label celli)
            {
                solveScalar val = wAPtr[celli];

                const label fStart = losortStartPtr[celli];
                const label fEnd = losortStartPtr[celli + 1];

                for (label i = fStart; i < fEnd; ++i)
                {
                    const label facei = losortPtr[i];
                    val -= rDPtr[celli]*upperPtr[facei]*wAPtr[lPtr[facei]];
                }

                wAtmpPtr[celli] = val;
            }
        );
    });

    wA = wAtmp;

    // Sweep 2: gather contributions onto each cell from the faces for
    // which it is the owner (reads only the neighbour-side wA values,
    // as just updated by sweep 1).
    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            labelRange(nCells).begin(),
            labelRange(nCells).end(),
            [=](const label celli)
            {
                solveScalar val = wAPtr[celli];

                const label fStart = ownStartPtr[celli];
                const label fEnd = ownStartPtr[celli + 1];

                for (label facei = fStart; facei < fEnd; ++facei)
                {
                    val -= rDPtr[celli]*upperPtr[facei]*wAPtr[uPtr[facei]];
                }

                wAtmpPtr[celli] = val;
            }
        );
    });

    wA = wAtmp;
}


void Foam::aDICPreconditioner::preconditionT
(
    solveScalarField& wT,
    const solveScalarField& rT,
    const direction
) const
{
    solveScalar* __restrict__ wTPtr = wT.begin();
    const solveScalar* __restrict__ rTPtr = rT.cbegin();
    const solveScalar* __restrict__ rDPtr = rD_.cbegin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().cbegin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().cbegin();

    const scalar* const __restrict__ upperPtr =
        solver_.matrix().upper().cbegin();

    // For a symmetric matrix (this preconditioner only registers for
    // symmetric matrices) the transpose is identical to the matrix
    // itself, so preconditionT mirrors precondition() exactly.
    const label* const __restrict__ ownStartPtr =
        solver_.matrix().lduAddr().ownerStartAddr().cbegin();
    const label* const __restrict__ losortStartPtr =
        solver_.matrix().lduAddr().losortStartAddr().cbegin();
    const label* const __restrict__ losortPtr =
        solver_.matrix().lduAddr().losortAddr().cbegin();

    const label nCells = wT.size();

    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            labelRange(nCells).begin(),
            labelRange(nCells).end(),
            [=](const label celli)
            {
                wTPtr[celli] = rDPtr[celli]*rTPtr[celli];
            }
        );
    });

    tmp<solveScalarField> wTTmp = tmp<solveScalarField>::New(wT);
    solveScalarField& wTtmp = wTTmp.ref();
    solveScalar* __restrict__ wTtmpPtr = wTtmp.begin();

    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            labelRange(nCells).begin(),
            labelRange(nCells).end(),
            [=](const label celli)
            {
                solveScalar val = wTPtr[celli];

                const label fStart = losortStartPtr[celli];
                const label fEnd = losortStartPtr[celli + 1];

                for (label i = fStart; i < fEnd; ++i)
                {
                    const label facei = losortPtr[i];
                    val -= rDPtr[celli]*upperPtr[facei]*wTPtr[lPtr[facei]];
                }

                wTtmpPtr[celli] = val;
            }
        );
    });

    wT = wTtmp;

    apply_policy(nCells, [&](auto policy)
    {
        std::for_each
        (
            policy,
            labelRange(nCells).begin(),
            labelRange(nCells).end(),
            [=](const label celli)
            {
                solveScalar val = wTPtr[celli];

                const label fStart = ownStartPtr[celli];
                const label fEnd = ownStartPtr[celli + 1];

                for (label facei = fStart; facei < fEnd; ++facei)
                {
                    val -= rDPtr[celli]*upperPtr[facei]*wTPtr[uPtr[facei]];
                }

                wTtmpPtr[celli] = val;
            }
        );
    });

    wT = wTtmp;
}


// ************************************************************************* //
