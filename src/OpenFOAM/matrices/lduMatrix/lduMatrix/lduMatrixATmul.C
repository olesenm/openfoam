/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2026 OpenCFD Ltd.
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

    Some parts here added have been inspired by the work done at the following
    repository https://github.com/SimFlowCFD/RapidCFD-dev

Description
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

#ifdef FOAM_OFFLOAD
#include "lduMatrixFunctors.H"
#include "lduAddressingFunctors.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef FOAM_OFFLOAD

namespace Foam
{

template<int nUnroll>
struct matrixMultiplyFunctor
{
    const solveScalarField& psi_;
    const scalarField& diag_;
    const scalarField& lower_;
    const scalarField& upper_;

    const labelUList& own_;
    const labelUList& nei_;
    const labelUList& ownStart_;
    const labelUList& losortStart_;
    const labelUList& losort_;

    matrixMultiplyFunctor
    (
        const solveScalarField& psi,
        const scalarField& diag,
        const scalarField& lower,
        const scalarField& upper,
        const labelUList& own,
        const labelUList& nei,
        const labelUList& ownStart,
        const labelUList& losortStart,
        const labelUList& losort
    )
    :
        psi_(psi),
        diag_(diag),
        lower_(lower),
        upper_(upper),
        own_(own),
        nei_(nei),
        ownStart_(ownStart),
        losortStart_(losortStart),
        losort_(losort)
    {}

    solveScalar operator()(const label celli) const
    {
        solveScalar tmpSum[2*nUnroll] = {Zero};
        solveScalar nbrSum = Zero;

        const label ownStart = ownStart_[celli];
        const label ownSize = ownStart_[celli + 1] - ownStart;

        const label neiStart = losortStart_[celli];
        const label neiSize = losortStart_[celli + 1] - neiStart;

        solveScalar result = diag_[celli]*psi_[celli];

        for (label i = 0; i < nUnroll; ++i)
        {
            if (i < ownSize)
            {
                const label facei = ownStart + i;
                tmpSum[i] = upper_[facei]*psi_[nei_[facei]];
            }
        }

        for (label i = 0; i < nUnroll; ++i)
        {
            if (i < neiSize)
            {
                const label facei = losort_[neiStart + i];
                tmpSum[i + nUnroll] = lower_[facei]*psi_[own_[facei]];
            }
        }

        for (label i = 0; i < 2*nUnroll; ++i)
        {
            result += tmpSum[i];
        }

        for (label i = nUnroll; i < ownSize; ++i)
        {
            const label facei = ownStart + i;
            result += upper_[facei]*psi_[nei_[facei]];
        }

        for (label i = nUnroll; i < neiSize; ++i)
        {
            const label facei = losort_[neiStart + i];
            nbrSum += lower_[facei]*psi_[own_[facei]];
        }

        return result + nbrSum;
    }
};

} // namespace Foam
#endif

void Foam::lduMatrix::Amul
(
    solveScalarField& Apsi,
    const tmp<solveScalarField>& tpsi,
    const solveScalarField& interfaceBouCoeffsData,
    const labelUList& offsets,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    const auto& addr = lduAddr();

    solveScalar* __restrict__ ApsiPtr = Apsi.begin();

    const solveScalarField& psi = tpsi();
    const solveScalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = addr.upperAddr().begin();
    const label* const __restrict__ lPtr = addr.lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfaceBouCoeffsData,
        offsets,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    const label nCells = diag().size();

    if (hasLowerCSR())
    {
        // Use cell-based looping
        if (debug == 2) PoutInFunction<< "cell-based looping" << endl;

        const label* const __restrict__ oStartPtr =
            addr.ownerStartAddr().begin();
        const label* const __restrict__ loStartPtr =
            addr.losortStartAddr().begin();
            const label* const __restrict__ lcsrPtr =
                addr.lowerCSRAddr().begin();

        // Note: lowerCSR constructed from lower if available, upper otherwise
        //       so is handling symmetric()
        const scalar* const __restrict__ lowercsrPtr = lowerCSR().begin();

        const labelRange cellIndices(nCells);

        apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
        {
            std::for_each
            (
                policy,
                cellIndices.begin(),
                cellIndices.end(),
                [&](const label cell)
                {
                    auto& val = ApsiPtr[cell];

                    val = diagPtr[cell]*psiPtr[cell];

                    // Add lower contributions
                    {
                        const label start = loStartPtr[cell];
                        const label end = loStartPtr[cell+1];

                        for (label i = start; i < end; i++)
                        {
                            const label nbrCell = lcsrPtr[i];
                            val += lowercsrPtr[i]*psiPtr[nbrCell];
                        }
                    }
                    // Add upper contributions
                    {
                        const label start = oStartPtr[cell];
                        const label end = oStartPtr[cell+1];

                        for (label i = start; i < end; i++)
                        {
                            const label nbrCell = uPtr[i];
                            val += upperPtr[i]*psiPtr[nbrCell];
                        }
                    }
                }
            );
        });
    }
    else
    {
        #ifdef FOAM_OFFLOAD
        const labelRange indices(psi.size());

        apply_policy(static_cast<std::size_t>(psi.size()), [&](auto policy)
        {
            std::transform
            (
                policy,
                indices.begin(),
                indices.end(),
                Apsi.begin(),
                matrixMultiplyFunctor<3>
                (
                    tpsi,
                    diag(),
                    lower(),
                    upper(),
                    lduAddr().lowerAddr(),
                    lduAddr().upperAddr(),
                    lduAddr().ownerStartAddr(),
                    lduAddr().losortStartAddr(),
                    lduAddr().losortAddr()
                )
            );
        });
        #else
        for (label cell=0; cell<nCells; cell++)
        {
            ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
        }

        const label nFaces = lduAddr().lowerAddr().size();
        for (label face=0; face<nFaces; face++)
        {
            ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
            ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
        }
        #endif
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffsData,
        offsets,
        interfaces,
        psi,
        Apsi,
        cmpt,
        startRequest
    );

    tpsi.clear();
}


void Foam::lduMatrix::Tmul
(
    solveScalarField& Tpsi,
    const tmp<solveScalarField>& tpsi,
    const solveScalarField& interfaceBouCoeffsData,
    const labelUList& offsets,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ TpsiPtr = Tpsi.begin();

    const solveScalarField& psi = tpsi();
    const solveScalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ lowerPtr = lower().begin();
    const scalar* const __restrict__ upperPtr = upper().begin();

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfaceBouCoeffsData,
        offsets,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    const label nCells = diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        TpsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    const label nFaces = lduAddr().lowerAddr().size();
    for (label face=0; face<nFaces; face++)
    {
        TpsiPtr[uPtr[face]] += upperPtr[face]*psiPtr[lPtr[face]];
        TpsiPtr[lPtr[face]] += lowerPtr[face]*psiPtr[uPtr[face]];
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffsData,
        offsets,
        interfaces,
        psi,
        Tpsi,
        cmpt,
        startRequest
    );

    tpsi.clear();
}


void Foam::lduMatrix::sumA
(
    solveScalarField& sumA,
    const solveScalarField& interfaceBouCoeffsData,
    const labelUList& offsets,
    const lduInterfaceFieldPtrsList& interfaces
) const
{
    solveScalar* __restrict__ sumAPtr = sumA.begin();

    const scalar* __restrict__ diagPtr = diag().begin();

    const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* __restrict__ lowerPtr = lower().begin();
    const scalar* __restrict__ upperPtr = upper().begin();

    const label nCells = diag().size();
    const label nFaces = lduAddr().lowerAddr().size();

    #ifdef FOAM_OFFLOAD
    matrixOperation
    (
        diag().begin(),
        sumA,
        lduAddr(),
        matrixCoeffsFunctor<scalar, unityOp<scalar>>
        (
            upper(),
            unityOp<scalar>()
        ),
        matrixCoeffsFunctor<scalar, unityOp<scalar>>
        (
            lower(),
            unityOp<scalar>()
        )
    );
    #else
    for (label cell=0; cell<nCells; cell++)
    {
        sumAPtr[cell] = diagPtr[cell];
    }

    for (label face=0; face<nFaces; face++)
    {
        sumAPtr[uPtr[face]] += lowerPtr[face];
        sumAPtr[lPtr[face]] += upperPtr[face];
    }
    #endif

    // Subtract interface boundary coefficients from row sum
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            const labelUList& pa = lduAddr().patchAddr(patchi);
            const label start = offsets[patchi];
            const SubField<solveScalar> pCoeffs(interfaceBouCoeffsData, pa.size(), start);

            const labelRange faceIndices(pa.size());

            apply_policy(static_cast<std::size_t>(pa.size()), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    faceIndices.begin(),
                    faceIndices.end(),
                    [&](const auto& face)
                    {
                        sumAPtr[pa[face]] -= pCoeffs[face];
                    }
                );
            });
        }
    }
}


void Foam::lduMatrix::residual
(
    solveScalarField& rA,
    const solveScalarField& psi,
    const scalarField& source,
    const solveScalarField& interfaceBouCoeffsData,
    const labelUList& offsets,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ rAPtr = rA.begin();

    const solveScalar* const __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = diag().begin();
    const scalar* const __restrict__ sourcePtr = source.begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        false,
        interfaceBouCoeffsData,
        offsets,
        interfaces,
        psi,
        rA,
        cmpt
    );

    const label nCells = diag().size();

    #ifdef FOAM_OFFLOAD
    const labelRange cellIndices(psi.size());

    apply_policy(static_cast<std::size_t>(nCells), [&](auto policy)
    {
        std::transform
        (
            policy,
            cellIndices.begin(),
            cellIndices.end(),
            rAPtr,
            [&](const label i)
            {
                return source[i] - diagPtr[i]*psi[i];
            }
        );
    });

    //const label nFaces = upper().size();
    matrixOperation
    (
        rA.begin(),
        rA,
        lduAddr(),
        matrixCoeffsMultiplyFunctor<solveScalar, scalar,
        negateUnaryOperatorFunctor<solveScalar, scalar>>
        (
            psi,
            upper(),
            lduAddr().upperAddr(),
            negateUnaryOperatorFunctor<solveScalar,scalar>()
        ),
        matrixCoeffsMultiplyFunctor<solveScalar, scalar,
        negateUnaryOperatorFunctor<solveScalar, scalar>>
        (
            psi,
            lower(),
            lduAddr().lowerAddr(),
            negateUnaryOperatorFunctor<solveScalar,scalar>()
        )
    );
    #else
    for (label cell=0; cell<nCells; cell++)
    {
        rAPtr[cell] = sourcePtr[cell] - diagPtr[cell]*psiPtr[cell];
    }

    const label nFaces = lduAddr().lowerAddr().size();

    for (label face=0; face<nFaces; face++)
    {
        rAPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
        rAPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
    }
    #endif

    // Update interface interfaces
    updateMatrixInterfaces
    (
        false,
        interfaceBouCoeffsData,
        offsets,
        interfaces,
        psi,
        rA,
        cmpt,
        startRequest
    );
}


Foam::tmp<Foam::Field<Foam::solveScalar>> Foam::lduMatrix::residual
(
    const solveScalarField& psi,
    const scalarField& source,
    const solveScalarField& interfaceBouCoeffsData,
    const labelUList& offsets,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    auto trA = tmp<solveScalarField>::New(psi.size());

    residual
    (
        trA.ref(),
        psi,
        source,
        interfaceBouCoeffsData,
        offsets,
        interfaces,
        cmpt
    );

    return trA;
}


Foam::tmp<Foam::scalarField> Foam::lduMatrix::H1() const
{
    auto tH1 = tmp<scalarField>::New(lduAddr().size(), Zero);

    if (lowerPtr_ || upperPtr_)
    {
        scalar* __restrict__ H1Ptr = tH1.ref().begin();

        const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* __restrict__ lowerPtr = lower().begin();
        const scalar* __restrict__ upperPtr = upper().begin();

        const label nFaces = lduAddr().lowerAddr().size();

        #ifdef FOAM_OFFLOAD
        matrixOperation
        (
            H1Ptr,
            tH1.ref(),
            lduAddr(),
            matrixCoeffsFunctor<scalar, negateUnaryOperatorFunctor<scalar,
                scalar>>
            (
                upper(),
                negateUnaryOperatorFunctor<scalar, scalar>()
            ),
            matrixCoeffsFunctor<scalar, negateUnaryOperatorFunctor<scalar,
                scalar>>
            (
                lower(),
                negateUnaryOperatorFunctor<scalar, scalar>()
            )
        );
        #else
        for (label face=0; face<nFaces; face++)
        {
            H1Ptr[uPtr[face]] -= lowerPtr[face];
            H1Ptr[lPtr[face]] -= upperPtr[face];
        }
        #endif
    }

    return tH1;
}


// ************************************************************************* //
