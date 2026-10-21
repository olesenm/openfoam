/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
    Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
    Copyright (C) 2025 Science and Technology Facility Council. All rights
                       reserved.
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

#include "twoStageGaussSeidelSmoother.H"
#include "PrecisionAdaptor.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoStageGaussSeidelSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable
    <
        twoStageGaussSeidelSmoother
    > addtwoStageGaussSeidelSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable
    <
        twoStageGaussSeidelSmoother
    > addtwoStageGaussSeidelSmootherAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoStageGaussSeidelSmoother::twoStageGaussSeidelSmoother
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
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffsData,
        interfaceIntCoeffsData,
        coeffsOffsets,
        interfaces
	//solverControls
    ),
    polynomialExpansionOrder_(solverControls.getCheckOrDefault<label>
        ("polynomialExpansionOrder", 1, labelMinMax::ge(0)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoStageGaussSeidelSmoother::smooth
(
    const word& fieldName_,
    solveScalarField& psi,
    const lduMatrix& matrix_,
    const solveScalarField& source,
    const scalarField& interfaceBouCoeffsData_,
    const labelUList& coeffsOffsets,
    const lduInterfaceFieldPtrsList& interfaces_,
    const direction cmpt,
    const label nSweeps,
    const label polynomialExpansionOrder
)
{
    const label nCells = psi.size();
    const labelRange indices(nCells);

    //solveScalarField bPrime(nCells);
    solveScalarField& bPrime = matrix_.work(nCells);

    const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();

    const label* const __restrict__ losortStartAddrPtr =
         matrix_.lduAddr().losortStartAddr().begin();

    matrix_.lowerCSR();

    //temporary field
    solveScalarField Z(psi.size()); //temporary array
    solveScalarField R(psi.size()); //residual

    solveScalar* __restrict__ rhs_ptr = bPrime.begin();
    solveScalar* __restrict__ r_ptr = R.begin();
    solveScalar* __restrict__ u_ptr = psi.begin();
    solveScalar* __restrict__ z_ptr = Z.begin();


    // Parallel boundary initialisation.  The parallel boundary is treated
    // as an effective jacobi interface in the boundary.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        bPrime = source;

        const label startRequest = UPstream::nRequests();

        matrix_.initMatrixInterfaces
        (
            false,
            interfaceBouCoeffsData_,
            coeffsOffsets,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        matrix_.updateMatrixInterfaces
        (
            false,
            interfaceBouCoeffsData_,
            coeffsOffsets,
            interfaces_,
            psi,
            bPrime,
            cmpt,
            startRequest
        );

        //r = Matrix * solution
        /******    WILL BE MORE ELEGANT IF WE HAD A FUNCITON TO COMPUTE LOCAL
        y = Matrix * x    ******/

        const label* const __restrict__ lcsrPtr =
            matrix_.lduAddr().lowerCSRAddr().begin();
        const scalar* const __restrict__ lowercsrPtr =
            matrix_.lowerCSR().begin();

        apply_policy(nCells, [&](auto policy)
        {
            std::for_each
            (
                policy,
                indices.begin(),
                indices.end(),
                [&](label cell)
                {
                    auto val = diagPtr[cell] * u_ptr[cell];

                    // Lower contributions
                    {
                        const label start = losortStartAddrPtr[cell];
                        const label end = losortStartAddrPtr[cell + 1];
                        #pragma unroll 2
                        for (label i = start; i < end; ++i)
                        {
                            val += lowercsrPtr[i] * u_ptr[lcsrPtr[i]];
                        }
                    }

                    // Upper contributions
                    {
                        const label start = ownStartPtr[cell];
                        const label end = ownStartPtr[cell + 1];
                        #pragma unroll 4
                        for (label i = start; i < end; ++i)
                        {
                            val += upperPtr[i] * u_ptr[uPtr[i]];
                        }
                    }

                    r_ptr[cell] = val;
                }
            );
        });

        if (0 == polynomialExpansionOrder)
        {
            scalar relax_weight = 0.8;

            // 1) z = r/D, u = u + relax_weight * z
            apply_policy(nCells, [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    indices.begin(),
                    indices.end(),
                    [&](size_t i)
                    {
                        scalar r = rhs_ptr[i] - r_ptr[i];
                        u_ptr[i] += relax_weight * r / diagPtr[i];
                    }
                );
            });
        }
        else
        {
            // 1) z = r/D, u = u + z
            apply_policy(nCells, [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    indices.begin(),
                    indices.end(),
                    [&](size_t i)
                    {
                        scalar r = rhs_ptr[i] - r_ptr[i]; // residual
                        z_ptr[i] = r / diagPtr[i];
                        u_ptr[i] += z_ptr[i];
                    }
                );
            });

            scalar multiplier = -1.0;
            for (label sweepID = 0; sweepID < polynomialExpansionOrder;
                sweepID++)
            {
                // 2) r = U * z
                apply_policy(nCells, [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        indices.begin(),
                        indices.end(),
                        [&](label celli)
                        {
                            const label start = ownStartPtr[celli];
                            const label end   = ownStartPtr[celli + 1];

                            scalar tmp = 0.0;
                            #pragma unroll 4
                            for (label facei = start; facei < end; ++facei)
                            {
                                tmp += upperPtr[facei] * z_ptr[uPtr[facei]];
                            }

                            r_ptr[celli] = tmp;
                        }
                    );
                });

                // 3) z = r/D, u = u + m * z
                if (sweepID < polynomialExpansionOrder-1)
                {
                    apply_policy(nCells, [&](auto policy)
                    {
                        std::for_each
                        (
                            policy,
                            indices.begin(),
                            indices.end(),
                            [&](label celli)
                            {
                                z_ptr[celli] = r_ptr[celli] / diagPtr[celli];
                                u_ptr[celli] += multiplier * z_ptr[celli];
                            }
                        );
                    });
                }
                else
                {
                    //no need to save "Z"
                    apply_policy(nCells, [&](auto policy)
                    {
                        std::for_each
                        (
                            policy,
                            indices.begin(), indices.end(),
                            [&](label celli)
                            {
                                u_ptr[celli] += multiplier*r_ptr[celli]
                                    /diagPtr[celli];
                            }
                        );
                    });
                }
                multiplier *= -1.0;
            }
        }
    }
}


void Foam::twoStageGaussSeidelSmoother::smooth
(
    solveScalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        fieldName_,
        psi,
        matrix_,
        ConstPrecisionAdaptor<solveScalar, scalar>(source),
        interfaceBouCoeffsData_,
        coeffsOffsets_,
        interfaces_,
        cmpt,
        nSweeps,
	    polynomialExpansionOrder_
    );
}


void Foam::twoStageGaussSeidelSmoother::scalarSmooth
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        fieldName_,
        psi,
        matrix_,
        source,
        interfaceBouCoeffsData_,
        coeffsOffsets_,
        interfaces_,
        cmpt,
        nSweeps,
	    polynomialExpansionOrder_
    );
}


// ************************************************************************* //
