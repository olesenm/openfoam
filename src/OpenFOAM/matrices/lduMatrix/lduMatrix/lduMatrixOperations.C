/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    lduMatrix member operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

#ifdef FOAM_OFFLOAD
#include "lduMatrixFunctors.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::sumDiag()
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
    scalarField& Diag = diag();

    #ifdef FOAM_OFFLOAD
    matrixOperation
    (
        Diag.begin(),
        Diag,
        lduAddr(),
        matrixCoeffsFunctor<scalar, unityUnaryOperatorFunctor<scalar, scalar>>
        (
            Lower,
            unityUnaryOperatorFunctor<scalar, scalar>()
        ),
        matrixCoeffsFunctor<scalar, unityUnaryOperatorFunctor<scalar, scalar>>
        (
            Upper,
            unityUnaryOperatorFunctor<scalar, scalar>()
        )
    );
    #else
    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face=0; face<l.size(); face++)
    {
        Diag[l[face]] += Lower[face];
        Diag[u[face]] += Upper[face];
    }
    #endif
}


void Foam::lduMatrix::negSumDiag()
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
    scalarField& Diag = diag();
    #ifdef FOAM_OFFLOAD
    matrixOperation
    (
        Diag.begin(),
        Diag,
        lduAddr(),
        matrixCoeffsFunctor<scalar, negateUnaryOperatorFunctor<scalar, scalar>>
        (
            Lower,
            negateUnaryOperatorFunctor<scalar, scalar>()
        ),
        matrixCoeffsFunctor<scalar, negateUnaryOperatorFunctor<scalar, scalar>>
        (
            Upper,
            negateUnaryOperatorFunctor<scalar, scalar>()
        )
    );
    #else
    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face=0; face<l.size(); face++)
    {
        Diag[l[face]] -= Lower[face];
        Diag[u[face]] -= Upper[face];
    }
    #endif
}


void Foam::lduMatrix::sumMagOffDiag
(
    scalarField& sumOff
) const
{
    #ifdef FOAM_OFFLOAD
    matrixOperation
    (
        sumOff.begin(),
        sumOff,
        lduAddr(),
        matrixCoeffsFunctor<scalar, magUnaryFunctionFunctor<scalar, scalar>>
        (
            upper(),
            magUnaryFunctionFunctor<scalar,scalar>()
        ),
        matrixCoeffsFunctor<scalar, magUnaryFunctionFunctor<scalar, scalar>>
        (
            lower(),
            magUnaryFunctionFunctor<scalar, scalar>()
        )
    );
    #else
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face = 0; face < l.size(); face++)
    {
        sumOff[u[face]] += mag(Lower[face]);
        sumOff[l[face]] += mag(Upper[face]);
    }
    #endif
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::lduMatrix::operator=(const lduMatrix& A)
{
    if (this == &A)
    {
        return;  // Self-assignment is a no-op
    }

    if (A.hasLower())
    {
        lower() = A.lower();
    }
    else
    {
        lowerPtr_.reset(nullptr);
    }

    if (A.hasUpper())
    {
        upper() = A.upper();
    }
    else
    {
        upperPtr_.reset(nullptr);
    }

    if (A.hasDiag())
    {
        diag() = A.diag();
    }
}


void Foam::lduMatrix::operator=(lduMatrix&& A)
{
    if (this == &A)
    {
        return;  // Self-assignment is a no-op
    }

    diagPtr_ = std::move(A.diagPtr_);
    upperPtr_ = std::move(A.upperPtr_);
    lowerPtr_ = std::move(A.lowerPtr_);
}


void Foam::lduMatrix::negate()
{
    if (diagPtr_)
    {
        diagPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }
}


void Foam::lduMatrix::operator+=(const lduMatrix& A)
{
    if (A.hasDiag())
    {
        diag() += A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        upper() += A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() += A.upper();
        lower() += A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.hasUpper())
        {
            lower() += A.upper();
            upper() += A.upper();
        }
        else
        {
            lower() += A.lower();
            upper() += A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() += A.lower();
        upper() += A.upper();
    }
    else if (diagonal())
    {
        if (A.hasUpper())
        {
            upper() = A.upper();
        }

        if (A.hasLower())
        {
            lower() = A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        if (debug > 1)
        {
            WarningInFunction
                << "Unknown matrix type combination" << nl
                << "    this : " << this->matrixTypeName()
                << "    A    : " << A.matrixTypeName() << endl;
        }
    }
}


void Foam::lduMatrix::operator-=(const lduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() -= A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        upper() -= A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() -= A.upper();
        lower() -= A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.hasUpper())
        {
            lower() -= A.upper();
            upper() -= A.upper();
        }
        else
        {
            lower() -= A.lower();
            upper() -= A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() -= A.lower();
        upper() -= A.upper();
    }
    else if (diagonal())
    {
        if (A.hasUpper())
        {
            upper() = -A.upper();
        }

        if (A.hasLower())
        {
            lower() = -A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        if (debug > 1)
        {
            WarningInFunction
                << "Unknown matrix type combination" << nl
                << "    this : " << this->matrixTypeName()
                << "    A    : " << A.matrixTypeName() << endl;
        }
    }
}


void Foam::lduMatrix::operator*=(const scalarField& sf)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    // Non-uniform scaling causes a symmetric matrix
    // to become asymmetric
    if (symmetric() || asymmetric())
    {
        scalarField& upper = this->upper();
        scalarField& lower = this->lower();

        const labelUList& l = lduAddr().lowerAddr();
        const labelUList& u = lduAddr().upperAddr();

        const labelRange upperFaces(upper.size());

        apply_policy(static_cast<std::size_t>(upper.size()), [&](auto policy)
        {
            std::for_each
            (
                policy,
                upperFaces.begin(),
                upperFaces.end(),
                [&](const label face)
                {
                    upper[face] *= sf[l[face]];
                }
            );
        });

        const labelRange lowerFaces(lower.size());

        apply_policy(static_cast<std::size_t>(lower.size()), [&](auto policy)
        {
            std::for_each
            (
                policy,
                lowerFaces.begin(),
                lowerFaces.end(),
                [&](const label face)
                {
                    lower[face] *= sf[u[face]];
                }
            );
        });
    
    }
}


void Foam::lduMatrix::operator*=(scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }
}


// ************************************************************************* //
