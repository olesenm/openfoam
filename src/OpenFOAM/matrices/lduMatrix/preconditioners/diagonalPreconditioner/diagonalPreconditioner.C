/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "diagonalPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diagonalPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<diagonalPreconditioner>
        adddiagonalPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<diagonalPreconditioner>
        adddiagonalPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diagonalPreconditioner::diagonalPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD(sol.matrix().diag().size())
{
    const auto& diag = solver_.matrix().diag();

    apply_policy(static_cast<std::size_t>(diag.size()), [&](auto policy)
    {
        std::transform
        (
            policy,
            diag.cbegin(),
            diag.cend(),
            rD.begin(),
            [](const solveScalar value)
            {
                return solveScalar(1)/value;
            }
        );
    });
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diagonalPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction
) const
{
    apply_policy(static_cast<std::size_t>(rD.size()), [&](auto policy)
    {
        std::transform
        (
            policy,
            rD.cbegin(),
            rD.cend(),
            rA.cbegin(),
            wA.begin(),
            [](const solveScalar dValue, const solveScalar aValue)
            {
                return dValue*aValue;
            }
        );
    });
}


// ************************************************************************* //
