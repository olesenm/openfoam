/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
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

Class
    Foam::twoStageGaussSeidelSmoother

Group
    grpLduMatrixSmoothers

Description
    A lduMatrix::smoother for Two Stage Gauss-Seidel

    https://www.researchgate.net/publication/356375724_Scaled_Smoothers_for_Navier-Stokes_Pressure_Projection
    DOI:10.48550/arXiv.2111.09512
    Scaled Smoothers for Navier-Stokes Pressure Projection

SourceFiles
    twoStageGaussSeidelSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef twoStageGaussSeidelSmoother_H
#define twoStageGaussSeidelSmoother_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class twoStageGaussSeidelSmoother Declaration
\*---------------------------------------------------------------------------*/

class twoStageGaussSeidelSmoother
:
    public lduMatrix::smoother
{

private:
    label polynomialExpansionOrder_;

public:

    //- Runtime type information
    TypeName("twoStageGaussSeidel");

    // Constructors

        //- Construct from components
        twoStageGaussSeidelSmoother
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
	        const dictionary& solverControls
        );


    // Member Functions

        //- Smooth for the given number of sweeps
        static void smooth
        (
            const word& fieldName,
            solveScalarField& psi,
            const lduMatrix& matrix,
            const solveScalarField& source,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const direction cmpt,
            const label nSweeps,
	        const label polynomialExpansionOrder
        );


        //- Smooth the solution for a given number of sweeps
        virtual void smooth
        (
            solveScalarField& psi,
            const scalarField& Source,
            const direction cmpt,
            const label nSweeps
        ) const;

        //- Smooth the solution for a given number of sweeps
        virtual void scalarSmooth
        (
            solveScalarField& psi,
            const solveScalarField& Source,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
