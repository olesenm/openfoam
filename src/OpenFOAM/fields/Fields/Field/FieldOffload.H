/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Header for template-like definitions for Field operations that can
    be offloaded to GPU.

SourceFiles
    FieldOffload.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_FieldOffload_H
#define Foam_FieldOffload_H

#include "UList.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type>
    Type maxOffload(const UList<Type>& f1);

    template<class Type>
    Type minOffload(const UList<Type>& f1);

    template<class Type>
    typename typeOfSolve<Type>::type
    sumOffload(const UList<Type>& f);

    template<class Type>
    Type maxMagSqrOffload(const UList<Type>& f1);

    template<class Type>
    Type minMagSqrOffload(const UList<Type>& f1);

    template<class Type>
    typename scalarProduct<Type, Type>::type
    sumProdOffload(const UList<Type>& f1, const UList<Type>& f2);

    float sumProdScalarOffload(const UList<float>& f1, const UList<float>& f2);
    double sumProdScalarOffload(const UList<double>& f1, const UList<double>& f2);

    template<class Type>
    Type sumCmptProdOffload(const UList<Type>& f1, const UList<Type>& f2);

    template<class Type>
    typename outerProduct1<Type>::type
    sumSqrOffload(const UList<Type>& f1);

    template<class Type>
    typename typeOfMag<Type>::type
    sumMagOffload(const UList<Type>& f1);

    template<class Type>
    Type sumCmptMagOffload(const UList<Type>& f1);

    template<class Type>
    void gWeightedAverageOffload
    (const UList<scalar>& weights, const UList<Type>& fld, scalar& weight, Type& result);

    template<class Type>
    void gWeightedSumOffload
    (const UList<scalar>& weights, const UList<Type>& fld, Type& result);
}

#endif

// ************************************************************************* //

