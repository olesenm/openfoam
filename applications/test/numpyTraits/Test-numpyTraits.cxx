/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "pTraits.H"
#include "contiguous.H"
#include "boolVector.H"  // A FixedList pretending to be a vector
#include "complex.H"
#include "vector.H"
#include "tensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "Switch.H"

#include "foamNumpyCore.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

template<class T>
void printTypeName()
{
    // Both float and double have pTraits typeName = "scalar"!
    if constexpr (std::is_same_v<float, std::remove_cv_t<T>>)
    {
        Info<< "<float>";
    }
    else if constexpr (std::is_same_v<double, std::remove_cv_t<T>>)
    {
        Info<< "<double>";
    }
    else if constexpr (pTraits_has_typeName<T>)
    {
        Info<< pTraits<std::remove_cv_t<T>>::typeName
            << " ncmpt=" << int(pTraits_nComponents_v<T>);
    }
    else
    {
        Info<< error::demangle<T>();
    }
}


template<class T>
void printTraits()
{
    printTypeName<T>();

    Info<< ": numpy"
        << " char='" << numpy::dtypeChar<T>() << "'"
        << " size='" << numpy::dtypeSize<T>() << "'"
        << " descr=" << numpy::dtypeDescr<T>() << endl;
}


int main()
{
    printTraits<bool>();
    printTraits<label>();
    printTraits<scalar>();
    printTraits<int16_t>();
    printTraits<uint64_t>();
    printTraits<const complex>();     // Uses specialized pTraits_...
    printTraits<floatVector>();
    printTraits<doubleVector>();
    printTraits<sphericalTensor>();
    printTraits<symmTensor>();
    printTraits<tensor>();
    printTraits<boolVector>();  // Uses specialized pTraits_...
    printTraits<word>();
    printTraits<std::string>();

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
