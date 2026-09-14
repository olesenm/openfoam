/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "error.H"
#include "label.H"
#include "stdFoam.H"
#include "Istream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if WM_LABEL_SIZE == 32
const char* const Foam::pTraits<int32_t>::typeName = "label";
const char* const Foam::pTraits<int64_t>::typeName = "int64";
#elif WM_LABEL_SIZE == 64
const char* const Foam::pTraits<int32_t>::typeName = "int32";
const char* const Foam::pTraits<int64_t>::typeName = "label";
#endif


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Binary reading with narrowing/widening
template<class readType, class dataType>
void reading(Foam::Istream& is, dataType* data, size_t nElem)
{
    if constexpr (sizeof(dataType) == sizeof(readType))
    {
        // Read uses the native data size
        is.readRaw(reinterpret_cast<char*>(data), nElem*sizeof(dataType));
    }
    else
    {
        if (data == nullptr)
        {
            // Even handle 'seek' behaviour
            is.readRaw(nullptr, nElem*sizeof(readType));
            return;
        }

        for (const dataType* last = data + nElem; data != last; ++data)
        {
            readType val;
            is.readRaw(reinterpret_cast<char*>(&val), sizeof(readType));

            if constexpr (sizeof(dataType) < sizeof(readType))
            {
                // Narrowing: currently only need (int32 <- int64)
                *data = Foam::narrowInt32(val);
            }
            else
            {
                // Type widening
                *data = dataType(val);
            }
        }
    }
}


// Dispatch for binary reading with narrowing/widening
// No check for binary vs ascii, the caller knows what they are doing
template<class dataType>
void dispatch_reading(Foam::Istream& is, dataType* data, size_t nElem)
{
    if constexpr
    (
        (sizeof(dataType) == sizeof(int32_t))
     || (sizeof(dataType) == sizeof(int64_t))
    )
    {
        switch (is.labelByteSize())
        {
            case sizeof(int32_t) : reading<int32_t>(is, data, nElem); break;
            case sizeof(int64_t) : reading<int64_t>(is, data, nElem); break;

            default:
            {
                // Cannot recover from this
                using namespace Foam;

                FatalIOErrorInFunction(is)
                    << "Currently no code to read int"
                    << (8*is.labelByteSize())
                    << " as int" << (8*sizeof(dataType)) << nl
                    << abort(FatalIOError);
                break;
            }
        }
    }
    else
    {
        static_assert
        (
            stdFoam::dependent_false_v<dataType>,
            "Only int32 or int64 are currently supported"
        );
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::readRawLabel(Foam::Istream& is, Foam::label* data, size_t nElem)
{
    dispatch_reading(is, data, nElem);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::pow(label a, label b)
{
    label ans = 1;
    for (label i=0; i<b; i++)
    {
        ans *= a;
    }

    #ifdef FULLDEBUG
    if (b < 0)
    {
        FatalErrorInFunction
            << "negative value for b is not supported"
            << abort(FatalError);
    }
    #endif

    return ans;
}


Foam::label Foam::factorial(label n)
{
    static label factTable[13] =
    {
        1, 1, 2, 6, 24, 120, 720, 5040, 40320,
        362880, 3628800, 39916800, 479001600
    };

    #ifdef FULLDEBUG
    if (n < 0 || n > 12)
    {
        FatalErrorInFunction
            << "n value out of range"
            << abort(FatalError);
    }
    #endif

    return factTable[n];
}


// ************************************************************************* //
