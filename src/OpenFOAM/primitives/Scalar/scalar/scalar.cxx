/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "scalar.H"
#include "stdFoam.H"
#include "IOstreams.H"

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
                // Narrowing: currently only need (float <- double)
                *data = Foam::narrowFloat(val);
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
        (sizeof(dataType) == sizeof(float))
     || (sizeof(dataType) == sizeof(double))
    )
    {
        switch (is.scalarByteSize())
        {
            case sizeof(float) : reading<float>(is, data, nElem); break;
            case sizeof(double): reading<double>(is, data, nElem); break;

            default:
            {
                // Cannot recover from this
                using namespace Foam;

                FatalIOErrorInFunction(is)
                    << "Currently no code to read float"
                    << (8*is.scalarByteSize())
                    << " as float" << (8*sizeof(dataType)) << nl
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
            "Only float or double are currently supported"
        );
    }
}


// Read scalar from Istream (if good and correct token type)
template<class dataType>
bool try_read(Foam::Istream& is, dataType& value)
{
    if (is.good())
    {
        Foam::token tok(is);

        // NB: does not handle separated '-' (or '+') prefixes
        // like operator>>(Istream&, scalar&) does

        if (tok.isNumber())
        {
            value = tok.number();
            return true;
        }

        is.putBack(tok);
    }

    return false;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

#if defined(WM_SP) || defined(WM_SPDP) || defined(WM_DP)

void Foam::readRawScalar(Istream& is, scalar* data, size_t nElem)
{
    dispatch_reading(is, data, nElem);
}


Foam::scalar Foam::readScalarOrDefault(Istream& is, scalar defaultValue)
{
    if (scalar val(0); try_read(is, val))
    {
        return val;
    }
    else
    {
        return defaultValue;
    }
}

#endif

// ************************************************************************* //
