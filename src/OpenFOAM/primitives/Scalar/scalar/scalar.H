/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Typedef
    Foam::scalar

Description
    A floating-point number identical to float or double depending on
    whether WM_SP, WM_SPDP or WM_DP is defined.

SourceFiles
    scalar.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_scalar_H
#define Foam_scalar_H

#include "floatScalar.H"
#include "doubleScalar.H"

// ------------------------------------------------------------------------- //

// Typedefs (floatScalar, doubleScalar, scalar, solveScalar)
// as well as their limits are defined in scalarFwd.H

// Two-parameter forms of readScalar :
//   * readScalar(const char*, scalar&)
//   * readScalar(const std::string&, scalar&)
// are defined in floatScalar.H and doubleScalar.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(WM_SP) || defined(WM_SPDP)

// With scalar == (float), solveScalar == (float | double)

namespace Foam
{

    //- Read scalar from c-string and return value
    inline scalar readScalar(const char* buf)
    {
        return readFloat(buf);
    }

    //- Read scalar from string and return value
    inline scalar readScalar(const std::string& str)
    {
        return readFloat(str);
    }

    //- Read scalar from stream and return value
    inline scalar readScalar(Istream& is)
    {
        return readFloat(is);
    }

    //- Read scalar from stream if present or return default value
    scalar readScalarOrDefault(Istream& is, scalar defaultValue);

    //- Read raw scalar(s) from binary stream.
    //  \note No check for binary vs ascii, the caller knows what they are doing
    void readRawScalar(Istream& is, scalar* data, size_t nElem = 1);

    //- Read raw scalar from binary stream.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    inline void readRawScalar(Istream& is, scalar& val)
    {
        readRawScalar(is, &val, 1);
    }

    //- Read raw scalar from binary stream and return value.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    inline scalar readRawScalar(Istream& is)
    {
        scalar val(0);
        readRawScalar(is, &val, 1);
        return val;
    }

} // End namespace Foam

#elif defined(WM_DP)

// With scalar == (double), solveScalar == (double)

namespace Foam
{

    //- Read scalar from c-string and return value
    inline scalar readScalar(const char* buf)
    {
        return readDouble(buf);
    }

    //- Read scalar from string and return value
    inline scalar readScalar(const std::string& str)
    {
        return readDouble(str);
    }

    //- Read scalar from stream and return value
    inline scalar readScalar(Istream& is)
    {
        return readDouble(is);
    }

    //- Read scalar from stream if present or return default value
    scalar readScalarOrDefault(Istream& is, scalar defaultValue);

    //- Read raw scalar(s) from binary stream.
    //  \note No check for binary vs ascii, the caller knows what they are doing
    void readRawScalar(Istream& is, scalar* data, size_t nElem = 1);

    //- Read raw scalar from binary stream.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    inline void readRawScalar(Istream& is, scalar& val)
    {
        readRawScalar(is, &val, 1);
    }

    //- Read raw scalar from binary stream and return value.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    inline scalar readRawScalar(Istream& is)
    {
        scalar val(0);
        readRawScalar(is, &val, 1);
        return val;
    }

} // End namespace Foam

#else
//    #error "PRECISION must be set to WM_SP, WM_SPDP or WM_DP"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Type conversions (narrowing)

namespace Foam
{

//- Type narrowing from double to float
//  Underflow: round to zero.
//  Overflow: fix silently.
inline constexpr float narrowFloat(const double val) noexcept
{
    return
    (
        (val <= -floatScalarVGREAT) ? -floatScalarVGREAT
      : (val >=  floatScalarVGREAT) ?  floatScalarVGREAT
      : (val > -floatScalarVSMALL && val < floatScalarVSMALL) ? 0
      : static_cast<float>(val)
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Type to use for extended precision
    template<>
    class typeOfSolve<scalar>
    {
    public:

        typedef solveScalar type;
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
