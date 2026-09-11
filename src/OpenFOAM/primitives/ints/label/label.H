/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
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
    Foam::label

Description
    A label is an int32_t or int64_t as specified by the pre-processor macro
    WM_LABEL_SIZE.

    A readLabel function is defined so that label can be constructed from
    Istream.

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_label_H
#define Foam_primitives_label_H

#include "int.H"
#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef INT_ADD_SIZE
#undef INT_ADD_DEF_SIZE
#undef INT_SIZE

#define INT_ADD_SIZE(x,s,y) x ## s ## y
#define INT_ADD_DEF_SIZE(x,s,y) INT_ADD_SIZE(x,s,y)
#define INT_SIZE(x,y) INT_ADD_DEF_SIZE(x,WM_LABEL_SIZE,y)

// ------------------------------------------------------------------------- //

// Size checks and typedefs (label) in labelFwd.H

// Two-parameter forms of readLabel :
//   * readLabel(const char*, label&)
//   * readLabel(const std::string&, label&)
// are defined in int32.H and int64.H

// ------------------------------------------------------------------------- //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Parse entire buffer as a label, skipping leading/trailing whitespace.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
//  \return Parsed value or FatalIOError on any problem
inline label readLabel(const char* buf)
{
    return INT_SIZE(readInt,) (buf);
}

//- Parse entire string as a label, skipping leading/trailing whitespace.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
//  \return Parsed value or FatalIOError on any problem
inline label readLabel(const std::string& str)
{
    return INT_SIZE(readInt,) (str);
}

//- Read label from stream.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
inline label readLabel(Istream& is)
{
    return INT_SIZE(readInt,) (is);
}


//- Read raw label(s) from binary stream.
//  \note No check for binary vs ascii, the caller knows what they are doing
void readRawLabel(Istream& is, label* data, size_t nElem = 1);

//- Read raw label from binary stream.
//  \note No internal check for binary vs ascii,
//        the caller knows what they are doing
inline void readRawLabel(Istream& is, label& val)
{
    readRawLabel(is, &val, 1);
}

//- Read raw label from binary stream and return value.
//  \note No internal check for binary vs ascii,
//        the caller knows what they are doing
inline label readRawLabel(Istream& is)
{
    label val(0);
    readRawLabel(is, &val, 1);
    return val;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Raise one label to the power of another
label pow(label a, label b);

//- Evaluate n! : 0 < n <= 12
label factorial(label n);

//- Non-const access to integer-type (has no components)
inline label& setComponent(label& val, const direction) noexcept
{
    return val;
}

//- Return integer value (has no components)
inline constexpr label component(const label val, const direction) noexcept
{
    return val;
}


// * * * * * * * * * * * * * * * General Functions * * * * * * * * * * * * * //

//- Compare two values for equality
template<class T>
inline bool equal(const T& a, const T& b)
{
    return (a == b);
}

//- Return value clamped between upper and lower limits.
//  Unlike the std::clamp, which selects between references, this version
//  wraps the min/max free functions for component-wise clamping
template<class T>
inline T clamp(const T& val, const T& lower, const T& upper)
{
    return min(max(val, lower), upper);
}


/*---------------------------------------------------------------------------*\
                          Struct labelOp Declaration
\*---------------------------------------------------------------------------*/

//- Conversion/extraction to label operation
//  Specialization of this shall provide a corresponding \c operator().
template<class> struct labelOp;

//- Convert (likely promote) from int32_t to label
template<>
struct labelOp<int32_t>
{
    constexpr label operator()(const int32_t& val) const noexcept
    {
        return val;
    }
};


//- Convert (possibly truncate) from int64_t to label
template<>
struct labelOp<int64_t>
{
    constexpr label operator()(const int64_t& val) const noexcept
    {
        if constexpr (sizeof(label) == sizeof(int32_t))
        {
            return label(val);
        }
        else
        {
            return val;
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Type conversions (narrowing)

//- Type narrowing from int64_t to int32_t
//  Overflow: fix silently.
inline constexpr int32_t narrowInt32(const int64_t val) noexcept
{
    return
    (
        (val <= INT32_MIN) ? INT32_MIN
      : (val >= INT32_MAX) ? INT32_MAX
      : static_cast<int32_t>(val)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Label specific functions for which there is no uLabel equivalent

//- Return 1 if s is greater_equal zero, or otherwise -1
inline constexpr label sign(const label s) noexcept
{
    return (s >= 0)? 1: -1;
}


//- Return 1 if s is greater_equal zero, or otherwise 0
inline constexpr label pos0(const label s) noexcept
{
    return (s >= 0)? 1: 0;
}


//- Return 1 if s is less than zero, or otherwise 0
inline constexpr label neg(const label s) noexcept
{
    return (s < 0)? 1: 0;
}


//- Return 1 if s is less_equal zero, or otherwise 0
inline constexpr label neg0(const label s) noexcept
{
    return (s <= 0)? 1: 0;
}


//- Return the positive part of s, otherwise zero. Same as max(0, s).
inline constexpr label posPart(const label s) noexcept
{
    return (s > 0)? s: 0;
}


//- Return the negative part of s, otherwise zero. Same as min(0, s).
//- Does not change the sign.
inline constexpr label negPart(const label s) noexcept
{
    return (s < 0)? s: 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef INT_ADD_SIZE
#undef INT_ADD_DEF_SIZE
#undef INT_SIZE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
