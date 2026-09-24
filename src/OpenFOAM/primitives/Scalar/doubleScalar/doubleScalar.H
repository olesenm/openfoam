/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::doubleScalar

Description
    Floating-point double precision scalar type.

SourceFiles
    doubleScalar.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_double_H
#define Foam_primitives_double_H

#include "scalarFwd.H"
#include "doubleFloat.H"
#include "direction.H"
#include "pTraits.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Typedef (doubleScalar) and its limits in scalarFwd.H

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define Scalar double
#define ScalarVGREAT doubleScalarVGREAT
#define ScalarVSMALL doubleScalarVSMALL
#define ScalarROOTVGREAT doubleScalarROOTVGREAT
#define ScalarROOTVSMALL doubleScalarROOTVSMALL
#define ScalarRead readDouble

#if defined(FOAM_OFFLOAD) && defined(__HIP_DEVICE_COMPILE__)
__host__ __device__
#endif
inline Scalar mag(const Scalar s)
{
    return std::fabs(s);
}

#if defined(FOAM_OFFLOAD) && defined(__HIP_DEVICE_COMPILE__)
__host__ __device__
#endif
inline Scalar hypot(const Scalar x, const Scalar y)
{
    return std::hypot(x, y);
}

#if defined(FOAM_OFFLOAD) && defined(__HIP_DEVICE_COMPILE__)
__host__ __device__
#endif
inline Scalar atan2(const Scalar y, const Scalar x)
{
    return std::atan2(y, x);
}


// Normal (double-precision) transcendental functions
#if defined(FOAM_OFFLOAD) && defined(__HIP_DEVICE_COMPILE__)
// GPU version — use ROCm OCML intrinsics
#define transFunc(func)                                                        \
__host__ __device__ inline Scalar func(const Scalar s)                         \
{                                                                              \
    return __ocml_##func##_f64(static_cast<double>(s));                        \
}
#else
// CPU version — normal math library call
#define transFunc(func)                                                        \
inline Scalar func(const Scalar s)                                             \
{                                                                              \
    return ::func(s);                                                          \
}
#endif

// Normal (double-precision) bessel functions.
// May not be available on all systems
#ifdef Foam_no_besselFunc
    // Not available
    #define besselFunc(func)                                    \
    inline Scalar func(const Scalar s)                          \
    {                                                           \
        std::cerr<< "No '" << #func << "' function\n";          \
        return 0;                                               \
    }
    #define besselFunc2(func)                                   \
    inline Scalar func(const int n, const Scalar s)             \
    {                                                           \
        std::cerr<< "No '" << #func << "' function\n";          \
        return 0;                                               \
    }
#elif defined(__MINGW32__)
    // Mingw: with '_' prefix
    #define besselFunc(func)                                    \
    inline Scalar func(const Scalar s)                          \
    {                                                           \
        return _##func(s);                                      \
    }
    #define besselFunc2(func)                                   \
    inline Scalar func(const int n, const Scalar s)             \
    {                                                           \
        return _##func(n, s);                                   \
    }
#else
    #define besselFunc(func)                                    \
    inline Scalar func(const Scalar s)                          \
    {                                                           \
        return ::func(s);                                       \
    }
    #define besselFunc2(func)                                   \
    inline Scalar func(const int n, const Scalar s)             \
    {                                                           \
        return ::func(n, s);                                    \
    }
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// scalarImpl.H is used for template-like substitution (but using macros)

#define Foam_use_scalarImpl_header
#include "scalarImpl.H"
#undef Foam_use_scalarImpl_header

#undef Scalar
#undef ScalarVGREAT
#undef ScalarVSMALL
#undef ScalarROOTVGREAT
#undef ScalarROOTVSMALL
#undef ScalarRead
#undef transFunc
#undef besselFunc
#undef besselFunc2


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
