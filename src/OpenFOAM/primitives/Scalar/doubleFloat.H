/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef Foam_doubleFloat_H
#define Foam_doubleFloat_H

#include "label.H"
#include "products.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#if defined(FOAM_OFFLOAD) && defined(__HIP_DEVICE_COMPILE__)
    // Device compilation (ROCm)
    #define FOAM_POW_CALL(base, exp)                                           \
        __ocml_pow_f64(static_cast<double>(base), static_cast<double>(exp))
#else
    // Host (CPU) path
    #define FOAM_POW_CALL(base, exp)                                           \
    (std::pow(double(base), double(exp)))
#endif

#if defined(__HIP_DEVICE_COMPILE__)
    #define FOAM_HOST_DEVICE __host__ __device__
#else
    #define FOAM_HOST_DEVICE
#endif

#undef  MAXMINPOW
#define MAXMINPOW(RetType, Type1, Type2)                                       \
                                                                               \
/*! \brief Floating/integral min. Use std::min() to preserve references */     \
inline constexpr RetType min(Type1 s1, Type2 s2) noexcept                      \
{                                                                              \
    return (s1 < s2) ? s1 : s2;                                                \
}                                                                              \
                                                                               \
/*! \brief Floating/integral max. Use std::max() to preserve references */     \
inline constexpr RetType max(Type1 s1, Type2 s2) noexcept                      \
{                                                                              \
    return (s2 < s1) ? s1 : s2;                                                \
}                                                                              \
                                                                               \
/*! \brief Raise base to the power expon, returning double */                  \
FOAM_HOST_DEVICE                                                               \
inline RetType pow(const Type1 base, const Type2 expon)                        \
{                                                                              \
    return FOAM_POW_CALL(base, expon);                                         \
}

MAXMINPOW(double, double, double)
MAXMINPOW(double, double, float)
MAXMINPOW(double, float, double)
MAXMINPOW(double, double, int)
MAXMINPOW(double, int, double)
MAXMINPOW(double, double, long)
MAXMINPOW(double, long, double)
MAXMINPOW(float, float, float)
MAXMINPOW(float, float, int)
MAXMINPOW(float, int, float)
MAXMINPOW(float, float, long)
MAXMINPOW(float, long, float)
#if (WM_LABEL_SIZE == 64) && (defined(__APPLE__) || defined(_WIN32))
MAXMINPOW(double, double, int64_t)
MAXMINPOW(double, int64_t, double)
MAXMINPOW(float, float, int64_t)
MAXMINPOW(float, int64_t, float)
#endif

#undef MAXMINPOW
#undef FOAM_HOST_DEVICE


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
