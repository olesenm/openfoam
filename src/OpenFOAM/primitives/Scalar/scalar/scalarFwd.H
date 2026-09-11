/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Typedefs for float/double/scalar and their numerical limits
    without requiring scalar.H or <limits>

    The definition of scalar as a floating-point number depends
    on pre-processor macros WM_SP, WM_SPDP or WM_DP being defined.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Foam_scalarFwd_H
#define Foam_scalarFwd_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A typedef for float
typedef float floatScalar;

//- A typedef for double
typedef double doubleScalar;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Largest and smallest scalar values allowed in certain parts of the code.
// See std::numeric_limits max(), min(), epsilon()

// limits <float>
constexpr floatScalar floatScalarGREAT = 1.0e+6;
constexpr floatScalar floatScalarROOTGREAT = 1.0e+3;
constexpr floatScalar floatScalarVGREAT = 1.0e+37;
constexpr floatScalar floatScalarROOTVGREAT = 1.0e+18;
constexpr floatScalar floatScalarSMALL = 1.0e-6;
constexpr floatScalar floatScalarROOTSMALL = 1.0e-3;
constexpr floatScalar floatScalarVSMALL = 1.0e-37;
constexpr floatScalar floatScalarROOTVSMALL = 1.0e-18;

// limits <double>
constexpr doubleScalar doubleScalarGREAT = 1.0e+15;
constexpr doubleScalar doubleScalarROOTGREAT = 3.0e+8;
constexpr doubleScalar doubleScalarVGREAT = 1.0e+300;
constexpr doubleScalar doubleScalarROOTVGREAT = 1.0e+150;
constexpr doubleScalar doubleScalarSMALL = 1.0e-15;
constexpr doubleScalar doubleScalarROOTSMALL = 3.0e-8;
constexpr doubleScalar doubleScalarVSMALL = 1.0e-300;
constexpr doubleScalar doubleScalarROOTVSMALL = 1.0e-150;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// OpenFOAM scalar and solveScalar typedefs, depending on compilation precision

#if defined(WM_SP)

typedef floatScalar scalar;
typedef floatScalar solveScalar;

#elif defined(WM_SPDP)

typedef floatScalar scalar;
typedef doubleScalar solveScalar;

#elif defined(WM_DP)

typedef doubleScalar scalar;
typedef doubleScalar solveScalar;

#else
//    #error "PRECISION must be set to WM_SP, WM_SPDP or WM_DP"
// Provide fallback?
// typedef doubleScalar scalar;
// typedef doubleScalar solveScalar;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(WM_SP) || defined(WM_SPDP)

// With scalar == (float), solveScalar == (float | double)

constexpr auto GREAT = floatScalarGREAT;
constexpr auto ROOTGREAT = floatScalarROOTGREAT;
constexpr auto VGREAT = floatScalarVGREAT;
constexpr auto ROOTVGREAT = floatScalarROOTVGREAT;
constexpr auto SMALL = floatScalarSMALL;
constexpr auto ROOTSMALL = floatScalarROOTSMALL;
constexpr auto VSMALL = floatScalarVSMALL;
constexpr auto ROOTVSMALL = floatScalarROOTVSMALL;

#elif defined(WM_DP)

// With scalar == (double), solveScalar == (double)

constexpr auto GREAT = doubleScalarGREAT;
constexpr auto ROOTGREAT = doubleScalarROOTGREAT;
constexpr auto VGREAT = doubleScalarVGREAT;
constexpr auto ROOTVGREAT = doubleScalarROOTVGREAT;
constexpr auto SMALL = doubleScalarSMALL;
constexpr auto ROOTSMALL = doubleScalarROOTSMALL;
constexpr auto VSMALL = doubleScalarVSMALL;
constexpr auto ROOTVSMALL = doubleScalarROOTVSMALL;

#endif

#if defined(COMPAT_OPENFOAM_ORG) && \
    (defined(WM_SP) || defined(WM_SPDP) || defined(WM_DP))
//! \cond compat_openfoam_org
//  Accommodate name changes from 2018-01
constexpr auto great = GREAT;
constexpr auto rootGreat = ROOTGREAT;
constexpr auto vGreat = VGREAT;
constexpr auto rootVGreat = ROOTVGREAT;
constexpr auto small = SMALL;
constexpr auto rootSmall = ROOTSMALL;
constexpr auto vSmall = VSMALL;
constexpr auto rootVSmall = ROOTVSMALL;
//! \endcond
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
