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

Namespace
    Foam::numpy

Description
    Namespace for handling NUMPY related input/output.
    https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html

Class
    Foam::fileFormats::numpyCore

Description
    Core routines for writing NUMPY files.

SourceFiles
    foamNumpyCore.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_numpyCore_H
#define Foam_numpyCore_H

#include "pTraits.H"
#include "foamEndianFwd.H"  // For big/little endian defines
#include "UList.H"

#include <cstdint>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace numpy
{

// Data Types

//- Enumeration of basic numpy data types, encoded with the type identifier
//- and nbytes. No complex data types for now
enum class dataType : std::uint16_t
{
    type_bool  = (('b' << 8) | 1),
    type_int8  = (('i' << 8) | 1),
    type_int16 = (('i' << 8) | 2),
    type_int32 = (('i' << 8) | 4),
    type_int64 = (('i' << 8) | 8),
    type_uint8  = (('u' << 8) | 1),
    type_uint16 = (('u' << 8) | 2),
    type_uint32 = (('u' << 8) | 4),
    type_uint64 = (('u' << 8) | 8),
    type_float8  = (('f' << 8) | 1),
    type_float16 = (('f' << 8) | 2),
    type_float32 = (('f' << 8) | 4),
    type_float64 = (('f' << 8) | 8),
    invalid = 0
};

//- The numpy endian char for this system
#ifdef WM_BIG_ENDIAN
static constexpr char endian_char = '>';
#else
static constexpr char endian_char = '<';
#endif


//- The numpy char name for given numpy::dataType
inline constexpr char dtypeChar(numpy::dataType dtype) noexcept
{
    return ((std::uint16_t(dtype) >> 8) & 0xFF);   // eg, 'f', 'i' etc
}

//- The size (bytes) for given numpy::dataType
inline constexpr int dtypeSize(numpy::dataType dtype) noexcept
{
    return (std::uint16_t(dtype) & 0xFF);
}

//- The numpy char name for given scalar/vector-space component type
template<class T>
constexpr char dtypeChar() noexcept
{
    using cmpt = typename pTraits_cmptType<std::remove_cv_t<T>>::type;

    /**/ if constexpr (std::is_same_v<bool,cmpt>) { return 'b'; }
    else if constexpr (std::is_unsigned_v<cmpt>)  { return 'u'; }
    else if constexpr (std::is_integral_v<cmpt>)  { return 'i'; }
    else if constexpr (std::is_floating_point_v<cmpt>) { return 'f'; }
    else
    {
        // TBD:
        // static_assert(stdFoam::dependent_false_v<T>, "Unmapped data type");
        return '\0';  // Fallback as 'invalid'
    }
}

//- The size (bytes) for given scalar/vector-space component type
template<class T>
constexpr int dtypeSize() noexcept
{
    using cmpt = typename pTraits_cmptType<std::remove_cv_t<T>>::type;

    /**/ if constexpr (std::is_same_v<bool, cmpt>) { return 1; }
    else if constexpr (std::is_arithmetic_v<cmpt>) { return sizeof(cmpt); }
    else
    {
        // TBD:
        // static_assert(stdFoam::dependent_false_v<T>, "Unmapped data type");
        return 0;  // Fallback as 'invalid'
    }
}

//- The numpy string representation for given scalar/vector-space component type
template<class T>
std::string dtypeDescr()
{
    using cmpt = typename pTraits_cmptType<std::remove_cv_t<T>>::type;

    if constexpr (std::is_same_v<bool, cmpt>)
    {
        return "|b1";   // No endian for bool (1-byte)
    }
    else if constexpr (std::is_arithmetic_v<cmpt>)
    {
        constexpr int nbytes = numpy::dtypeSize<T>();

        if constexpr (nbytes < 10)
        {
            // Normal case (eg, "<f8")
            std::string descr(3, numpy::endian_char);
            descr[1] = numpy::dtypeChar<T>();   // ('f', 'i', ...)
            descr[2] = char('0' + nbytes);
            return descr;
        }
        else
        {
            // Future-proof case (eg, "<f16")
            std::string descr(2, numpy::endian_char);
            descr[1] = numpy::dtypeChar<T>();   // ('f', 'i', ...)
            descr += std::to_string(nbytes);
            return descr;
        }
    }
    else
    {
        return std::string();  // or return "|O" (object) ?
    }
}

//- The string representation for given numpy::dataType.
inline std::string dtypeDescr(numpy::dataType dtype)
{
    if (dtype == numpy::dataType::invalid)
    {
        return std::string();  // or return "|O" (object) ?
    }
    else if (dtype == numpy::dataType::type_bool)
    {
        return "|b1";   // No endian for bool (1-byte)
    }
    else
    {
        // NB: the current enumeration always fits within 3 chars
        std::string descr(3, numpy::endian_char);
        descr[1] = numpy::dtypeChar(dtype);   // ('f', 'i', ...)
        descr[2] = char('0') + numpy::dtypeSize(dtype);  // nbytes
        return descr;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numpy
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::numpyCore Declaration
\*---------------------------------------------------------------------------*/

class numpyCore
{
public:

    // Helpers (low-level)

        //- Numpy shape description for data
        template<class T>
        static std::string dataShape(label count)
        {
            constexpr auto nCmpts = int(pTraits_nComponents_v<T>);
            if constexpr (is_vectorspace_v<T> || nCmpts > 1)
            {
                // vector-space or complex as (count,ncmpt)
                return
                (
                    "(" + std::to_string(count)
                  + "," + std::to_string(nCmpts) + ")"
                );
            }
            else
            {
                // scalar types as (count,)
                return "(" + std::to_string(count) + ",)";
            }
        }

        //- Numpy shape description for data content
        template<class T>
        static std::string dataShape(const UList<T>& content)
        {
            return dataShape<T>(content.size());
        }

        //- Create a numpy header dictionary.
        //- Includes padding and trailing newline
        static std::string headerDict
        (
            //! The numpy data type (eg, \c "f4")
            const std::string& dtype,
            //! The numpy data shape
            const std::string& shape,
            //! Output in fortran order
            bool fortranOrder = false
        );

        //- Write numpy header to stream
        static std::size_t writeHeader
        (
            std::ostream& os,
            //! The numpy data type (eg, \c "<f8")
            const std::string& dtype,
            //! The numpy data shape
            const std::string& shape,
            //! Output in fortran order
            bool fortranOrder = false
        );


    // Helpers (high-level)

        //- Write numpy header to stream for given OpenFOAM data
        template<class T>
        static std::size_t writeHeader
        (
            std::ostream& os,
            //! Data content (for the shape description)
            const UList<T>& content,
            //! Output in fortran order
            bool fortranOrder = false
        )
        {
            return writeHeader
            (
                os,
                numpy::dtypeDescr<T>(),
                dataShape(content),
                fortranOrder
            );
        }

        //- From OpenFOAM data type (scalar, vector, etc)
        //- to numpy shape description
        template<class T>
        static void writeData(std::ostream& os, const UList<T>& content)
        {
            os.write(content.cdata_bytes(), content.size_bytes());
        }


protected:

    // Constructors

        //- Default construct
        numpyCore() noexcept = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
