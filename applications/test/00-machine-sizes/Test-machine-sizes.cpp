/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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
    Test the sizeof for basic types.
    Also tests how the data mapping of OpenFOAM types to UPstream (MPI)
    type ids are handled.

    Can be compiled and run without any OpenFOAM libraries.

        g++ -std=c++17 -oTest-machine-sizes Test-machine-sizes.cpp

\*---------------------------------------------------------------------------*/

#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Partial copy from UPstream.H

//- Some MPI data types
//
//- Mapping of some fundamental and aggregate types to MPI data types
enum class dataTypes : int
{
    // Fundamental Types [10]:
    Basic_begin,
    type_byte = Basic_begin,
    type_int16,
    type_int32,
    type_int64,
    type_uint16,
    type_uint32,
    type_uint64,
    type_float,
    type_double,
    type_long_double,

    // Tuple Types [4]:
    type_int_int,
    type_long_int,
    type_float_int,
    type_double_int,

    invalid,
    Basic_end = invalid
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Partial copy from UPstreamTraits.H

//- UPstream data type corresponding to an intrinsic (MPI) type
template<class T>
struct UPstream_mpi_dataType : std::false_type
{
    static constexpr auto datatype_id = dataTypes::invalid;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations to match elements of UPstream::dataTypes
#undef  defineUPstreamDataTraits
#define defineUPstreamDataTraits(TypeId, Type)                                \
    template<> struct UPstream_mpi_dataType<Type> : std::true_type            \
    {                                                                         \
        static constexpr auto datatype_id = dataTypes::TypeId;                \
    };


defineUPstreamDataTraits(type_byte,   char);
defineUPstreamDataTraits(type_byte,   unsigned char);
defineUPstreamDataTraits(type_int16,  int16_t);
defineUPstreamDataTraits(type_int32,  int32_t);
defineUPstreamDataTraits(type_int64,  int64_t);
defineUPstreamDataTraits(type_uint16, uint16_t);
defineUPstreamDataTraits(type_uint32, uint32_t);
defineUPstreamDataTraits(type_uint64, uint64_t);
defineUPstreamDataTraits(type_float,  float);
defineUPstreamDataTraits(type_double, double);
defineUPstreamDataTraits(type_long_double, long double);

#undef defineUPstreamDataTraits


//- Explicit handling of data type aliases. This is necessary since
//- different systems map things like 'unsigned long' differently but we
//- restrict ourselves to int32/int64 types
template<class T>
struct UPstream_alias_dataType
:
    std::bool_constant
    <
        // Basic MPI type
        UPstream_mpi_dataType<std::remove_cv_t<T>>::value ||
        (
            // Or some int 32/64 type to re-map
            std::is_integral_v<T>
         && (sizeof(T) == sizeof(int32_t) || sizeof(T) == sizeof(int64_t))
        )
    >
{
    using base = std::conditional_t
    <
        UPstream_mpi_dataType<std::remove_cv_t<T>>::value,
        std::remove_cv_t<T>,  // <- using mpi type (no alias)
        std::conditional_t    // <- using alias
        <
            (
                std::is_integral_v<T>
             && (sizeof(T) == sizeof(int32_t) || sizeof(T) == sizeof(int64_t))
            ),
            std::conditional_t
            <
                (sizeof(T) == sizeof(int32_t)),
                std::conditional_t<std::is_signed_v<T>, int32_t, uint32_t>,
                std::conditional_t<std::is_signed_v<T>, int64_t, uint64_t>
            >,
            char  // Fallback is a byte (eg, arbitrary contiguous data)
        >
    >;

    static constexpr auto datatype_id =
        UPstream_mpi_dataType<base>::datatype_id;
};


// Handle int8_t/uint8_t as aliases since 'signed char' etc may be
// ambiguous

//- Map \c int8_t to UPstream::dataTypes::type_byte
template<>
struct UPstream_alias_dataType<int8_t> : std::true_type
{
    using base = char;
    static constexpr auto datatype_id = dataTypes::type_byte;
};

//- Map \c uint8_t to UPstream::dataTypes::type_byte
template<>
struct UPstream_alias_dataType<uint8_t> : std::true_type
{
    using base = unsigned char;
    static constexpr auto datatype_id = dataTypes::type_byte;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void print(const char* name, bool showLimits = true)
{
    std::cout
        << "name=\"" << name << "\" sizeof=" << sizeof(T);

    if (showLimits)
    {
        std::cout
            << " max=<";

        if constexpr (sizeof(T) == 1)
        {
            std::cout << int(std::numeric_limits<T>::max());
        }
        else
        {
            std::cout << std::numeric_limits<T>::max();
        }
        std::cout << '>';
    }

    // A declared or deduced MPI type, or aliased
    if constexpr (UPstream_mpi_dataType<T>::value)
    {
        std::cout
            << " is_mpi=("
            << int(UPstream_mpi_dataType<T>::datatype_id) << ')';
    }
    else
    {
        std::cout << " is_mpi=(null)";
    }

    // Any aliases?
    if constexpr (UPstream_alias_dataType<T>::value)
    {
        if constexpr (UPstream_mpi_dataType<T>::value)
        {
            std::cout << " alias=base";
        }
        else
        {
            std::cout
                << " alias=("
                << int(UPstream_alias_dataType<T>::datatype_id) << ')';
        }
    }

    std::cout<< '\n';
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    std::cout<< "c++ = " << __cplusplus << '\n';
    std::cout<< "machine sizes (and some MPI traits)\n---\n\n";

    print<int8_t>("int8_t");
    print<uint8_t>("uint8_t");
    print<int16_t>("int16_t");
    print<uint16_t>("uint16_t");
    print<int32_t>("int32_t");
    print<uint32_t>("uint32_t");
    print<int64_t>("int64_t");
    print<uint64_t>("uint64_t");

    std::cout << '\n';
    print<char>("char");
    print<signed char>("signed char");
    print<unsigned char>("unsigned char");
    print<short>("short");
    print<int>("int");
    print<unsigned>("unsigned");
    print<long>("long");
    print<unsigned long>("unsigned long");
    print<long long>("long long");

    std::cout << '\n';
    print<std::size_t>("std::size_t");
    print<std::streamsize>("std::streamsize");

    std::cout << '\n';
    print<float>("float");
    print<double>("double");
    print<long double>("long double");

    std::cout << "\n---\nEnd\n\n";

    return 0;
}


// ************************************************************************* //
