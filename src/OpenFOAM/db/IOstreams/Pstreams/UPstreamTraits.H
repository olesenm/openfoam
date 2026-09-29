/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    A set of traits associated with UPstream communication

    \par Front-facing traits

    - \c UPstream_basic_dataType trait :<br>
      The main entry point for reduction operations
      (requires fundamental types).<br>
      Checks against fundamental/aliased (excludes user-defined),
      or is a component-wise aggregate of the same.

    - \c UPstream_dataType trait :<br>
      The main entry point for transmission (broadcast, send/recv, ...).<br>
      Checks against fundamental/aliased/user-defined,
      or is a component-wise aggregate of the same.

    - \c UPstream_opType trait :<br>
      Mapping of OpenFOAM ops to their MPI equivalent.
      The \c opcode_id is the corresponding internal representation.

    - \c UPstream_data_opType trait :<br>
      Combination of UPstream_opType and UPstream_basic_dataType.
    .

    \par Additional helper traits (not normally used directly):

    - \c UPstream_mpi_dataType trait :<br>
      Tests true and provides valid \c datatype_id for MPI fundamental
      data types. This trait will should not normally be used directly:
      use UPstream_alias_dataType for 'low-level' purposes (see below).

    - \c UPstream_user_dataType trait :<br>
      Tests true and provides valid \c datatype_id for user-defined
      data types.

    - \c UPstream_alias_dataType trait :<br>
      Use this in preference to UPstream_mpi_dataType.<br>
      A pass-through to UPstream_mpi_dataType, but provides additional
      mappings for <int/long/long long,...> to the fundamental
      32/64 bit integrals, since <int/long/long long,...> may not otherwise
      map directly on all systems.

    - \c UPstream_any_dataType trait :<br>
      Used as a building block for uniform aggregate types.<br>
      Combines UPstream_user_dataType and UPstream_alias_dataType into a
      single trait.
    .

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstreamTraits_H
#define Foam_UPstreamTraits_H

#include "UPstream.H"
#include <cstdint>
#include <ios>  // For streamsize
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

// Some vector-space types
// -----------------------
//! \cond
template<class T> class Vector;
template<class T> class SymmTensor;
template<class T> class Tensor;
//! \endcond

// -------------------------
// Some binary operators (as per ops.H), but since ListOps.H is included
// by UPstream.H, don't need to forward declare
// -------------------------
// template<class T> struct minOp;
// template<class T> struct maxOp;
// template<class T> struct plusOp;
// template<class T> struct sumOp;
// template<class T> struct multiplyOp;
// template<class T> struct bitAndOp;
// template<class T> struct bitOrOp;
// template<class T> struct bitXorOp;

//! \cond
template<class T> struct UPstream_basic_dataType;
template<class T> struct UPstream_dataType;
//! \endcond


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- UPstream data type corresponding to a fundamental (MPI) type
template<class T>
struct UPstream_mpi_dataType : std::false_type
{
    static constexpr auto datatype_id = UPstream::dataTypes::invalid;
};

//- Disallow \c void
template<> struct UPstream_mpi_dataType<void> : std::false_type
{
    static constexpr auto datatype_id = UPstream::dataTypes::invalid;
};

//- UPstream data type corresponding to user-defined type
template<class T>
struct UPstream_user_dataType : std::false_type
{
    static constexpr auto datatype_id = UPstream::dataTypes::invalid;
};

//- Disallow \c void
template<> struct UPstream_user_dataType<void> : std::false_type
{
    static constexpr auto datatype_id = UPstream::dataTypes::invalid;
};

//- UPstream data type for a fundamental (MPI) location tuple,
//- based on the first tuple element (second is always \c int).
//- Limited to a very few pre-defined types!
template<class T>
struct UPstream_location_dataType : std::false_type
{
    //- The basic (C) type for the tuple first element
    using type = T;
    static constexpr auto datatype_id = UPstream::dataTypes::invalid;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Trait specializations (data types)

// Specializations to match elements of UPstream::dataTypes.

#undef  defineUPstreamDataTraits
#define defineUPstreamDataTraits(TypeId, Type)                                \
                                                                              \
    /*! \brief Map \c Type to UPstream::dataTypes::TypeId */                  \
    template<> struct UPstream_mpi_dataType<Type> : std::true_type            \
    {                                                                         \
        static constexpr auto datatype_id = UPstream::dataTypes::TypeId;      \
    };                                                                        \
    /*! \brief Map \c const \c Type to \c UPstream::dataTypes::TypeId */      \
    template<> struct UPstream_mpi_dataType<const Type> : std::true_type      \
    {                                                                         \
        static constexpr auto datatype_id = UPstream::dataTypes::TypeId;      \
    };

// Fundamental Types [10]:
// Note: uses 'int32_t,int64_t,...' instead of 'int,long,...' to minimize
// the possibility of duplicate types. However, 'int8_t,uint8_t' are treated
// as aliases (char,unsigned char) to avoid possible compilation issues.
//
// OpenFOAM defines Foam::label as either int32_t,int64_t (not int,long) too.
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

// ------------------------------------------------------------------------- //

#undef  defineUPstreamDataTraits
#define defineUPstreamDataTraits(TypeId, Type)                                \
                                                                              \
    /*! \brief Map \c Type to UPstream::dataTypes::TypeId */                  \
    template<> struct UPstream_user_dataType<Type> : std::true_type           \
    {                                                                         \
        static constexpr auto datatype_id = UPstream::dataTypes::TypeId;      \
    };                                                                        \
    /*! \brief Map \c const \c Type to \c UPstream::dataTypes::TypeId */      \
    template<> struct UPstream_user_dataType<const Type> : std::true_type     \
    {                                                                         \
        static constexpr auto datatype_id = UPstream::dataTypes::TypeId;      \
    };


// User Types [6]:
defineUPstreamDataTraits(type_3float,  Vector<float>);
defineUPstreamDataTraits(type_3double, Vector<double>);
defineUPstreamDataTraits(type_6float,  SymmTensor<float>);
defineUPstreamDataTraits(type_6double, SymmTensor<double>);
defineUPstreamDataTraits(type_9float,  Tensor<float>);
defineUPstreamDataTraits(type_9double, Tensor<double>);

#undef defineUPstreamDataTraits

// ------------------------------------------------------------------------- //

#undef  defineUPstreamDataTraits
#define defineUPstreamDataTraits(TypeId, Type, BasicType)                     \
                                                                              \
    /*! \brief Map tuple \c (Type, int) to UPstream::dataTypes::TypeId */     \
    template<> struct UPstream_location_dataType<Type> : std::true_type       \
    {                                                                         \
        using type = BasicType;                                               \
        static constexpr auto datatype_id = UPstream::dataTypes::TypeId;      \
    };

// Location tuple Types [4]:
defineUPstreamDataTraits(type_int_int,  int32_t, int);
defineUPstreamDataTraits(type_long_int, int64_t, long);
defineUPstreamDataTraits(type_float_int, float, float);
defineUPstreamDataTraits(type_double_int, double, double);

#undef defineUPstreamDataTraits

// ------------------------------------------------------------------------- //

//- Explicit handling of data type aliases. This is necessary since
//- different systems map things like 'unsigned long' differently but we
//- restrict ourselves to int32/int64 types.
//
//  Note that this trait serves as the single pass-through point when needing
//  to reference UPstream_mpi_dataType elsewhere
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
    static constexpr auto datatype_id = UPstream::dataTypes::type_byte;
};

//- Map \c uint8_t to UPstream::dataTypes::type_byte
template<>
struct UPstream_alias_dataType<uint8_t> : std::true_type
{
    using base = unsigned char;
    static constexpr auto datatype_id = UPstream::dataTypes::type_byte;
};

// ------------------------------------------------------------------------- //

//- UPstream data type (fundamental or user-defined),
//- after resolving any aliases
template<class T>
struct UPstream_any_dataType
:
    std::bool_constant
    <
        UPstream_user_dataType<std::remove_cv_t<T>>::value
     || UPstream_alias_dataType<T>::value
    >
{
    using base = std::conditional_t
    <
        UPstream_user_dataType<std::remove_cv_t<T>>::value,
        std::remove_cv_t<T>,
        typename UPstream_alias_dataType<T>::base
    >;

    //- The corresponding UPstream::dataTypes enumeration
    static constexpr auto datatype_id = []() constexpr noexcept
    {
        if constexpr (UPstream_user_dataType<std::remove_cv_t<T>>::value)
        {
            // A user-defined type
            return UPstream_user_dataType<std::remove_cv_t<T>>::datatype_id;
        }
        else if constexpr (UPstream_alias_dataType<T>::value)
        {
            // Fundamental type or alias to a fundamental type
            return UPstream_alias_dataType<T>::datatype_id;
        }
        else
        {
            return UPstream::dataTypes::invalid;
        }
    }();
};


// ------------------------------------------------------------------------- //

//- UPstream fundamental/aliased (excludes user-defined) data type
//- or a component aggregate of the same.
//
// True for the following conditions:
// - The \c Type is directly supported
// - The \c cmptType (eg, from VectorSpace) exists and is directly supported
// - Fallback to byte-wise representation (ie, for contiguous)
// .
template<class T>
struct UPstream_basic_dataType
:
    std::bool_constant
    <
        // non-aggregate type
        UPstream_alias_dataType<T>::value
        // aggregate type
     || UPstream_alias_dataType<typename pTraits_cmptType<T>::type>::value
    >
{
    //- The underlying data type (if supported) or byte
    using base = std::conditional_t
    <
        UPstream_alias_dataType<T>::value,
        typename UPstream_alias_dataType<T>::base,  // <- non-aggregate
        typename UPstream_alias_dataType
        <typename pTraits_cmptType<T>::type>::base  // <- aggregate
    >;

    //- The corresponding UPstream::dataTypes enumeration
    static constexpr auto datatype_id =
        UPstream_alias_dataType<base>::datatype_id;

    //- The size in terms of the number of underlying data elements
    static constexpr std::streamsize size(std::streamsize n) noexcept
    {
        if constexpr (UPstream_alias_dataType<T>::value)
        {
            // non-aggregate: no multiplier
            return n;
        }
        else
        {
            // aggregate: with multiplier
            return n*(sizeof(T)/sizeof(base));
        }
    }
};

//- Disallow \c void
template<> struct UPstream_basic_dataType<void> : UPstream_mpi_dataType<void>
{
    using base = void;
    static constexpr std::streamsize size(std::streamsize n) noexcept
    {
        return n;
    }
};


// ------------------------------------------------------------------------- //

//- UPstream fundamental/aliased/user-defined data type
//- or a component aggregate of the same.
//
// True for the following conditions:
// - The \c Type is directly supported
// - The \c cmptType (eg, from VectorSpace) exists and is directly supported
// - Fallback to byte-wise representation (ie, for contiguous)
// .
template<class T>
struct UPstream_dataType
:
    std::bool_constant
    <
        UPstream_any_dataType<T>::value
     || UPstream_any_dataType<typename pTraits_cmptType<T>::type>::value
    >
{
    //- The underlying data type (if supported) or byte
    using base = std::conditional_t
    <
        UPstream_any_dataType<T>::value,
        typename UPstream_any_dataType<T>::base,   // <- non-aggregate
        typename UPstream_any_dataType
        <typename pTraits_cmptType<T>::type>::base  // <- aggregate
    >;

    //- The corresponding UPstream::dataTypes enumeration
    static constexpr auto datatype_id =
        UPstream_any_dataType<base>::datatype_id;

    //- The size in terms of the number of base data elements
    static constexpr std::streamsize size(std::streamsize n) noexcept
    {
        if constexpr (UPstream_any_dataType<T>::value)
        {
            // non-aggregate: no multiplier
            return n;
        }
        else
        {
            // aggregate: with multiplier
            return n*(sizeof(T)/sizeof(base));
        }
    }
};

//- Disallow \c void
template<> struct UPstream_dataType<void> : UPstream_basic_dataType<void> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Reduction op-codes

//- A supported UPstream (MPI) reduce/window operation type
template<class BinaryOp>
struct UPstream_opType : std::false_type
{
    static constexpr auto opcode_id = UPstream::opCodes::invalid;
};

//- Disallow \c void
template<> struct UPstream_opType<void> : std::false_type
{
    static constexpr auto opcode_id = UPstream::opCodes::invalid;
};

//- Map minOp\<T\> to \c UPstream::opCodes::op_min
template<class T>
struct UPstream_opType<Foam::minOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_min;
};

//- Map maxOp\<T\> to \c UPstream::opCodes::op_max
template<class T>
struct UPstream_opType<Foam::maxOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_max;
};

//- Map sumOp\<T\> to \c UPstream::opCodes::op_sum
template<class T>
struct UPstream_opType<Foam::sumOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_sum;
};

//- Map plusOp\<T\> (same as sumOp\<T\>) to \c UPstream::opCodes::op_sum
template<class T>
struct UPstream_opType<Foam::plusOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_sum;
};

//- Map multiplyOp\<T\> to \c UPstream::opCodes::op_prod
template<class T>
struct UPstream_opType<Foam::multiplyOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_prod;
};

// NOTE (2025-02):
// currently no mappings provided for
// (op_bool_and, op_bool_or, op_bool_xor) until the calling semantics
// have been properly defined


// These are only viable for unsigned integral types,
// probably not for signed integral types.
// Be extra restrictive for now

//- Map bitAndOp\<T\> to \c UPstream::opCodes::op_bit_and
//- for integrals (signed or unsigned), but also allow void as "generic"
template<class T>
struct UPstream_opType<Foam::bitAndOp<T>>
:
    std::bool_constant<std::is_integral_v<T> || std::is_void_v<T>>
{
    static constexpr auto opcode_id = []() constexpr noexcept
    {
        if constexpr (std::is_integral_v<T> || std::is_void_v<T>)
            return UPstream::opCodes::op_bit_and;
        else
            return UPstream::opCodes::invalid;
    }();
};

//- Map bitOrOp\<T\> to \c UPstream::opCodes::op_bit_or
//- for integrals (signed or unsigned), but also allow void as "generic"
template<class T>
struct UPstream_opType<Foam::bitOrOp<T>>
:
    std::bool_constant<std::is_integral_v<T> || std::is_void_v<T>>
{
    static constexpr auto opcode_id = []() constexpr noexcept
    {
        if constexpr (std::is_integral_v<T> || std::is_void_v<T>)
            return UPstream::opCodes::op_bit_or;
        else
            return UPstream::opCodes::invalid;
    }();
};

//- Map bitXorOp\<T\> to \c UPstream::opCodes::op_bit_xor
//- for integrals (signed or unsigned), but also allow void as "generic"
template<class T>
struct UPstream_opType<Foam::bitXorOp<T>>
:
    std::bool_constant<std::is_integral_v<T> || std::is_void_v<T>>
{
    static constexpr auto opcode_id = []() constexpr noexcept
    {
        if constexpr (std::is_integral_v<T> || std::is_void_v<T>)
            return UPstream::opCodes::op_bit_xor;
        else
            return UPstream::opCodes::invalid;
    }();
};


//- Combined query of opType and the underlying basic data type
//  This handling may be simplified in the future...
template<class BinaryOp, class T>
struct UPstream_data_opType
:
    std::bool_constant
    <
        UPstream_opType<BinaryOp>::value
     && UPstream_basic_dataType<T>::value
    >
{
    static constexpr auto opcode_id = []() constexpr noexcept
    {
        if constexpr
        (
            UPstream_opType<BinaryOp>::value
         && UPstream_basic_dataType<T>::value
        )
            return UPstream_opType<BinaryOp>::opcode_id;
        else
            return UPstream::opCodes::invalid;
    }();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convenience Functions (FUTURE?)

// inline bool is_UPstream_mpi_dataTypeCode(UPstream::dataTypes id) noexcept
// {
//     return
//     (
//         int(id) >= int(UPstream::opCodes::Basic_begin)
//      && int(id)  < int(UPstream::opCodes::Basic_end)
//     );
// }
//
// inline bool is_UPstream_reduceOpCode(UPstream::opCodes id) noexcept
// {
//     return
//     (
//         int(id) >= int(UPstream::opCodes::Basic_begin)
//      && int(id)  < int(UPstream::opCodes::Basic_end)
//     );
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
