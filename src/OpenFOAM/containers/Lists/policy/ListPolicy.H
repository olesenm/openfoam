/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    Foam::ListPolicy

Description
    Additional compile-time controls of List behaviour

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListPolicy_H
#define Foam_ListPolicy_H

#include "MemoryPool.H"  // Also includes <cstddef>
#include "contiguous.H"  // Also includes <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class keyType;
class word;
class wordRe;

namespace ListPolicy
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Number of items before requiring line-breaks in the list output.
//
//  Default definition: 10
template<class T>
struct short_length : std::integral_constant<int,10> {};

// Can override on a per-type basis
// Eg,
// template<> struct short_length<label> : std::integral_constant<int,20> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Can suppress additional line breaks separate ASCII data content
//- when the data elements are primitives, or contiguous
//
//  Default definition: (integral | floating-point) are contiguous and thus
//  never need any line breaks
template<class T>
struct no_linebreak : std::is_arithmetic<std::remove_cv_t<T>> {};

// Specialization for word-like classes
// These elements are normally fairly short, so ok to output a few (eg, 10)
// of them on a single line.

//- Suppress line-breaks for keyType
template<> struct no_linebreak<keyType> : std::true_type {};

//- Suppress line-breaks for word
template<> struct no_linebreak<word> : std::true_type {};

//- Suppress line-breaks for wordRe
template<> struct no_linebreak<wordRe> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Memory allocation/deallocation handling - primarily used by List and Matrix
//
// - is_aligned_type() :
//   Defines which types may be aligned
//
// - use_alignment(n) :
//   Lower threshold for using memory alignment
//
// - use_memory_pool(n) :
//   Lower threshold for using a memory pool.
//   Must be larger than use_alignment() value.
//
// - use_offload(n) :
//   Lower threshold for switching to device offloading
//
//
// Use of the memory-pool is controlled by the 'is_aligned_type()' test
// and the minimum field size, controlled by the 'use_memory_pool()' test.
//
// If the memory-pool is not enabled or not required according to the two
// above tests, the allocation falls back to either an aligned or unaligned
// allocation.
//
// The decision about when to choose aligned vs unaligned allocation
// is still a compile-time option. Made by direct edit of the
// appropriate functions.
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Consider aligned allocation for the given type?
//  Benefits for field data (floating-point, ints, vectorspace),
//  but avoid for char data, strings, pointers etc
template<class T>
inline constexpr bool is_aligned_type() noexcept
{
    #if defined(_WIN32)
    // Placement alignment seems to fail on MINGW (#3488), so disable there.
    // Note that returning false also effectively disables MemoryPool too,
    // which is not currently supported on MINGW anyhow.
    return false;
    #else
    return
    (
        !std::is_enum_v<T>
     && !std::is_pointer_v<T>
     && !std::is_union_v<T>
     && (sizeof(T) >= 4)    // skip small data (eg, char)
     && is_contiguous_v<T>
    );
    #endif
}


//- True if size exceeds the min-size for using memory alignment
template<class IntType>
inline constexpr bool use_alignment(IntType n) noexcept
{
    return (n >= IntType(200));
}


//- True if size exceeds the min-size for using the memory pool
template<class IntType>
inline constexpr bool use_memory_pool(IntType n) noexcept
{
    return (n >= IntType(3000));
}


//- True if size exceeds the min-size for offloading
template<class IntType>
inline constexpr bool use_offload(IntType n) noexcept
{
    return (n >= IntType(1000));
}


//- Default alignment for larger fields
inline constexpr std::align_val_t default_alignment() noexcept
{
    return std::align_val_t(256);
}


//- Allocate from memory pool (if active), or aligned, or normal
template<class T, class IntType>
inline T* allocate(IntType n)
{
    if constexpr (ListPolicy::is_aligned_type<T>())
    {
        // Note: threshold for use_memory_pool() >= use_alignment()

        if (ListPolicy::use_alignment(n))
        {
            if
            (
                void *pool_ptr
                (
                    // Consider memory pool for large amounts of data
                    ListPolicy::use_memory_pool(n)
                  ? Foam::MemoryPool::try_allocate(sizeof(T)*n)
                  : nullptr
                );
                pool_ptr
            )
            {
                // Placement new
                return new (pool_ptr) T[n];
            }
            else
            {
                return new (ListPolicy::default_alignment()) T[n];
            }
        }
        else
        {
            // Plain new
            return new T[n];
        }
    }
    else
    {
        // Plain new
        return new T[n];
    }
}


//- Deallocate from memory pool, or normal
template<class T, class IntType>
inline void deallocate(T* ptr)
{
    if constexpr (ListPolicy::is_aligned_type<T>())
    {
        if (ptr && !Foam::MemoryPool::try_deallocate(ptr))
        {
            // Plain new
            delete[] ptr;
        }
    }
    else
    {
        // Plain new
        delete[] ptr;
    }
}


//- Deallocate from memory pool, or aligned, or normal
template<class T, class IntType>
inline void deallocate(T* ptr, [[maybe_unused]] IntType n)
{
    if constexpr (ListPolicy::is_aligned_type<T>())
    {
        // Note: threshold for use_memory_pool() >= use_alignment()

        if (ListPolicy::use_alignment(n))
        {
            if (ptr && !Foam::MemoryPool::try_deallocate(ptr))
            {
                // Alignment depends on the number of elements
                ::operator delete[](ptr, ListPolicy::default_alignment());
            }
        }
        else
        {
            // Plain new
            delete[] ptr;
        }
    }
    else
    {
        // Plain new
        delete[] ptr;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Calculate a reserve size (eg, doubling) based on the requested length
//- and the current capacity
template<int SizeMin, int Numerator, class IntType>
inline IntType reserve_size(IntType requested, IntType capacity) noexcept
{
    static_assert(Numerator > 1, "Invalid numerator");

    // The caller already checks this:
    // if (requested < capacity) { return capacity; }

    IntType size(capacity*Numerator);
    if (size < requested)
    {
        size = requested;
    }
    if constexpr (SizeMin > 0)  // The min size is optional
    {
        if (size < SizeMin)
        {
            size = SizeMin;
        }
    }
    return size;
}


//- Calculate a reserve size based on the requested length
//- and the current capacity
template<int SizeMin, int Numerator, int Denominator, class IntType>
inline IntType reserve_size(IntType requested, IntType capacity) noexcept
{
    static_assert(Numerator > Denominator, "Invalid numerator");
    static_assert(Denominator > 0, "Invalid denominator");

    // The caller already checks this:
    // if (requested < capacity) { return capacity; }

    // Very unlikely that capacity is less than Denominator,
    // so divide before multiply to avoid overflow
    IntType size((capacity/Denominator)*Numerator);
    if (size < requested)
    {
        size = requested;
    }
    if constexpr (SizeMin > 0)  // The min size is optional
    {
        if (size < SizeMin)
        {
            size = SizeMin;
        }
    }
    return size;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Classification of list/container uniformity.
//- The values can be used with bit-wise \c or reduction
enum uniformity : unsigned char
{
    EMPTY = 0,          //!< An empty container
    UNIFORM = 0x1,      //!< Container (non-empty) with identical values
    NONUNIFORM = 0x2,   //!< Container (non-empty) with different values
    MIXED = 0x3         //!< Mixed uniform/non-uniform (eg, after reduction)
};


//- Algorithm to determine list/container uniformity
template<class InputIt>
enum uniformity check_uniformity(InputIt first, InputIt last)
{
    if (first == last) return uniformity::EMPTY;

    // Like std::all_of() with checking against element 0,
    // but without using a lambda with auto type (pre C++14) etc.

    const auto& elem0 = *first;

    for ((void)++first; (first != last); (void)++first)
    {
        if (elem0 != *first)
        {
            return uniformity::NONUNIFORM;
        }
    }

    return uniformity::UNIFORM;
}


//- Checking for List uniformity requires the "operator!=()" :
template<class T, class = void>
struct has_uniformity : std::false_type {};

template<class T>
struct has_uniformity
<
    T,
    std::void_t<decltype(std::declval<const T&>() != std::declval<const T&>())>
> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ListPolicy
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
