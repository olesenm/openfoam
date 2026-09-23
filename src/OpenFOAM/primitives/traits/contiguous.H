/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::is_contiguous

Description
    A template class to specify that a data type can be considered as being
    contiguous in memory.

    Normally only integral and floating-point types can be considered
    contiguous, but some other types (eg, FixedList, Pair, Vector etc)
    consisting purely of these fundamental types can be considered
    as having a contiguous memory layout as well.

Note
    In OpenFOAM 1906 and earlier, the contiguous trait was handled
    by templated \c contiguous global functions.

    While possible to mark this as deleted, this does not detect or
    prevent specializations. Thus omit the usual housekeeping.

Class
    Foam::is_contiguous_label

Description
    A template class to specify if a data type is composed solely of
    Foam::label elements.

Class
    Foam::is_contiguous_scalar

Description
    A template class to specify if a data type is composed solely of
    Foam::scalar elements.

Class
    Foam::contiguous_element

Description
    A template class providing the underlying element type used when
    marshalling a contiguous type (eg, for MPI transfers) as a sequence
    of some fundamental type.

    The base definition resolves as follows:
    - arithmetic types : pass-through (self)
    - \c is_contiguous_scalar : \c Foam::scalar
    - \c is_contiguous_label  : \c Foam::label
    - otherwise: pass-through (self)
    .

Note
    All four of \c is_contiguous, \c is_contiguous_label,
    \c is_contiguous_scalar and \c contiguous_element provide additional
    generic definition for derived types with \c cmptType
    (eg vector-space: vector, tensor, etc)

\*---------------------------------------------------------------------------*/

#ifndef Foam_contiguous_H
#define Foam_contiguous_H

#include "labelFwd.H"
#include "scalarFwd.H"
#include <memory>    // for addressof etc
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                              VectorSpace Traits
\*---------------------------------------------------------------------------*/

//- Detect a VectorSpace : default is false
template<class T, class = void>
struct is_vectorspace : std::false_type {};

//- Detect a VectorSpace : tests for T::rank != 0 static member directly
//- which is characteristic of a vector-space.
//- If type does not have a \c rank member it is not a vector-space.
//- If type does have a \c rank member, but its value is 0, it is a
//- scalar quantity and not a vector-space quantity.
//
//  \note does not need the pTraits layer, since rank is defined via
//  VectorSpace etc
template<class T>
struct is_vectorspace<T, std::void_t<decltype(std::remove_cv_t<T>::rank)>>
:
    std::bool_constant<std::remove_cv_t<T>::rank != 0>
{};

//- The is_vectorspace value of Type
template<class T>
inline constexpr bool is_vectorspace_v
    = is_vectorspace<T>::value;


//- Pointer to component data (VectorSpace) or the address (regular types)
template<class T>
inline constexpr const auto* vectorspace_cdata(const T& obj) noexcept
{
    if constexpr (is_vectorspace_v<T>)
    {
        return obj.cdata();
    }
    else
    {
        return std::addressof(obj);
    }
}

//- Pointer to component data (VectorSpace) or the address (regular types)
template<class T>
inline constexpr auto* vectorspace_data(T& obj) noexcept
{
    if constexpr (is_vectorspace_v<T>)
    {
        return obj.data();
    }
    else
    {
        return std::addressof(obj);
    }
}


//- Test for VectorSpace and multiple components (eg, not sphericalTensor)
//- which means it is not rotationally invariant
template<class T, class = void>
struct is_rotational_vectorspace : std::false_type {};

//- Test for VectorSpace and multiple components (eg, not sphericalTensor)
//  Do not need pTraits layer since rank and nComponents are both
//  defined via VectorSpace
template<class T>
struct is_rotational_vectorspace
<
    T,
    std::void_t
    <
        decltype(std::remove_cv_t<T>::rank),
        decltype(std::remove_cv_t<T>::nComponents)
    >
>
:
    std::bool_constant
    <
        (std::remove_cv_t<T>::rank != 0)
     && (std::remove_cv_t<T>::nComponents > 1)
    >
{};


//- The is_rotational_vectorspace value of Type
template<class T>
inline constexpr bool is_rotational_vectorspace_v
    = is_rotational_vectorspace<T>::value;


/*---------------------------------------------------------------------------*\
                                  contiguous
\*---------------------------------------------------------------------------*/

// Forward Declarations
template<class T> struct contiguous_element;
template<class T> struct is_contiguous_label;
template<class T> struct is_contiguous_scalar;
template<class T> struct is_contiguous;

namespace Detail
{

// Dispatch helpers for the generic vector-space forwarding, keyed on a
// plain bool (is_vectorspace_v<T>) rather than a SFINAE second parameter
// of the public trait itself. This keeps the public traits as
// single-specialization templates, so a user-provided override never
// competes/ambiguity-clashes with the generic forwarding - even if the
// type later happens to add a 'rank' member.

template<class T, bool = is_vectorspace_v<T>>
struct is_contiguous_label_impl
    : std::is_same<label, std::remove_cv_t<T>> {};

template<class T>
struct is_contiguous_label_impl<T, true>
    : is_contiguous_label<typename T::cmptType> {};


template<class T, bool = is_vectorspace_v<T>>
struct is_contiguous_scalar_impl
    : std::is_same<scalar, std::remove_cv_t<T>> {};

template<class T>
struct is_contiguous_scalar_impl<T, true>
    : is_contiguous_scalar<typename T::cmptType> {};


template<class T, bool = is_vectorspace_v<T>>
struct is_contiguous_impl : std::is_arithmetic<T> {};

template<class T>
struct is_contiguous_impl<T, true>
:
    std::disjunction
    <
        is_contiguous_label<T>,
        is_contiguous_scalar<T>,
        is_contiguous<typename T::cmptType>
    >
{};


template<class T, bool = is_vectorspace_v<T>>
struct contiguous_element_impl
{
    using type = std::conditional_t
    <
        std::is_arithmetic_v<T>,
        T,
        std::conditional_t
        <
            is_contiguous_scalar<T>::value,
            scalar,
            std::conditional_t
            <
                is_contiguous_label<T>::value,
                label,
                T  // Fallback: pass-through
            >
        >
    >;
};

template<class T>
struct contiguous_element_impl<T, true>
{
    using type = typename contiguous_element<typename T::cmptType>::type;
};

} // End namespace Detail


//- A template class to specify if a data type is composed solely of
//- Foam::label elements.
//  Types can provide their own explicit specialization here (eg, when
//  they lack 'rank' and thus do not qualify for the generic forwarding).
template<class T>
struct is_contiguous_label
:
    Detail::is_contiguous_label_impl<std::remove_cv_t<T>>
{};


//- A template class to specify if a data type is composed solely of
//- Foam::scalar elements.
template<class T>
struct is_contiguous_scalar
:
    Detail::is_contiguous_scalar_impl<std::remove_cv_t<T>>
{};


//- A template class to specify that a data type can be considered as
//- being contiguous in memory.
template<class T>
struct is_contiguous
:
    Detail::is_contiguous_impl<std::remove_cv_t<T>>
{};

//- The is_contiguous value of Type (after stripping of qualifiers)
template<class T>
inline constexpr bool is_contiguous_v = is_contiguous<T>::value;


//- A template class providing the underlying element type used when
//- marshalling a contiguous type as a sequence of some fundamental type.
template<class T>
struct contiguous_element
:
    Detail::contiguous_element_impl<std::remove_cv_t<T>>
{};


//- The contiguous_element type
template<class T>
using contiguous_element_t = typename contiguous_element<T>::type;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
