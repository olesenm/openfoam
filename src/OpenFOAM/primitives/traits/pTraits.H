/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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
    Foam::pTraits

Description
    A traits class, which is primarily used for primitives and vector-space.

    All primitives need a specialised version of this class. The
    specialised versions will normally also require a conversion
    method.

\*---------------------------------------------------------------------------*/

#ifndef Foam_pTraits_H
#define Foam_pTraits_H

#include "direction.H"
#include <memory>       // For addressof
#include <type_traits>  // For std::integral_constant, std::void_t
#include <utility>      // For declval

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class zero_one Declaration
\*---------------------------------------------------------------------------*/

//- Represents 0/1 range or concept. Used for tagged dispatch or clamping
class zero_one {};


/*---------------------------------------------------------------------------*\
                           Class pTraits Declaration
\*---------------------------------------------------------------------------*/

// The base implementation is a pass-through to the base class.
// Accordingly it inherits all static methods (eg, typeName etc).
template<class Base>
class pTraits
:
    public Base
{
public:

    // Constructors

        //- Copy construct from base class
        explicit pTraits(const Base& obj)
        :
            Base(obj)
        {}

        //- Construct from Istream
        explicit pTraits(Istream& is)
        :
            Base(is)
        {}
};


/*---------------------------------------------------------------------------*\
                         Robuster generalised pTraits
\*---------------------------------------------------------------------------*/

//- The vector-space rank: default is 0.
template<class T, class = void>
inline constexpr Foam::direction pTraits_rank_v = 0;

//- The vector-space rank:
//- using the pTraits \c rank static member.
template<class T>
inline constexpr Foam::direction pTraits_rank_v
<
    T,
    std::void_t<decltype(pTraits<std::remove_cv_t<T>>::rank)>
> = pTraits<std::remove_cv_t<T>>::rank;

//- Number of data or vector-space components:
//- default is 1.
template<class T, class = void>
struct pTraits_nComponents : std::integral_constant<Foam::direction, 1> {};

//- Number of data or vector-space components:
//- using the pTraits \c nComponents static member.
//
//  Note (2026-08-06):
//  slight inconsistency here. It tests the base type for nComponents
//  but uses nComponents from pTraits (may need further review)
template<class T>
struct pTraits_nComponents
<
    T,
    std::void_t<decltype(std::remove_cv_t<T>::nComponents)>
>
:
    std::integral_constant
    <
        Foam::direction,
        pTraits<std::remove_cv_t<T>>::nComponents
    >
{};

//- Number of data or vector-space components
template<class T>
inline constexpr Foam::direction pTraits_nComponents_v
    = pTraits_nComponents<T>::value;

//- The underlying component data type: default is pass-through.
template<class T, class = void>
struct pTraits_cmptType { typedef T type; };

//- The underlying component data type for vector-space (or complex).
//
//  Enabled when pTraits<T>:zero has a defined type (ie, exists),
//  such as for arithmetic primitives, vector-space etc. where the concept
//  of a component type also makes sense.
template<class T>
struct pTraits_cmptType
<
    T,
    std::void_t<decltype(pTraits<std::remove_cv_t<T>>::zero)>
>
{
    using type = typename pTraits<std::remove_cv_t<T>>::cmptType;
};


//- Test for pTraits zero : default is false
template<class T, class = void>
inline constexpr bool pTraits_has_zero = false;

//- Test for pTraits zero
template<class T>
inline constexpr bool pTraits_has_zero
<
    T,
    std::void_t<decltype(pTraits<std::remove_cv_t<T>>::zero)>
> = true;


//- Test for pTraits typeName : default is false
template<class T, class = void>
inline constexpr bool pTraits_has_typeName = false;

//- Test for pTraits typeName
template<class T>
inline constexpr bool pTraits_has_typeName
<
    T,
    std::void_t<decltype(pTraits<std::remove_cv_t<T>>::typeName)>
> = true;


/*---------------------------------------------------------------------------*\
                              VectorSpace Traits
\*---------------------------------------------------------------------------*/

//- Test for VectorSpace : default is false
template<class T, class = void>
struct is_vectorspace : std::false_type {};

//- Test for VectorSpace : test for T::rank != 0 static member directly
//  Do not need pTraits layer since rank is defined via VectorSpace etc
template<class T>
struct is_vectorspace
<
    T,
    std::void_t<decltype(std::remove_cv_t<T>::rank)>
>
:
    std::bool_constant<std::remove_cv_t<T>::rank != 0>
{};


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


//- The is_vectorspace value of Type
template<class T>
inline constexpr bool is_vectorspace_v
    = is_vectorspace<T>::value;

//- The is_rotational_vectorspace value of Type
template<class T>
inline constexpr bool is_rotational_vectorspace_v
    = is_rotational_vectorspace<T>::value;


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


/*---------------------------------------------------------------------------*\
                               Container Traits
\*---------------------------------------------------------------------------*/

//- Test for containers with begin/end range iterators.
template<class T, class = void>
struct is_range : std::false_type {};

//- Test for list containers with begin/end range iterators.
//- \attention May yield a false positive with HashTable, HashSet etc
template<class T>
struct is_range
<
    T,
    std::void_t
    <
        decltype(std::declval<T>().begin()),
        decltype(std::declval<T>().end())
    >
>
:
    std::true_type
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
