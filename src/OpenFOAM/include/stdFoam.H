/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2025 OpenCFD Ltd.
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
    Includes some common C++ headers, defines global macros and templates
    used in multiple places by OpenFOAM.

\*---------------------------------------------------------------------------*/

#ifndef Foam_stdFoam_H
#define Foam_stdFoam_H

#include <algorithm>
#include <initializer_list>
#include <iterator>  // for std::begin, std::end, ...
#include <memory>    // for addressof etc
#include <numeric>   // for std::iota, std::reduce, ...
#include <type_traits>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compile-time warning for use of deprecated methods (compiler-dependent).
// Use within the class declaration.

#define FOAM_DEPRECATED(since) [[deprecated("Since " #since)]]
#define FOAM_DEPRECATED_FOR(since, replacement) [[deprecated("Since " #since "; use " #replacement)]]

#ifdef FOAM_COMPILE_STRICT
#  define FOAM_DEPRECATED_STRICT(since, replacement) [[deprecated("Since " #since "; use " #replacement)]]
#  if (FOAM_COMPILE_STRICT > 1)
#  define FOAM_DEPRECATED_STRICTER(since, replacement) [[deprecated("Since " #since "; use " #replacement)]]
#  endif
#endif

#ifndef FOAM_DEPRECATED_STRICT
#define FOAM_DEPRECATED_STRICT(since, replacement)
#endif
#ifndef FOAM_DEPRECATED_STRICTER
#define FOAM_DEPRECATED_STRICTER(since, replacement)
#endif

// Branch prediction helpers. With C++20 can use [[likely]], [[unlikely]]
#if defined(__GNUC__) || defined(__clang__)
#  define FOAM_UNLIKELY(cond) __builtin_expect(!!(cond),0)
#  define FOAM_LIKELY(cond)   __builtin_expect(!!(cond),1)
#else
#  define FOAM_UNLIKELY(cond) (cond)
#  define FOAM_LIKELY(cond)   (cond)
#endif

// Shadow macro for __GNUC__, excluding compilers masquerading as gcc
#undef FOAM_REAL_GNUC
#if defined(__GNUC__) && !defined(__llvm__)
#  define FOAM_REAL_GNUC __GNUC__
#endif

// Suppress false positives from -Wdangling-reference (gcc >= 14)
#if (FOAM_REAL_GNUC >= 14)
#  define FOAM_NO_DANGLING_REFERENCE [[gnu::no_dangling]]
#endif

#ifndef FOAM_NO_DANGLING_REFERENCE
#define FOAM_NO_DANGLING_REFERENCE
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Namespace for OpenFOAM
namespace Foam
{

//- Implementation details for various OpenFOAM classes
namespace Detail {}

//- Additional OpenFOAM modules
namespace Module {}


//- A functor that returns its argument unchanged (cf. C++20 std::identity)
//- Should \em never be specialized.
struct identityOp
{
    using is_transparent = void;

    template<class T>
    constexpr T&& operator()(T&& val) const noexcept
    {
        return std::forward<T>(val);
    }

    // Allow use as an identity array/map
    template<class T>
    constexpr T&& operator[](T&& val) const noexcept
    {
        return std::forward<T>(val);
    }
};


// Dispatch tags

//- Dispatch tag type, the opposite of piecewise_construct
struct direct_construct_t { explicit direct_construct_t() = default; };

//- Dispatch tag type
struct identity_construct_t { explicit identity_construct_t() = default; };

//- Dispatch tag type (eg, for tuple types)
struct piecewise_construct_t { explicit piecewise_construct_t() = default; };

//- Dispatch tag for direct construct
constexpr direct_construct_t direct_construct{};

//- Dispatch tag for identity construct
constexpr identity_construct_t identity_construct{};

//- Dispatch tag for piecewise construct
constexpr piecewise_construct_t piecewise_construct{};

//- Represents 0/1 range or concept. Used for tagged dispatch or clamping
class zero_one {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Swap non-array types as per std::swap, but in Foam namespace.
//  \sa http://www.cplusplus.com/reference/utility/swap/
//
//  \note For complex structures, it is usual to provide a swap() member
//  function and specialize Swap()
template<class T>
void Swap(T& a, T& b)
{
    std::swap(a, b);
}


//- Swap array types as per std::swap example, but in Foam namespace.
//  \sa http://www.cplusplus.com/reference/utility/swap/
template<class T, size_t N>
void Swap(T (&a)[N], T (&b)[N])
{
    for (size_t i = 0; i < N; ++i)
    {
        Foam::Swap(a[i], b[i]);
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Namespace for std templates that are are part of future C++ standards
//- or that are in a state of change.
//
// SeeAlso
//     - https://en.cppreference.com/w/cpp/iterator/begin
//     - https://en.cppreference.com/w/cpp/iterator/end
//     - https://en.cppreference.com/w/cpp/iterator/rbegin
//     - https://en.cppreference.com/w/cpp/iterator/rend

namespace stdFoam
{

//- Map any dependent type to false (workaround before CWG2518)
template<typename...>
inline constexpr bool dependent_false_v = false;

//- Test for \c bool type (after removing cv)
template<typename T>
inline constexpr bool is_bool_v
    = std::is_same<bool, std::remove_cv_t<T>>::value;

//- Test for \c char type (after removing cv)
template<typename T>
inline constexpr bool is_char_v
    = std::is_same<char, std::remove_cv_t<T>>::value;

#if defined(__cpp_lib_remove_cvref) && (__cpp_lib_remove_cvref >= 201711L)
//- Forward to c++20 std::remove_cvref_t
template<class T>
using remove_cvref_t = std::remove_cvref_t<T>;
#else
//- C++17 fallback for std::remove_cvref_t (c++20)
template<class T>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
#endif

//- Provides the member typedef type that names T
//- (i.e., the identity transformation).
//  Available as C++20 std::type_identity
//  (__cpp_lib_type_identity >= 201806L)
template<class T>
struct type_identity
{
    using type = T;
};

//- The \c type from type_identity
//  Available as C++20 std::type_identity_t
//  (__cpp_lib_type_identity >= 201806L)
template<class T>
using type_identity_t = typename type_identity<T>::type;


//- Test for a \c is_base_of relationship for templated base classes
template<template<class...> class Base, class Derived>
struct is_template_base_of
{
private:

    // Derived from Base => true
    template<class... Args>
    static std::true_type test(const Base<Args...>*);

    // Fallback => false
    static std::false_type test(...);

    using decayed = std::remove_cv_t<std::remove_reference_t<Derived>>;

public:

    // OR: decltype(test(std::declval<decayed*>()))::value
    static constexpr bool value
        = decltype(test(static_cast<decayed*>(nullptr)))::value;
};

//- Test for a \c is_base_of relationship for templated base classes
//- with one non-type template parameter
template<template<class, auto> class Base, class Derived>
struct is_template_nttp_base_of
{
private:

    // Derived from Base => true
    template<class T, auto N>
    static std::true_type test(const Base<T, N>*);

    // Fallback => false
    static std::false_type test(...);

    using decayed = std::remove_cv_t<std::remove_reference_t<Derived>>;

public:

    static constexpr bool value
        = decltype(test(static_cast<decayed*>(nullptr)))::value;
};

//- Test for \c is_base_of relationship for templated base classes
template<template<class...> class Base, class Derived>
inline constexpr bool
is_template_base_of_v = is_template_base_of<Base, Derived>::value;

//- Test for \c is_base_of relationship for templated base classes
//- with one non-type template parameter (NNTP)
template<template<class, auto> class Base, class Derived>
inline constexpr bool
is_template_nttp_base_of_v = is_template_nttp_base_of<Base, Derived>::value;

} // End namespace stdFoam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Iterate across all elements in the \a container object.
//  \par Usage
//  \code
//  forAllIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIters, forAllIter, forAllConstIters
#define forAllIters(container,iter)                                            \
    for                                                                        \
    (                                                                          \
        auto iter = std::begin(container);                                     \
        iter != std::end(container);                                           \
        ++iter                                                                 \
    )


//- Iterate across all elements of the \a container object with const access.
//  \par Usage
//  \code
//  forAllConstIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIters, forAllIter, forAllConstIter
#define forAllConstIters(container,iter)                                       \
    for                                                                        \
    (                                                                          \
        auto iter = std::cbegin(container);                                    \
        iter != std::cend(container);                                          \
        ++iter                                                                 \
    )


//- Reverse iterate across elements in the \a container object of type
//  \a Container.
//  \par Usage
//  \code
//  forAllReverseIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstReverseIters
#define forAllReverseIters(container,iter)                                     \
    for                                                                        \
    (                                                                          \
        auto iter = std::rbegin(container);                                    \
        iter != std::rend(container);                                          \
        ++iter                                                                 \
    )


//- Reverse iterate across elements of \a container object with const access.
//  \par Usage
//  \code
//  forAllReverseConstIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllReverseIters
#define forAllConstReverseIters(container,iter)                                \
    for                                                                        \
    (                                                                          \
        auto iter = std::crbegin(container);                                   \
        iter != std::crend(container);                                         \
        ++iter                                                                 \
    )


//- Loop across all elements in \a list
// \par Usage
// \code
// forAll(anyList, i)
// {
//      statements;
// }
// \endcode
// \sa forAllReverse
#define forAll(list, i) \
    for (Foam::label i=0; i<(list).size(); ++i)


//- Reverse loop across all elements in \a list
//  \par Usage
//  \code
//  forAllReverse(anyList, i)
//  {
//       statements;
//  }
//  \endcode
//  \sa forAll
#define forAllReverse(list, i) \
    for (Foam::label i=(list).size()-1; i>=0; --i)


// Compatibility macros for pre C++11

//- Iterate across all elements in the \a container object
//  of type \a Container.
//  \par Usage
//  \code
//  forAllIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIter
#define forAllIter(Container,container,iter)                                   \
    for                                                                        \
    (                                                                          \
        Container::iterator iter = (container).begin();                        \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )


//- Iterate across all elements in the \a container object
//  of type \a Container with const access.
//  \par Usage
//  \code
//  forAllConstIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIter
#define forAllConstIter(Container,container,iter)                              \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).cbegin();                 \
        iter != (container).cend();                                            \
        ++iter                                                                 \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace stdFoam
{

/*---------------------------------------------------------------------------*\
                        Class stdFoam::span Declaration
\*---------------------------------------------------------------------------*/

//- Rudimentary functionality similar to std::span for holding memory view
template<class Type>
class span
{
    // Private Data

        //- The data pointer
        Type* data_;

        //- The size of the span
        std::size_t size_;


public:

    // STL type definitions

        using element_type = Type;
        using value_type = std::remove_cv_t<Type>;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using pointer = Type*;
        using const_pointer = const Type*;
        using reference = element_type&;
        using const_reference = const element_type&;
        using iterator = Type*;


    // Generated Methods

        //- Copy construct
        span(const span& other) noexcept = default;

        //- Copy assignment
        span& operator=(const span& other) noexcept = default;


    // Constructors

        //- Default construct
        constexpr span() noexcept
        :
            data_(nullptr),
            size_(0)
        {}

        //- Construct from pointer and size
        constexpr span(pointer ptr, size_type count) noexcept
        :
            data_(ptr),
            size_(count)
        {}

        //- Construct from begin/end pointers
        constexpr span(pointer first, pointer last) noexcept
        :
            data_(first),
            size_(last - first)
        {}


    //- Destructor
    ~span() noexcept = default;


    // Member Functions

    // Iterators

        //- Iterator to begin of span
        constexpr iterator begin() const noexcept { return data_; }

        //- Iterator to one-past end of span
        constexpr iterator end() const noexcept { return (data_ + size_); }


    // Element access

        //- Access the first element. Undefined if span is empty
        constexpr reference front() const { return *(data_); }

        //- Access the last element. Undefined if span is empty
        constexpr reference back() const { return *(data_ + size_ - 1); }

        //- Access an element of the sequence
        constexpr reference operator[](size_type idx) const
        {
           return *(data_ + idx);
        }

        //- Return a pointer to the beginning of the sequence
        constexpr pointer data() const noexcept { return data_; }


    // Observers

        //- Number of elements in the sequence
        constexpr size_type size() const noexcept { return size_; }

        //- The size of the sequence in bytes
        constexpr size_type size_bytes() const noexcept
        {
            return (size_*sizeof(Type));
        }

        //- True if the sequence is empty
        constexpr bool empty() const noexcept { return !size_; }


    // Subview

        //- Obtains a span of the first count elements
        span<Type> first(size_type count) const noexcept
        {
            return span<Type>(data_, count);
        }

        //- Obtains a span of the last count elements
        span<Type> last(size_type count) const noexcept
        {
            return span<Type>(data_ + (size_ - count), count);
        }

        //- Obtains a sub-span starting at pos until end of the sequence
        span<Type> subspan(size_type pos) const noexcept
        {
            return span<Type>(data_ + pos, size_ - pos);
        }

        //- Obtains a sub-span of length len starting at pos.
        //  Graciously handles excess lengths.
        span<Type> subspan(size_type pos, size_type len) const noexcept
        {
            return span<Type>(data_ + pos, std::min(size_ - pos, len));
        }


    // Additional members, similar to UList etc.
    // The std::span has as_bytes() and as_writeable_bytes() as free functions

        //- A readonly view as byte content
        constexpr const char* cdata_bytes() const noexcept
        {
            return reinterpret_cast<const char*>(data_);
        }

        //- A writable view as byte content (if the pointer type is non-const).
        //- Like data(), the access itself is const.
        template
        <
            class TypeT = Type,
            class = std::enable_if_t<!std::is_const_v<TypeT>>
        >
        char* data_bytes() const noexcept
        {
            return reinterpret_cast<char*>(data_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stdFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
