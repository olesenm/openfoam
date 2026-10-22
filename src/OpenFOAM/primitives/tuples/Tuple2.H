/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::Tuple2

Description
    A 2-tuple for storing two objects of dissimilar types.
    The container is similar in purpose to std::pair, but does not expose
    its members directly.
    The 2-tuple is considered contiguous if members are both contiguous.

See also
    Foam::Pair for storing two objects of identical types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_Tuple2_H
#define Foam_Tuple2_H

#include "Istream.H"
#include "Ostream.H"
#include "Pair.H"  // includes <tuple>, <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T1, class T2> class Tuple2;

template<class T1, class T2>
Istream& operator>>(Istream&, Tuple2<T1,T2>&);

template<class T1, class T2>
Ostream& operator<<(Ostream&, const Tuple2<T1,T2>&);


/*---------------------------------------------------------------------------*\
                           Class Tuple2 Declaration
\*---------------------------------------------------------------------------*/

template<class T1, class T2 = T1>
class Tuple2
{
    // Typedefs

        //- The data content type
        typedef std::pair<T1, T2> tuple_type;


    // Private Data

        //- The data content
        std::pair<T1, T2> tuple_;

public:

    // Typedefs (cf. std::pair)

        //- Type of member first, the first template parameter (T1)
        typedef T1 first_type;

        //- Type of member second, the second template parameter (T2)
        typedef T2 second_type;


    // Generated Methods

        //- Default construct
        Tuple2() = default;

        //- Copy construct
        constexpr Tuple2(const Tuple2&) = default;

        //- Move construct
        constexpr Tuple2(Tuple2&&) = default;

        //- Copy assignment
        Tuple2& operator=(const Tuple2&) = default;

        //- Move assignment
        Tuple2& operator=(Tuple2&&) = default;


    // Constructors

        //- Copy construct from components
        constexpr Tuple2(const T1& f, const T2& s)
        :
            tuple_(f, s)
        {}

        //- Move construct from components
        constexpr Tuple2(T1&& f, T2&& s)
        :
            tuple_(std::move(f), std::move(s))
        {}

        //- Copy construct from std::pair
        constexpr Tuple2(const std::pair<T1,T2>& tup)
        :
            tuple_(tup)
        {}

        //- Move construct from std::pair
        constexpr Tuple2(std::pair<T1,T2>&& tup)
        :
            tuple_(std::move(tup))
        {}

        //- Construct from Istream, as (first second) components
        explicit Tuple2(Istream& is) { is >> *this; }


    // Member Functions

        //- Access the first element
        constexpr const T1& first() const noexcept { return tuple_.first; }

        //- Access the first element
        constexpr T1& first() noexcept { return tuple_.first; }

        //- Access the second element
        constexpr const T2& second() const noexcept { return tuple_.second; }

        //- Access the second element
        constexpr T2& second() noexcept { return tuple_.second; }


    // Member Operators

        //- Compare for equality
        bool operator==(const Tuple2& other) const
        {
            return (this->tuple_ == other.tuple_);
        }

        #if __cpp_lib_three_way_comparison
        auto operator<=>(const Tuple2& other) const
        {
            return (this->tuple_ <=> other.tuple_);
        }
        #else  /* __cpp_lib_three_way_comparison */

        //- Compare lexicographically
        bool operator<(const Tuple2& other) const
        {
            return (this->tuple_ < other.tuple_);
        }

        bool operator!=(const Tuple2& other) const { return !(*this == other); }
        bool operator<=(const Tuple2& other) const { return !(other < *this); }
        bool operator>(const Tuple2& other) const { return (other < *this); }
        bool operator>=(const Tuple2& other) const { return !(*this < other); }
        #endif  /* __cpp_lib_three_way_comparison */
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

namespace ListPolicy
{

//- No compoundToken for List of arbitrary Tuple2,
//- so no implicit binary file IO (can still use MPI send/recv etc.)
template<class T1, class T2>
struct allow_contiguous<Tuple2<T1, T2>> : std::false_type {};

} // End namespace ListPolicy


//- Tuple2 is contiguous if both of its components are contiguous
template<class T1, class T2>
struct is_contiguous<Tuple2<T1, T2>>
:
    std::bool_constant<is_contiguous_v<T1> && is_contiguous_v<T2>>
{};

//- Is Tuple2 composed of labels
template<class T1, class T2>
struct is_contiguous_label<Tuple2<T1, T2>>
:
    std::bool_constant
    <is_contiguous_label<T1>::value && is_contiguous_label<T2>::value>
{};

//- Is Tuple composed of scalars
template<class T1, class T2>
struct is_contiguous_scalar<Tuple2<T1, T2>>
:
    std::bool_constant
    <is_contiguous_scalar<T1>::value && is_contiguous_scalar<T2>::value>
{};

//- If Tuple has common base element type
template<class T1, class T2>
struct contiguous_element<Tuple2<T1, T2>>
{
    using type = std::conditional_t
    <
        std::is_same_v<contiguous_element_t<T1>, contiguous_element_t<T2>>,
        contiguous_element_t<T1>,
        Tuple2<T1, T2>  // Fallback: pass-through
    >;
};

//- Hashing for Tuple2 data
template<class T1, class T2>
struct Hash<Tuple2<T1, T2>>
{
    unsigned operator()(const Tuple2<T1, T2>& obj, unsigned seed=0) const
    {
        if constexpr (is_contiguous_v<Tuple2<T1, T2>>)
        {
            return Foam::Hasher(&obj, sizeof(obj), seed);
        }
        else
        {
            return Hash<T2>()(obj.second(), Hash<T1>()(obj.first(), seed));
        }
    }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read Tuple2 from Istream
template<class T1, class T2>
inline Istream& operator>>(Istream& is, Tuple2<T1,T2>& t)
{
    is.readBegin("Tuple2");
    is >> t.first() >> t.second();
    is.readEnd("Tuple2");

    is.check(FUNCTION_NAME);
    return is;
}


//- Write Tuple2 to Ostream
template<class T1, class T2>
inline Ostream& operator<<(Ostream& os, const Tuple2<T1,T2>& t)
{
    os  << token::BEGIN_LIST
        << t.first() << token::SPACE << t.second()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return reverse of a Tuple2
template<class T1, class T2>
inline Tuple2<T2, T1> reverse(const Tuple2<T1,T2>& t)
{
    return Tuple2<T2, T1>(t.second(), t.first());
}

// * * * * * * * * * * * * * * * Comparators * * * * * * * * * * * * * * * * //

// Comparing first only

//- Compare tuple-like containers
//  \return reference to the container with the smaller value of first
struct minFirstOp
{
    template<class T>
    const Pair<T>& operator()(const Pair<T>& a, const Pair<T>& b) const
    {
        return (b.first() < a.first()) ? b : a;
    }

    template<class T1, class T2>
    const Tuple2<T1,T2>&
    operator()(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b) const
    {
        return (b.first() < a.first()) ? b : a;
    }

    template<class T1, class T2>
    const std::pair<T1,T2>&
    operator()(const std::pair<T1,T2>& a, const std::pair<T1,T2>& b) const
    {
        return (b.first < a.first) ? b : a;
    }
};


//- Assign tuple-like container to use the one with the smaller value of first
struct minFirstEqOp
{
    template<class T>
    void operator()(Pair<T>& x, const Pair<T>& y) const
    {
        if (y.first() < x.first()) x = y;
    }

    template<class T1, class T2>
    void operator()(Tuple2<T1,T2>& x, const Tuple2<T1,T2>& y) const
    {
        if (y.first() < x.first()) x = y;
    }

    template<class T1, class T2>
    void operator()(std::pair<T1,T2>& x, const std::pair<T1,T2>& y) const
    {
        if (y.first < x.first) x = y;
    }
};


//- Compare tuple-like containers
//  \return reference to the container with the larger value of first
struct maxFirstOp
{
    template<class T>
    const Pair<T>& operator()(const Pair<T>& a, const Pair<T>& b) const
    {
        return (a.first() < b.first()) ? b : a;
    }

    template<class T1, class T2>
    const Tuple2<T1,T2>&
    operator()(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b) const
    {
        return (a.first() < b.first()) ? b : a;
    }

    template<class T1, class T2>
    const std::pair<T1,T2>&
    operator()(const std::pair<T1,T2>& a, const std::pair<T1,T2>& b) const
    {
        return (a.first < b.first) ? b : a;
    }
};


//- Assign tuple-like container to use the one with the larger value of first
struct maxFirstEqOp
{
    template<class T>
    void operator()(Pair<T>& x, const Pair<T>& y) const
    {
        if (x.first() < y.first()) x = y;
    }

    template<class T1, class T2>
    void operator()(Tuple2<T1,T2>& x, const Tuple2<T1,T2>& y) const
    {
        if (x.first() < y.first()) x = y;
    }

    template<class T1, class T2>
    void operator()(std::pair<T1,T2>& x, const std::pair<T1,T2>& y) const
    {
        if (x.first < y.first) x = y;
    }
};


// * * * * * * * * * * * * * * * ADL Lookups * * * * * * * * * * * * * * * * //

template<std::size_t I, class T1, class T2>
constexpr decltype(auto) get(Tuple2<T1, T2>& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Tuple2");
    if constexpr (I == 0) { return t.first(); }
    else { return t.second(); }
}

template<std::size_t I, class T1, class T2>
constexpr decltype(auto) get(const Tuple2<T1, T2>& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Tuple2");
    if constexpr (I == 0) { return t.first(); }
    else { return t.second(); }
}

template<std::size_t I, class T1, class T2>
constexpr decltype(auto) get(Tuple2<T1, T2>&& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Tuple2");
    if constexpr (I == 0) { return std::move(t.first()); }
    else { return std::move(t.second()); }
}

template<std::size_t I, class T1, class T2>
constexpr decltype(auto) get(const Tuple2<T1, T2>&& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Tuple2");
    if constexpr (I == 0) { return std::move(t.first()); }
    else { return std::move(t.second()); }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations (structured binding)
namespace std
{

//- The tuple_size of Tuple2 == 2
template<class T1, class T2>
struct tuple_size<Foam::Tuple2<T1, T2>>
: std::integral_constant<std::size_t, 2> {};

template<std::size_t I, class T1, class T2>
struct tuple_element<I, Foam::Tuple2<T1, T2>>
{
    static_assert(I < 2, "Index out-of-bounds for Tuple2");
};

template<class T1, class T2>
struct tuple_element<0, Foam::Tuple2<T1, T2>> { using type = T1; };

template<class T1, class T2>
struct tuple_element<1, Foam::Tuple2<T1, T2>> { using type = T2; };

} // End namespace std


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
