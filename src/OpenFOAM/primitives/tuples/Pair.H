/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::Pair

Description
    An ordered pair of two objects of type \<T\> with first() and second()
    elements.

SourceFiles
    PairI.H

See also
    Foam::Tuple2 for storing two objects of dissimilar types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_Pair_H
#define Foam_Pair_H

#include "FixedList.H"
#include "Istream.H"
#include "Ostream.H"

#include <cstddef>
#include <tuple>
#include <utility>  // For std::move, std::pair

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class Pair;

template<class T> Istream& operator>>(Istream&, Pair<T>&);
template<class T> Ostream& operator<<(Ostream&, const Pair<T>&);

// Common pair types
typedef Pair<label> labelPair;      //!< A pair of labels
typedef Pair<word> wordPair;        //!< A pair of words


/*---------------------------------------------------------------------------*\
                            Class Pair Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Pair
:
    public FixedList<T, 2>
{
    // Typedefs

        //- The data content type
        typedef FixedList<T, 2> tuple_type;

public:

    // Generated Methods

        //- Default construct
        Pair() = default;

        //- Copy construct
        constexpr Pair(const Pair&) = default;

        //- Move construct
        constexpr Pair(Pair&&) = default;

        //- Copy assignment
        constexpr Pair& operator=(const Pair&) = default;

        //- Move assignment
        constexpr Pair& operator=(Pair&&) = default;

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Copy construct from components
        inline constexpr Pair(const T& f, const T& s);

        //- Move construct from components
        inline constexpr Pair(T&& f, T&& s);

        //- Copy construct from std::pair
        inline constexpr Pair(const std::pair<T,T>& vals);

        //- Move construct from std::pair
        inline constexpr Pair(std::pair<T,T>&& vals);

        //- Copy construct FixedList of two items
        inline constexpr Pair(const FixedList<T, 2>& list);

        //- Copy construct, optionally sorted with first less-than second
        inline constexpr Pair(const T& f, const T& s, bool doSort);

        //- Copy construct, optionally sorted with first less-than second
        inline constexpr Pair(const FixedList<T, 2>& list, bool doSort);

        //- Construct from Istream, as (first second) components
        explicit Pair(Istream& is) : FixedList<T, 2>() { is >> *this; }


    // Member Functions

    // Access

        //- Access the first element
        constexpr const T& first() const noexcept
        {
            return this->template get<0>();
        }

        //- Access the first element
        constexpr T& first() noexcept { return this->template get<0>(); }

        //- Access the second element
        constexpr const T& second() const noexcept
        {
            return this->template get<1>();
        }

        //- Access the second element
        constexpr T& second() noexcept { return this->template get<1>(); }

        //- Return other element
        inline const T& other(const T& a) const;


    // Queries

        //- True if first() is less-than-equal second()
        inline bool is_sorted() const;


    // Editing

        //- Flip the Pair in-place.
        inline void flip();

        //- Sort so that first() is less-than second()
        inline void sort();


    // Comparison

        //- Compare Pairs
        //  \return
        //  -  0: different
        //  - +1: identical values and order used
        //  - -1: identical values, but in reversed order
        static inline int compare(const Pair<T>& a, const Pair<T>& b);


    // Member Operators

        //- Compare for equality
        bool operator==(const Pair& other) const
        {
            return
            (
                static_cast<const tuple_type&>(*this)
             == static_cast<const tuple_type&>(other)
            );
        }

        //- Compare lexicographically
        bool operator<(const Pair& other) const
        {
            return
            (
                static_cast<const tuple_type&>(*this)
              < static_cast<const tuple_type&>(other)
            );
        }

        bool operator!=(const Pair& other) const { return !(*this == other); }
        bool operator<=(const Pair& other) const { return !(other < *this); }
        bool operator>(const Pair& other) const { return (other < *this); }
        bool operator>=(const Pair& other) const { return !(*this < other); }


    // Hashing

        //- Symmetric hashing functor for Pair, hashes lower value first
        //  Regular hasher inherited from FixedList
        struct symmHasher
        {
            unsigned operator()(const Pair<T>& obj, unsigned seed=0) const
            {
                Foam::Hash<T> op;
                if (obj.second() < obj.first())
                {
                    return op(obj.first(), op(obj.second(), seed));
                }
                else
                {
                    return op(obj.second(), op(obj.first(), seed));
                }
            }
        };
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return reverse of a Pair
template<class T>
Pair<T> reverse(const Pair<T>& p)
{
    return Pair<T>(p.second(), p.first());
}

// Note: Traits and IO Operators defined later in the included inline file...

// * * * * * * * * * * * * * * * ADL Lookups * * * * * * * * * * * * * * * * //

template<std::size_t I, class T>
constexpr T& get(Pair<T>& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Pair");
    return Foam::get<I>(static_cast<FixedList<T, 2>&>(t));
}

template<std::size_t I, class T>
constexpr const T& get(const Pair<T>& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Pair");
    return Foam::get<I>(static_cast<const FixedList<T, 2>&>(t));
}

template<std::size_t I, class T>
constexpr T&& get(Pair<T>&& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Pair");
    return Foam::get<I>(static_cast<FixedList<T, 2>&&>(std::move(t)));
}

template<std::size_t I, class T>
constexpr const T&& get(const Pair<T>&& t) noexcept
{
    static_assert(I < 2, "Index out-of-bounds for Pair");
    return Foam::get<I>(static_cast<const FixedList<T, 2>&&>(std::move(t)));
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations (structured binding)
namespace std
{

//- The tuple_size of Pair == 2
template<class T>
struct tuple_size<Foam::Pair<T>> : std::integral_constant<std::size_t, 2> {};

//- Identical type for all elements of Pair
template<std::size_t I, class T>
struct tuple_element<I, Foam::Pair<T>> { using type = T; };

} // End namespace std


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
