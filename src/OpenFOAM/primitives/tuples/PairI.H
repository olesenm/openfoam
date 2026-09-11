/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
inline int Foam::Pair<T>::compare(const Pair<T>& a, const Pair<T>& b)
{
    if (a.first() == b.first() && a.second() == b.second())
    {
        return 1;
    }
    if (a.first() == b.second() && a.second() == b.first())
    {
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::Pair<T>::Pair(const T& f, const T& s)
:
    tuple_type(Foam::direct_construct_t{}, f, s)
{}


template<class T>
inline constexpr Foam::Pair<T>::Pair(T&& f, T&& s)
:
    tuple_type(Foam::direct_construct_t{}, std::move(f), std::move(s))
{}


template<class T>
inline constexpr Foam::Pair<T>::Pair(const std::pair<T,T>& vals)
:
    tuple_type(Foam::direct_construct_t{}, vals.first, vals.second)
{}


template<class T>
inline constexpr Foam::Pair<T>::Pair(std::pair<T,T>&& vals)
:
    tuple_type
    (
        Foam::direct_construct_t{},
        std::move(vals.first),
        std::move(vals.second)
    )
{}


template<class T>
inline constexpr Foam::Pair<T>::Pair(const FixedList<T, 2>& list)
:
    tuple_type(list)
{}


template<class T>
inline constexpr Foam::Pair<T>::Pair(const T& f, const T& s, bool doSort)
:
    tuple_type()
{
    if (doSort && s < f)
    {
        first()  = s;
        second() = f;
    }
    else
    {
        first()  = f;
        second() = s;
    }
}


template<class T>
inline constexpr Foam::Pair<T>::Pair(const FixedList<T, 2>& list, bool doSort)
:
    Pair<T>(list.template get<0>(), list.template get<1>(), doSort)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::Pair<T>::other(const T& val) const
{
    if (first() == second())
    {
        FatalErrorInFunction
            << "Call to other only valid for Pair with differing elements:"
            << *this << abort(FatalError);
    }
    else if (val == first())
    {
        return second();
    }
    else if (val != second())
    {
        FatalErrorInFunction
            << "Pair " << *this
            << " does not contain " << val << abort(FatalError);
    }

    return first();
}


template<class T>
inline void Foam::Pair<T>::flip()
{
    Foam::Swap(first(), second());
}


template<class T>
inline bool Foam::Pair<T>::is_sorted() const
{
    return !(second() < first());
}


template<class T>
inline void Foam::Pair<T>::sort()
{
    if (second() < first())
    {
        flip();
    }
}


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

namespace Foam
{

//- The base element for Pair
template<class T>
struct contiguous_element<Pair<T>> : contiguous_element<T> {};

//- Pair is contiguous if the type is contiguous
template<class T>
struct is_contiguous<Pair<T>> : is_contiguous<T> {};

//- Check for Pair of labels
template<class T>
struct is_contiguous_label<Pair<T>> : is_contiguous_label<T> {};

//- Check for Pair of scalars
template<class T>
struct is_contiguous_scalar<Pair<T>> : is_contiguous_scalar<T> {};

//- Hashing for Pair of data
template<class T> struct Hash<Pair<T>> : Pair<T>::hasher {};

} // End namespace Foam


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

namespace Foam
{

//- std::pair is contiguous if both of its components are contiguous
template<class T1, class T2>
struct is_contiguous<std::pair<T1, T2>>
:
    std::bool_constant<is_contiguous_v<T1> && is_contiguous_v<T2>>
{};

//- Is std::pair composed of labels
template<class T1, class T2>
struct is_contiguous_label<std::pair<T1, T2>>
:
    std::bool_constant
    <
        is_contiguous_label<T1>::value
     && is_contiguous_label<T2>::value
    >
{};

//- Is std::pair composed of scalars
template<class T1, class T2>
struct is_contiguous_scalar<std::pair<T1, T2>>
:
    std::bool_constant
    <
        is_contiguous_scalar<T1>::value
     && is_contiguous_scalar<T2>::value
    >
{};

//- If std::pair has common base element type
template<class T1, class T2>
struct contiguous_element<std::pair<T1, T2>>
{
    using type = std::conditional_t
    <
        std::is_same_v<contiguous_element_t<T1>, contiguous_element_t<T2>>,
        contiguous_element_t<T1>,
        std::pair<T1, T2>  // Fallback: pass-through
    >;
};

//- Hashing for std::pair data
template<class T1, class T2>
struct Hash<std::pair<T1, T2>>
{
    unsigned operator()(const std::pair<T1, T2>& obj, unsigned seed=0) const
    {
        return Hash<T2>()(obj.second, Hash<T1>()(obj.first, seed));
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

namespace Foam
{

//- std::tuple is contiguous if all of its components are contiguous
template<class... Ts>
struct is_contiguous<std::tuple<Ts...>>
:
    std::bool_constant<(is_contiguous_v<Ts> && ...)>
{};

//- Is std::tuple composed of labels
template<class... Ts>
struct is_contiguous_label<std::tuple<Ts...>>
:
    std::bool_constant<(is_contiguous_label<Ts>::value && ...)>
{};

//- Is std::tuple composed of scalars
template<class... Ts>
struct is_contiguous_scalar<std::tuple<Ts...>>
:
    std::bool_constant<(is_contiguous_scalar<Ts>::value && ...)>
{};

//- If std::tuple has common base element type
template<class... Ts>
struct contiguous_element<std::tuple<Ts...>>
{
    using common_type = std::tuple_element_t<0, std::tuple<contiguous_element_t<Ts>...>>;

    using type = std::conditional_t
    <
        (std::is_same_v<common_type, contiguous_element_t<Ts>> && ...),
        common_type,
        std::tuple<Ts...>  // Fallback: pass-through
    >;
};

//- Hashing for std::tuple data
template<class... Ts>
struct Hash<std::tuple<Ts...>>
{
    unsigned operator()(const std::tuple<Ts...>& obj, unsigned seed=0) const
    {
        std::apply
        (
            [&seed](const auto&... x)
            {
                ((seed = Hash<std::decay_t<decltype(x)>>()(x, seed)), ...);
            },
            obj
        );
        return seed;
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read Foam::Pair components from Istream
template<class T>
inline Foam::Istream& Foam::operator>>(Istream& is, Pair<T>& t)
{
    is.readBegin("Pair");
    is >> t.first() >> t.second();
    is.readEnd("Pair");

    is.check(FUNCTION_NAME);
    return is;
}

//- Write Foam::Pair components to Ostream
template<class T>
inline Foam::Ostream& Foam::operator<<(Ostream& os, const Pair<T>& t)
{
    os  << token::BEGIN_LIST
        << t.first() << token::SPACE << t.second()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

namespace Foam
{

//- Read std::pair from Istream
template<class T1, class T2>
inline Istream& operator>>(Istream& is, std::pair<T1,T2>& t)
{
    is.readBegin("pair");
    is >> t.first >> t.second;
    is.readEnd("pair");

    is.check(FUNCTION_NAME);
    return is;
}

//- Write std::pair to Ostream
template<class T1, class T2>
inline Ostream& operator<<
(
    Ostream& os,
    const std::pair<T1,T2>& t
)
{
    os  << token::BEGIN_LIST
        << t.first << token::SPACE << t.second
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

namespace Foam
{

//- Read std::tuple from Istream
template<class... Ts>
inline Istream& operator>>(Istream& is, std::tuple<Ts...>& tup)
{
    is.readBegin("tuple");
    std::apply
    (
        [&is](auto&... x)
        {
            ((is >> x), ...);
        },
        tup
    );
    is.readEnd("tuple");

    is.check(FUNCTION_NAME);
    return is;
}

//- Write std::tuple to Ostream
template<class... Ts>
inline Ostream& operator<<
(
    Ostream& os,
    const std::tuple<Ts...>& tup
)
{
    os << token::BEGIN_LIST;

    std::apply
    (
        [&os](const auto&... x)
        {
            char separator(0);  // Separator *and* tracking 1st invocation
            (
                (
                    (separator ? void(os << separator) : void(0)),
                    (os << x),
                    (separator = token::SPACE)
                ), ...
            );
        },
        tup
    );

    os << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}

} // End namespace Foam


// ************************************************************************* //
