/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Application
    Test-Tuple2

Description
    Test construction, comparison etc for Tuple2 and Pair.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelPair.H"
#include "Tuple2.H"
#include "label.H"
#include "scalar.H"
#include "List.H"
#include "ListOps.H"
#include "Pstream.H"
#include "nearestFaceAMI.H"
#include "StringStream.H"

#include <functional>
#include <tuple>

namespace Foam
{
namespace Detail
{

//- Compile-time detection if nullary callable can be constexpr
template<auto Fn, int = (Fn(), 0)>
constexpr bool is_constexpr(int) noexcept { return true; }

template<auto Fn>
constexpr bool is_constexpr(...) noexcept { return false; }

} // End namespace Detail

//- Compile-time detection if nullary callable can be constexpr
template<auto Fn>
inline constexpr bool is_constexpr_v = Detail::is_constexpr<Fn>(0);

} // End namespace Foam


using namespace Foam;

// Test for special comparison operation (using std::tie)
// Normal sort on first member, reverse sort on second member
struct special1
{
    template<class T1, class T2>
    bool operator()(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b) const
    {
        return
        (
            std::tie(a.first(), b.second())
          < std::tie(b.first(), a.second())
        );
    }
};


// Test for special comparison operation (using std::tie)
// Reverse sort on first member, normal sort on second member
struct special2
{
    template<class T1, class T2>
    bool operator()(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b) const
    {
        return
        (
            std::tie(b.first(), a.second())
          < std::tie(a.first(), b.second())
        );
    }
};


// Lexicographically sort with flipped first/second members.
template<class Compare = std::less<>>
struct compare_second
{
    template<class T1, class T2>
    bool operator()(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b) const
    {
        return Compare{}
        (
            std::tie(a.second(), a.first()),
            std::tie(b.second(), b.first())
        );
    }
};


// Print content and info
template<class T1, class T2>
void printTuple2Info()
{
    using tuple_type = Tuple2<T1,T2>;

    Info<< error::demangle<tuple_type>()
        << " first:" << error::demangle<T1>()
        << " second:" << error::demangle<T2>()
        << " contiguous:" << Switch(is_contiguous_v<tuple_type>) << nl;
}


// Print content and info
void printTuple2(const word& f, const word& s)
{
    Info<< '(' << f << ' ' << s << ") @ "
        << name(f.data()) << ','
        << name(s.data()) << nl;
}


// Print content and info
void printTuple2(const Tuple2<word, word>& t)
{
    Info<< "tuple: " << t << " @ "
        << name(t.first().data()) << ','
        << name(t.second().data()) << nl;
}


// Print content and info
void printTuple2(const Pair<word>& t)
{
    Info<< "tuple: " << t << " @ "
        << name(t.first().data()) << ','
        << name(t.second().data()) << nl;
}

// Print contiguous information
template<class Type>
void print_contiguous()
{
    Info<< "  " << error::demangle<Type>()
        << " contiguous:" << is_contiguous_v<Type>
        << " label:" << is_contiguous_label<Type>::value
        << " scalar:" << is_contiguous_scalar<Type>::value << endl;
}


// Print is_contiguous/is_contiguous_label/is_contiguous_scalar
// traits for a std::tuple<Types...>
template<class... Types>
void printTupleTraits()
{
    print_contiguous<std::tuple<Types...>>();
}

template<class IndexedScalar>
void testIndexedTuple()
{
    Info<< "Default constructed Tuple: " << IndexedScalar() << nl;

    IndexedScalar t2(3.2, 1);

    Info<< error::demangle<IndexedScalar>()
        << ": " << t2;

    // Simple structured binding usage
    {
        auto& [f, s] = t2;
        f /= 2; s *= 10;
        Info<< " * [0.5, 10]";
    }
    Info<< " -> " << t2 << nl;

    // As list. Generated so that we have duplicate indices
    List<IndexedScalar> list1(3*4);
    for (label i = 0; i < 4; ++i)
    {
        const label j = (i+1);
        const label idx = ((i % 2) ? -1 : 1) * (j);

        list1[i]   = IndexedScalar(idx, (j*j));
        list1[i+4] = IndexedScalar(idx, 2*j);    // duplicate index
        list1[i+8] = IndexedScalar(idx+12, 2*j); // duplicate value
    }

    Info<< "Unsorted tuples:" << nl << list1 << nl;

    // Test min/max functors
    {
        IndexedScalar minIndexed(GREAT, 0);
        IndexedScalar maxIndexed(-GREAT, 0);

        for (const auto& item : list1)
        {
            minFirstEqOp{}(minIndexed, item);
            maxFirstEqOp{}(maxIndexed, item);
        }

        Pstream::combineReduce(minIndexed, minFirstEqOp{});
        Pstream::combineReduce(maxIndexed, maxFirstEqOp{});

        Info<< "Min indexed: " << minIndexed << nl
            << "Max indexed: " << maxIndexed << nl;
    }

    // Test minFirst, maxFirst functors
    {
        IndexedScalar minIndexed(GREAT, 0);
        IndexedScalar maxIndexed(-GREAT, 0);

        for (const auto& item : list1)
        {
            minIndexed = minFirstOp{}(minIndexed, item);
            maxIndexed = maxFirstOp{}(maxIndexed, item);
        }

        Pstream::combineReduce(minIndexed, minFirstEqOp{});
        Pstream::combineReduce(maxIndexed, maxFirstEqOp{});

        Info<< "Min indexed: " << minIndexed << nl
            << "Max indexed: " << maxIndexed << nl;
    }

    Foam::sort(list1, std::less<>());

    Info<< "sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, std::greater<>());

    Info<< "reverse sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, special1{});

    Info<< "special sorted tuples - "
        << "sort on value, reverse on index:"
        << nl << list1 << nl;

    Foam::sort(list1, special2{});

    Info<< "special sorted tuples - "
        << "reverse sort on value, normal on index:"
        << nl << list1 << nl;

    Foam::sort(list1, compare_second<std::less<>>{});

    Info<< "special sorted tuples - "
        << "sorted (index value) instead of (value index):"
        << nl << list1 << nl;

    Foam::sort(list1, compare_second<std::greater<>>{});

    Info<< "special sorted tuples - "
        << "reverse sorted (index value) instead of (value index):"
        << nl << list1 << nl;


    {
        Info<< nl << nl << "Foam::Pair" << nl;

        typedef Pair<label> indexedLabel;

        indexedLabel pr(1, 3);

        Info<< "pair: "
            << pr << " => " << pr.first() << ',' << pr.second();;

        // Simple structured binding usage
        {
            auto& [f, s] = pr;
            f *= 5; s *= 7;
            Info<< " * [5, 7]";
        }
        Info<< " -> " << pr << nl;

        List<indexedLabel> list2 = ListOps::create<indexedLabel>
        (
            list1,
            [](const auto& t2) -> indexedLabel
            {
                return { t2.first(), t2.second() };
            }
        );

        Info<< "Unsorted pairs:" << nl << list2 << nl;
    }


    {
        Info<< nl << nl << "std::pair" << nl;

        typedef std::pair<label, label> indexedLabel;

        indexedLabel pr(1, 3);

        Info<< "pair: "
            << pr << " => "
            << pr.first << ' ' << pr.second << nl;

        List<indexedLabel> list2 = ListOps::create<indexedLabel>
        (
            list1,
            [](const auto& t2) -> indexedLabel
            {
                return { t2.first(), t2.second() };
            }
        );

        Info<< "Unsorted pairs:" << nl << list2 << nl;
    }
}


void test_nearestFace()
{
    using tuple_type = nearestFaceAMI::nearestAndDist;

    Info<< "Default constructed nearestAndDist: " << tuple_type() << nl;

    Info<< "nearestAndDist contiguous:"
        << Switch(is_contiguous_v<tuple_type>) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::noBanner();

    #include "setRootCase.H"

    test_nearestFace();

    // constexpr type:
    {
        using Type = Tuple2<scalar, labelPair>;

        constexpr auto far = Type(Foam::sqr(GREAT), {-1, -1});

        Info<< "constexpr " << far;
        print_contiguous<Type>();
    }

    // constexpr type:
    {
        using Type = Tuple2<scalar, Pair<scalar>>;

        constexpr auto far = Type(Foam::sqr(GREAT), {-1, -1});

        Info<< "constexpr " << far;
        print_contiguous<Type>();
    }

    // constexpr type:
    {
        using Type = Tuple2<point, Pair<scalar>>;

        constexpr auto far = Type(point(), {-1, -1});

        Info<< "constexpr " << far;
        print_contiguous<Type>();
    }

    // constexpr type:
    {
        using Type = Tuple2<point, label>;

        constexpr auto near_point = Type(point::uniform(3), -1);

        Info<< "constexpr " << near_point;
        print_contiguous<Type>();
    }

    // constexpr type:
    {
        using Type = Tuple2<scalar, point>;

        constexpr auto distSqrPoint = Type(3.14159, Foam::zero{});

        Info<< "constexpr " << distSqrPoint;
        print_contiguous<Type>();
    }

    // non-constexpr type:
    {
        using Type = Tuple2<scalar, word>;

        // will not compile:
        // constexpr auto distName = Type(3.14159, "test");

        const auto distName = Type(3.14159, "test");

        Info<< "const " << distName;
        print_contiguous<Type>();
    }

    // Sanity check: the trait itself reports true for a type
    // that we know is constexpr-constructible.
    {
        using Type = Tuple2<scalar, labelPair>;

        constexpr auto fn = +[]{ return Type(1, {-1, -1}); };

        static_assert
        (
            is_constexpr_v<fn>,
            "Tuple2<scalar, labelPair> construction should be constexpr"
        );

        Info<< "Tuple2<scalar, labelPair>(...)) : is_constexpr = "
            << Switch::name(is_constexpr_v<fn>) << nl;
    }

    // Confirm (without a hard compile error) that word's non-constexpr
    // construction/copying prevents Tuple2<scalar, word> from being
    // usable in a constant expression.
    {
        using Type = Tuple2<scalar, word>;

        constexpr auto fn = +[]{ return Type(3.14159, "test"); };

        static_assert
        (
            !is_constexpr_v<fn>,
            "Tuple2<scalar, word> construction should not be constexpr"
        );

        Info<< "Tuple2<scalar, word>(...)) : is_constexpr = "
            << Switch::name(is_constexpr_v<fn>) << nl;
    }

    Info<< "Default constructed Pair: "  << Pair<scalar>() << nl;

    {
        typedef Tuple2<scalar, label> Type;
        Info<< "Tuple2<scalar, label>" << nl;
        testIndexedTuple<Type>();
    }

    {
        word word1("hello");
        word word2("word");

        Info<< "create with ";
        printTuple2(word1, word2);

        Tuple2<word> tup(std::move(word2), std::move(word1));

        printTuple2(tup);

        Info<< "input is now ";
        printTuple2(word1, word2);
    }

    {
        word word1("hello");
        word word2("word");

        Info<< "create with ";
        printTuple2(word1, word2);

        Pair<word> tup(std::move(word2), std::move(word1));

        printTuple2(tup);

        Info<< "input is now ";
        printTuple2(word1, word2);
    }

    {
        Info<< nl << nl << "std::tuple - IO and Hash" << nl;

        typedef std::tuple<label, scalar, word> tupleType;

        tupleType src(1, 2.5, "abc");

        OStringStream oss;
        oss << src;

        Info<< "tuple: " << src << " -> " << oss.str() << nl;

        tupleType dst;
        IStringStream iss(oss.str());
        iss() >> dst;

        Info<< "round-trip: " << dst
            << " equal:" << Switch(src == dst) << nl;

        Hash<tupleType> hasher;

        Info<< "hash(src):" << hasher(src)
            << " hash(dst):" << hasher(dst)
            << " equal:" << Switch(hasher(src) == hasher(dst)) << nl;
    }

    {
        Info<< nl << nl << "std::tuple - contiguous traits" << nl;

        printTupleTraits<scalar, label, label>();
        printTupleTraits<scalar, label, word>();
        printTupleTraits<label, label, label>();
        printTupleTraits<scalar, scalar, scalar>();
    }

    {
        Info<< nl << nl
            << "std::tuple<scalar, label, label> - List IO" << nl;

        typedef std::tuple<scalar, label, label> tupleType;

        List<tupleType> list1
        {
            { 3.1, 2, 5 },
            { 1.5, 4, 1 },
            { 1.5, 4, 1 },
            { 2.2, 1, 9 }
        };

        Info<< "list: " << list1
            << " uniform:" << Switch(list1.uniform()) << nl;

        OStringStream oss;
        oss << list1;

        List<tupleType> list2;
        IStringStream iss(oss.str());
        iss() >> list2;

        Info<< "round-trip: " << list2
            << " equal:" << Switch(list1 == list2) << nl;

        Foam::sort(list2);

        Info<< "sorted (operator<): " << list2
            << " uniform:" << Switch(list2.uniform()) << nl;

        List<tupleType> list3(4, list1[1]);

        Info<< "repeated: " << list3
            << " uniform:" << Switch(list3.uniform()) << nl;
    }

    {
        Info<< nl << nl
            << "std::tuple<scalar, label, word> - List IO" << nl;

        typedef std::tuple<scalar, label, word> tupleType;

        List<tupleType> list1
        {
            { 3.1, 2, "charlie" },
            { 1.5, 4, "alpha" },
            { 1.5, 4, "alpha" },
            { 2.2, 1, "bravo" }
        };

        Info<< "list: " << list1
            << " uniform:" << Switch(list1.uniform()) << nl;

        OStringStream oss;
        oss << list1;

        List<tupleType> list2;
        IStringStream iss(oss.str());
        iss() >> list2;

        Info<< "round-trip: " << list2
            << " equal:" << Switch(list1 == list2) << nl;

        Foam::sort(list2);

        Info<< "sorted (operator<): " << list2
            << " uniform:" << Switch(list2.uniform()) << nl;

        List<tupleType> list3(4, list1[1]);

        Info<< "repeated: " << list3
            << " uniform:" << Switch(list3.uniform()) << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
