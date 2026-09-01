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

#include <functional>
#include <tuple>

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

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
