/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Test-CSRAddressing1

Description
    Basic test for CSRAddressing

\*---------------------------------------------------------------------------*/

#include "CSRAddressing.H"
#include "argList.H"
#include "error.H"
#include "IndirectList.H"
#include "IOstreams.H"
#include "Random.H"

using namespace Foam;

static label nFailed = 0;

#define CHECK(cond)                                                           \
    do                                                                        \
    {                                                                         \
        if (!(cond))                                                          \
        {                                                                     \
            ++nFailed;                                                        \
            Info<< "FAILED (line " << __LINE__ << "): " #cond << endl;        \
        }                                                                     \
    } while (false)


template<class IntType>
void printInfo(const CSRAddressing<IntType>& addr)
{
    Info<< "CSR<" << error::demangle<IntType>()
        << "> empty:" << addr.empty()
        << " single:" << addr.single()
        << " length:" << addr.length()
        << " total:" << addr.totalSize()
        << " can_narrow:"
        << pTraits<int32_t>::can_narrow(addr.totalSize()) << nl
        << "    offsets: ";
    addr.offsets().writeList(Info) << nl;
}


template<class IntType>
void testAddressing(const word& name)
{
    typedef CSRAddressing<IntType> Addr;

    Info<< "\n# " << name << nl;

    // Default: empty, no offsets
    {
        Addr addr;
        CHECK(addr.empty());
        CHECK(addr.size() == 0);
        CHECK(addr.offsets().empty());
        CHECK(addr.totalSize() == 0);
        CHECK(addr.findRow(0) == -1);
        CHECK(addr.range(0).empty());
    }

    // From row sizes (including an empty row)
    const labelList rowSizes({3, 0, 2, 5});
    Addr addr(rowSizes);

    CHECK(!addr.empty());
    CHECK(addr.size() == 4);
    CHECK(addr.length() == 4);
    CHECK(!addr.single());
    CHECK(addr.offsets().size() == 5);
    CHECK(addr.begin_value() == 0);
    CHECK(addr.end_value() == 10);
    CHECK(addr.totalSize() == 10);
    CHECK(addr.maxSize() == 5);
    CHECK(pTraits<int32_t>::can_narrow(addr.totalSize()));
    CHECK(addr.sizes() == rowSizes);
    CHECK(addr.localSizes() == rowSizes);
    CHECK(addr.localStarts().size() == 4);

    {
        const List<IntType> expected({0, 3, 3, 5, 10});
        CHECK(addr.offsets() == expected);
    }

    // Row accessors
    CHECK(addr.localStart(2) == 3);
    CHECK(addr.localEnd(2) == 5);
    CHECK(addr.localSize(1) == 0);
    CHECK(addr.localSize(3) == 5);
    CHECK(addr.range(3).start() == 5);
    CHECK(addr.range(3).size() == 5);
    CHECK(addr.range(-1).empty());
    CHECK(addr.range(4).empty());

    // Containment
    CHECK(addr.contains(0));
    CHECK(addr.contains(9));
    CHECK(!addr.contains(10));
    CHECK(!addr.contains(-1));
    CHECK(addr.contains(0, 2));
    CHECK(!addr.contains(0, 3));
    CHECK(!addr.contains(1, 3));  // empty row
    CHECK(addr.contains(3, 9));
    CHECK(!addr.contains(4, 9));  // invalid row

    // Local/global conversion
    CHECK(addr.toGlobal(2, 1) == 4);
    CHECK(addr.toLocal(2, 4) == 1);
    CHECK(addr.toLocal(3, 5) == 0);

    // Row searching: every index must map back to a row containing it
    for (IntType i = 0; i < addr.totalSize(); ++i)
    {
        const label row = addr.findRow(i);
        CHECK(row >= 0 && addr.contains(row, i));

        for (label hint = 0; hint < addr.size(); ++hint)
        {
            CHECK(addr.findRow(hint, i) == row);
        }
    }
    CHECK(addr.findRow(3) == 2);  // skips the empty row
    CHECK(addr.findRow(10) == -1);
    CHECK(addr.findRow(-1) == -1);

    CHECK(addr.findRowAbove(0, 4) == 2);
    CHECK(addr.findRowAbove(0, 2) == -1);
    CHECK(addr.findRowAbove(3, 9) == -1);
    CHECK(addr.findRowBelow(3, 4) == 2);
    CHECK(addr.findRowBelow(3, 7) == -1);
    CHECK(addr.findRowBelow(0, 0) == -1);

    // Ranges
    {
        const List<IntRange<IntType>> ranges = addr.ranges();
        CHECK(ranges.size() == 4);
        CHECK(ranges[2].start() == 3 && ranges[2].size() == 2);
        CHECK(ranges[1].empty());
    }

    // Comparison and copy
    {
        Addr copy(addr);
        CHECK(copy == addr);
        CHECK(!(copy != addr));
        copy.setLocalSize(1, 4);
        CHECK(copy != addr);
        CHECK(copy.localSize(1) == 4);
        CHECK(copy.totalSize() == 14);
    }

    // Reset / clear / resize
    {
        Addr work(addr);
        work.reset(labelList({1, 1}));
        CHECK(work.size() == 2 && work.totalSize() == 2);

        work.clear();
        CHECK(work.empty() && work.offsets().empty());

        work.resize(3);
        CHECK(work.size() == 3 && work.totalSize() == 0);
        CHECK(work.empty());

        work.resize(0);
        CHECK(work.offsets().empty());
    }

    // single / sentinel
    {
        Addr one;
        one.single(7);
        CHECK(one.single());
        CHECK(one.size() == 1 && one.totalSize() == 7);

        Info<< "single "; printInfo(one);

        one.single(2, 3, 10);
        CHECK(one.offsets().size() == 3);
        CHECK(one.localStart(0) == 2 && one.localEnd(0) == 5);
        CHECK(one.totalSize() == 10);

        Info<< "offset-tuple "; printInfo(one);

        Addr sent;
        sent.sentinel(6);
        CHECK(sent.totalSize() == 6);
        CHECK(sent.size() == 0);
        CHECK(sent.offsets().size() == 1);

        Info<< "sentinel "; printInfo(sent);
    }

    // Direct construction: a single offset is not valid (treated as empty)
    {
        Addr direct(Foam::direct_construct, List<IntType>({0, 2, 4}));
        CHECK(direct.size() == 2 && direct.totalSize() == 4);

        Addr bad(Foam::direct_construct, List<IntType>({5}));
        CHECK(bad.offsets().empty());
    }

    // Packing from a list of lists
    {
        List<labelList> lists(3);
        lists[0] = labelList({1, 2});
        lists[2] = labelList({3, 4, 5});

        Addr packed = Addr::pack(lists, true);
        const List<IntType> expected({0, 2, 2, 5});
        CHECK(packed.offsets() == expected);

        List<labelList> none;
        CHECK(Addr::pack(none).offsets().empty());
    }

    // Indirect lists: constructor from row sizes and pack from lists
    {
        const labelList allSizes({9, 3, 0, 2, 7, 5});
        const labelList select({3, 1, 5, 2});  // sizes 2, 3, 5, 0

        IndirectList<label> indirectSizes(allSizes, select);

        Addr fromIndirect(indirectSizes);
        const List<IntType> expected({0, 2, 5, 10, 10});
        CHECK(fromIndirect.offsets() == expected);
        CHECK(fromIndirect.size() == 4);
        CHECK(fromIndirect.totalSize() == 10);

        Addr checked(indirectSizes, true);
        CHECK(checked == fromIndirect);

        Addr fromEmpty
        (
            IndirectList<label>(allSizes, labelList())
        );
        CHECK(fromEmpty.offsets().empty());

        List<labelList> lists(6);
        lists[0] = labelList({1});
        lists[1] = labelList({1, 2, 3});
        lists[3] = labelList({4, 5});
        lists[5] = labelList({6, 7, 8, 9, 10});

        IndirectList<labelList> indirectLists(lists, select);
        Addr packed = Addr::pack(indirectLists, true);
        CHECK(packed.offsets() == expected);

        IndirectList<labelList> noLists(lists, labelList());
        CHECK(Addr::pack(noLists).offsets().empty());
    }

    // null() refers to a nullObject
    {
        const Addr& nul = Addr::null();
        CHECK(isNull(nul));
        CHECK(&nul == &Addr::null());
        CHECK(nul.empty());

        Info<< "null ";
        printInfo(nul);
    }

    // Negative row sizes are ignored; row count is retained
    {
        Addr neg(labelList({2, -1, 3}));
        CHECK(neg.size() == 3);
        CHECK(neg.localSize(1) == 0);
        CHECK(neg.totalSize() == 5);
    }

    // IO round trip
    {
        OStringStream os;
        os << addr;
        IStringStream is(os.str());
        Addr in(is);
        CHECK(in == addr);
    }

    // Random consistency check
    {
        Random rnd(123456);
        labelList sizes(50);
        label total = 0;
        forAll(sizes, i)
        {
            sizes[i] = rnd.position<label>(0, 6);
            total += sizes[i];
        }

        Addr rand(sizes, true);
        CHECK(rand.size() == sizes.size());
        CHECK(rand.totalSize() == total);
        CHECK(rand.sizes() == sizes);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();

    #include "setRootCase.H"

    testAddressing<int32_t>("int32");
    testAddressing<int64_t>("int64");

    // Mixed: 64-bit offsets beyond 32-bit range
    {
        CSRAddressing<int64_t> big;
        big.single(int64_t(INT32_MAX) + 10);
        CHECK(!pTraits<int32_t>::can_narrow(big.totalSize()));
        CHECK(big.findRow(int64_t(INT32_MAX) + 5) == 0);
    }

    if (nFailed)
    {
        Info<< "\nFAILED: " << nFailed << " check(s)" << endl;
        return 1;
    }

    Info<< "\nAll checks passed\n\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
