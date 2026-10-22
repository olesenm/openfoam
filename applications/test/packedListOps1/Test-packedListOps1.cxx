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
    Test-packedListOps1

Description
    PackedList<Width> counterpart of Test-bitSetOps1: exercises the
    range-based operations that were hoisted from bitSet up into the
    generic PackedList<Width> base (labelRange constructor, extract(),
    assign()/assign_at(), set(labelRange), unset(labelRange)) across
    several Width instantiations - including non-word-aligned widths
    (3, 5, 7) where blocks contain unused padding bits and the block-wise
    fast paths must fall back to element-wise handling.

\*---------------------------------------------------------------------------*/

#include "PackedList.H"
#include "labelRange.H"
#include "Random.H"

using namespace Foam;

label nFailTotal = 0;

template<unsigned Width>
inline unsigned int maxVal()
{
    return PackedList<Width>::max_value;
}


// Brute-force reference: bit/element-by-element copy (the old
// implementation, prior to the block-wise rewrite).
template<unsigned Width>
PackedList<Width> naiveExtract
(
    const PackedList<Width>& src,
    const labelRange& range
)
{
    const labelRange clamped(range.subset0(src.size()));
    PackedList<Width> out(clamped.size());
    label i = 0;
    for (label pos = clamped.start(); i < clamped.size(); ++pos, ++i)
    {
        out.set(i, src.get(pos));
    }
    return out;
}


// Brute-force reference for assign(offset, source) - mirrors the
// documented contract: always grows to (offset + source.size()) when
// source is non-empty, regardless of source content.
template<unsigned Width>
PackedList<Width> naiveAssign
(
    PackedList<Width> dest,
    const label offset,
    const PackedList<Width>& source
)
{
    if (source.empty())
    {
        return dest;
    }

    const label lower = Foam::max(0, offset);
    const label upper = (offset + source.size());

    if (upper <= lower)
    {
        return dest;
    }

    if (dest.size() < upper)
    {
        dest.resize(upper);
    }

    for (label i = 0; i < source.size(); ++i)
    {
        dest.set(i + offset, source.get(i));
    }
    return dest;
}


// Brute-force reference for set(range, val)/unset(range):
// val == 0 exercises the unset() contract (never grows beyond the
// original size); val != 0 exercises the set() contract (auto-vivifies).
template<unsigned Width>
PackedList<Width> naiveSetRange
(
    PackedList<Width> dest,
    const labelRange& range,
    unsigned int val
)
{
    if (!val)
    {
        const labelRange slice(range.subset0(dest.size()));
        for (label i = slice.start(); i < slice.end_value(); ++i)
        {
            dest.set(i, 0u);
        }
        return dest;
    }

    labelRange slice(range);
    slice.adjust();
    if (slice.empty())
    {
        return dest;
    }
    if (slice.end_value() >= dest.size())
    {
        dest.resize(slice.end_value(), 0u);
    }
    for (label i = slice.start(); i < slice.end_value(); ++i)
    {
        dest.set(i, val);
    }
    return dest;
}


template<unsigned Width>
PackedList<Width> makePattern(label n, label mod, label thresh, unsigned int val)
{
    PackedList<Width> out(n);
    for (label i = 0; i < n; ++i)
    {
        out.set(i, (((i*5 + 2) % mod) < thresh) ? val : 0u);
    }
    return out;
}


template<unsigned Width>
void test_extraction()
{
    Info<< "\n[Width=" << Width << "] range-extraction constructor / extract()\n";

    const unsigned int val = maxVal<Width>();

    PackedList<Width> src(200);
    for (label i = 0; i < src.size(); ++i)
    {
        src.set(i, (((i*7 + 3) % 11) < 5) ? val : 0u);
    }

    label nTest = 0;
    label nFail = 0;

    for (label start : {0, 1, 5, 31, 32, 33, 63, 64, 65, 100, 199, 200})
    {
        for (label len : {0, 1, 2, 5, 31, 32, 33, 64, 65, 150, 500})
        {
            ++nTest;
            const labelRange range(start, len);

            const PackedList<Width> expected(naiveExtract(src, range));
            const PackedList<Width> got1(src, range);          // constructor
            const PackedList<Width> got2(src.extract(range));       // extract(range)
            const PackedList<Width> got3(src.extract(start, len));  // extract(pos,len)

            const bool ok =
            (
                expected.size() == got1.size() && expected == got1
             && got1.size() == got2.size() && got1 == got2
             && got1.size() == got3.size() && got1 == got3
            );

            if (!ok)
            {
                ++nFail;
                Info<< "fail: start=" << start << " len=" << len << nl;
            }
        }
    }

    Info<< "extraction tests: " << nTest << " total, "
        << nFail << " failed" << nl;

    nFailTotal += nFail;

    // Corner cases
    {
        const PackedList<Width> none;
        const PackedList<Width> e1(none.extract(0, 10));
        Info<< "extract from empty source: size=" << e1.size()
            << " (expect 0)" << nl;
        if (e1.size() != 0) ++nFailTotal;
    }
    {
        const PackedList<Width> full(src.extract(labelRange(0, src.size())));
        const bool ok = (full == src);
        Info<< "extract full range equals source: "
            << (ok ? "pass" : "fail") << nl;
        if (!ok) ++nFailTotal;
    }
}


template<unsigned Width>
void test_assign_offset()
{
    Info<< "\n[Width=" << Width << "] assign(offset, source)\n";

    const unsigned int val = maxVal<Width>();

    PackedList<Width> base(150);
    for (label i = 0; i < base.size(); ++i)
    {
        base.set(i, (((i*3 + 1) % 7) < 4) ? val : 0u);
    }

    label nTest = 0;
    label nFail = 0;

    for
    (
        label offset
      : {-40, -32, -5, -1, 0, 1, 5, 31, 32, 33, 64, 90, 120, 149, 200}
    )
    {
        const PackedList<Width> sources[] =
        {
            PackedList<Width>(),                        // empty
            PackedList<Width>(20, val),                  // all-on
            PackedList<Width>(20),                       // all-off
            makePattern<Width>(20, 7, 4, val),            // sparse, small
            makePattern<Width>(70, 11, 6, val),           // crosses many blocks
            makePattern<Width>(300, 3, 2, val)            // extends well beyond base
        };

        for (const PackedList<Width>& source : sources)
        {
            ++nTest;

            const PackedList<Width> expected(naiveAssign(base, offset, source));

            PackedList<Width> got(base);
            got.assign(offset, source);

            const bool ok = (expected == got && expected.size() == got.size());

            if (!ok)
            {
                ++nFail;
                Info<< "fail: offset=" << offset
                    << " source.size=" << source.size() << nl
                    << "  expected size=" << expected.size()
                    << " got size=" << got.size() << nl;
            }
        }
    }

    Info<< "assign(offset,source) tests: " << nTest << " total, "
        << nFail << " failed" << nl;

    nFailTotal += nFail;

    // Randomised fuzz test
    {
        Random rndGen(20260908);

        label nFuzz = 0;
        label nFuzzFail = 0;

        for (label trial = 0; trial < 500; ++trial)
        {
            const label baseSize = rndGen.position<label>(0, 300);
            const label sourceSize = rndGen.position<label>(0, 300);
            const label offset = rndGen.position<label>(-150, 350);

            PackedList<Width> b(baseSize);
            for (label i = 0; i < baseSize; ++i)
            {
                b.set(i, (rndGen.sample01<scalar>() < 0.5) ? val : 0u);
            }

            PackedList<Width> s(sourceSize);
            for (label i = 0; i < sourceSize; ++i)
            {
                s.set(i, (rndGen.sample01<scalar>() < 0.35) ? val : 0u);
            }

            ++nFuzz;

            const PackedList<Width> expected(naiveAssign(b, offset, s));

            PackedList<Width> got(b);
            got.assign(offset, s);

            if (!(expected == got) || expected.size() != got.size())
            {
                ++nFuzzFail;
                Info<< "fuzz fail: trial=" << trial
                    << " baseSize=" << baseSize
                    << " sourceSize=" << sourceSize
                    << " offset=" << offset << nl;
            }
        }

        Info<< "assign(offset,source) fuzz: " << nFuzz << " total, "
            << nFuzzFail << " failed" << nl;

        nFailTotal += nFuzzFail;
    }
}


template<unsigned Width>
void test_set_unset_range()
{
    Info<< "\n[Width=" << Width << "] set(labelRange,val) / unset(labelRange)\n";

    const unsigned int val = maxVal<Width>();

    PackedList<Width> base(150);
    for (label i = 0; i < base.size(); ++i)
    {
        base.set(i, (((i*3 + 1) % 7) < 4) ? val : 0u);
    }

    label nTest = 0;
    label nFail = 0;

    for (label start : {-40, -5, -1, 0, 1, 5, 31, 32, 33, 64, 90, 149, 200})
    {
        for (label len : {-1, 0, 1, 2, 5, 31, 32, 64, 150, 400})
        {
            ++nTest;
            const labelRange range(start, len);

            // set()
            {
                const PackedList<Width> expected
                (
                    naiveSetRange(base, range, val)
                );

                PackedList<Width> got(base);
                got.set(range, val);

                const bool ok =
                    (expected == got && expected.size() == got.size());

                if (!ok)
                {
                    ++nFail;
                    Info<< "set fail: start=" << start << " len=" << len
                        << " expected size=" << expected.size()
                        << " got size=" << got.size() << nl;
                }
            }

            // unset()
            {
                const PackedList<Width> expected
                (
                    naiveSetRange(base, range, 0u)
                );

                PackedList<Width> got(base);
                got.unset(range);

                const bool ok =
                    (expected == got && expected.size() == got.size());

                if (!ok)
                {
                    ++nFail;
                    Info<< "unset fail: start=" << start << " len=" << len
                        << " expected size=" << expected.size()
                        << " got size=" << got.size() << nl;
                }
            }
        }
    }

    Info<< "set/unset(labelRange) tests: " << nTest << " total, "
        << nFail << " failed" << nl;

    nFailTotal += nFail;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    // Widths 1,2,4,8 are word-aligned (exercise the block-wise fast
    // paths); widths 3,5,7 are not (exercise the element-wise fallback)
    test_extraction<1>();
    test_extraction<2>();
    test_extraction<3>();
    test_extraction<4>();
    test_extraction<5>();
    test_extraction<7>();
    test_extraction<8>();

    test_assign_offset<1>();
    test_assign_offset<2>();
    test_assign_offset<3>();
    test_assign_offset<4>();
    test_assign_offset<5>();
    test_assign_offset<7>();
    test_assign_offset<8>();

    test_set_unset_range<1>();
    test_set_unset_range<2>();
    test_set_unset_range<3>();
    test_set_unset_range<4>();
    test_set_unset_range<5>();
    test_set_unset_range<7>();
    test_set_unset_range<8>();

    Info<< nl;
    if (nFailTotal)
    {
        Info<< nFailTotal << " FAILURES" << nl;
        return 1;
    }

    Info<< "All tests passed" << nl;
    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
