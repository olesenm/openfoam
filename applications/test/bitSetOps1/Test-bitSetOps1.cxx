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
    Test-bitSetOps1

Description
    More bitSet operations

\*---------------------------------------------------------------------------*/

#include "bitSet.H"
#include "FlatOutput.H"
#include "Random.H"

using namespace Foam;

void test_extraction()
{
    Info<< "\nrange-extraction constructor / extract()\n";

    // Brute-force reference: bit-by-bit copy (the old implementation)
    auto naiveExtract = [](const bitSet& src, const labelRange& range)
    {
        const labelRange clamped(range.subset0(src.size()));
        bitSet out(clamped.size());
        label i = 0;
        for (label pos = clamped.start(); i < clamped.size(); ++pos, ++i)
        {
            out.set(i, src.get(pos));
        }
        return out;
    };

    // Build a source pattern that is not simply periodic so that
    // shifted (mis-aligned) block copies are properly exercised.
    bitSet src(200);
    forAll(src, i)
    {
        src[i] = ((i*7 + 3) % 11) < 5;
    }

    label nFail = 0;
    label nTest = 0;

    // Exercise block-aligned and non-aligned starts, short and long
    // lengths, as well as ranges that overrun the end of the source.
    for (label start : {0, 1, 5, 31, 32, 33, 63, 64, 65, 100, 199, 200})
    {
        for (label len : {0, 1, 2, 5, 31, 32, 33, 64, 65, 150, 500})
        {
            ++nTest;
            const labelRange range(start, len);

            const bitSet expected(naiveExtract(src, range));
            const bitSet got1(src, range);       // constructor
            const bitSet got2(src.extract(range));  // extract(range)
            const bitSet got3(src.extract(start, len)); // extract(pos,len)

            const bool ok =
            (
                expected.size() == got1.size()
             && expected == got1
             && got1.size() == got2.size() && got1 == got2
             && got1.size() == got3.size() && got1 == got3
            );

            if (!ok)
            {
                ++nFail;
                Info<< "fail: start=" << start << " len=" << len << nl
                    << "  expected: " << flatOutput(expected) << nl
                    << "  ctor    : " << flatOutput(got1) << nl
                    << "  extract : " << flatOutput(got2) << nl
                    << "  extract2: " << flatOutput(got3) << nl;
            }
        }
    }

    Info<< "extraction tests: " << nTest << " total, "
        << nFail << " failed" << nl;

    // Simple corner cases
    {
        const bitSet none;
        const bitSet e1(none.extract(0, 10));
        Info<< "extract from empty source: size=" << e1.size()
            << " (expect 0)" << nl;
    }
    {
        const bitSet full(src.extract(labelRange(0, src.size())));
        Info<< "extract full range equals source: "
            << (full == src ? "pass" : "fail") << nl;
    }
}


void test_unset_offset()
{
    Info<< "\nunset(offset, other)\n";

    // Regression: 'this' is small (1 block) but 'other' extends far
    // beyond this->size(), while this->find_last() is still small.
    // With only "lower != offset" as the fallback trigger, the
    // block-wise path is taken even though offset+thatExtent
    // exceeds this->size() - out-of-bounds write into blocks_.
    {
        bitSet small(10);
        small.set(3);   // find_last() == 3  -> thisExtent == 4

        bitSet far(100);
        far.set(90);    // find_last() == 90 -> thatExtent == 91

        Info<< "regression (oob): this.size()=" << small.size()
            << " other.size()=" << far.size() << " ... ";

        small.unset(0, far);   // <-- potential OOB write here

        Info<< "survived, this=" << small << nl;
    }

    // Brute-force reference: unset every 'on' bit of 'other', shifted
    // by 'offset', bit-by-bit (silently ignoring out-of-range).
    auto naiveUnset = [](bitSet base, const label offset, const bitSet& other)
    {
        for (label i : other)
        {
            base.unset(i + offset);
        }
        return base;
    };

    // Non-trivial patterns so that both 0 and 1 bits, and block
    // (mis-)alignment are properly exercised.
    bitSet base(150);
    forAll(base, i)
    {
        base[i] = ((i*3 + 1) % 7) < 4;
    }

    auto makeOther = [](label n, label mod, label thresh)
    {
        bitSet b(n);
        forAll(b, i)
        {
            b[i] = ((i*5 + 2) % mod) < thresh;
        }
        return b;
    };

    label nTest = 0;
    label nFail = 0;

    for
    (
        label offset
      : {-40, -32, -5, -1, 0, 1, 5, 31, 32, 33, 64, 90, 120, 149, 200}
    )
    {
        for
        (
            const bitSet& other
          : {
                bitSet(),                    // empty
                bitSet(20, true),             // all-on
                bitSet(20),                   // all-off (none)
                makeOther(20, 7, 4),          // sparse, small
                makeOther(70, 11, 6),         // sparse, crosses many blocks
                makeOther(300, 3, 2)          // extends well beyond base
            }
        )
        {
            ++nTest;

            bitSet expected(naiveUnset(base, offset, other));

            bitSet got(base);
            got.unset(offset, other);

            if (!(expected == got) || expected.size() != got.size())
            {
                ++nFail;
                Info<< "fail: offset=" << offset
                    << " other.size=" << other.size() << nl
                    << "  expected: " << flatOutput(expected) << nl
                    << "  got     : " << flatOutput(got) << nl;
            }
        }
    }

    Info<< "unset(offset,other) tests: " << nTest << " total, "
        << nFail << " failed" << nl;

    // Corner cases
    {
        bitSet a(base);
        a.unset(0, bitSet());   // other empty: no-op
        Info<< "unset with empty other is no-op: "
            << (a == base ? "pass" : "fail") << nl;
    }
    {
        bitSet a;                // this->none() (empty)
        a.unset(5, makeOther(20, 7, 4));
        Info<< "unset on empty *this* stays empty: "
            << (a.size() == 0 ? "pass" : "fail") << nl;
    }
    {
        // offset == 0 shortcut path vs bit-by-bit reference
        bitSet a(base);
        const bitSet other(makeOther(70, 11, 6));
        a.unset(0, other);

        bitSet expect(naiveUnset(base, 0, other));
        Info<< "unset(0, other) matches unset(other): "
            << (a == expect ? "pass" : "fail") << nl;
    }

    // Randomised fuzz test: random sizes/densities for both bitsets,
    // and random (incl. negative and overrunning) offsets - to catch
    // any block-alignment edge case not covered by the hand-picked
    // cases above.
    {
        Random rndGen(20260904);

        label nFuzz = 0;
        label nFuzzFail = 0;

        for (label trial = 0; trial < 2000; ++trial)
        {
            const label baseSize = rndGen.position<label>(0, 300);
            const label otherSize = rndGen.position<label>(0, 300);
            const label offset = rndGen.position<label>(-150, 350);

            bitSet b(baseSize);
            forAll(b, i)
            {
                b[i] = (rndGen.sample01<scalar>() < 0.5);
            }

            bitSet o(otherSize);
            forAll(o, i)
            {
                o[i] = (rndGen.sample01<scalar>() < 0.35);
            }

            ++nFuzz;

            bitSet expected(naiveUnset(b, offset, o));

            bitSet got(b);
            got.unset(offset, o);

            if (!(expected == got) || expected.size() != got.size())
            {
                ++nFuzzFail;
                Info<< "fuzz fail: trial=" << trial
                    << " baseSize=" << baseSize
                    << " otherSize=" << otherSize
                    << " offset=" << offset << nl
                    << "  expected: " << flatOutput(expected) << nl
                    << "  got     : " << flatOutput(got) << nl;
            }
        }

        Info<< "unset(offset,other) fuzz: " << nFuzz << " total, "
            << nFuzzFail << " failed" << nl;
    }

    // Targeted fuzz: 'this' has a large declared size() but a genuinely
    // sparse tail (find_last() well below size()-1). This specifically
    // exercises the "upper = min(thisExtent, offset+thatExtent)" bound
    // and the "upper != (offset+thatExtent)" fallback trigger in the
    // 'Uncommon case' branch, which the earlier dense-pattern tests
    // above never really hit (their find_last() was almost always
    // size()-1).
    {
        Random rndGen(9876543);

        label nFuzz = 0;
        label nFuzzFail = 0;

        for (label trial = 0; trial < 3000; ++trial)
        {
            const label baseSize = rndGen.position<label>(1, 300);

            // Only allow 'on' bits within the leading part of 'base',
            // leaving a genuinely all-zero tail of random length.
            const label activeLen = rndGen.position<label>(0, baseSize);

            bitSet b(baseSize);
            for (label i = 0; i < activeLen; ++i)
            {
                b[i] = (rndGen.sample01<scalar>() < 0.5);
            }
            // (positions [activeLen, baseSize) remain 'off')

            const label otherSize = rndGen.position<label>(0, 300);
            bitSet o(otherSize);
            forAll(o, i)
            {
                o[i] = (rndGen.sample01<scalar>() < 0.35);
            }

            const label offset = rndGen.position<label>(-150, 350);

            ++nFuzz;

            bitSet expected(naiveUnset(b, offset, o));

            bitSet got(b);
            got.unset(offset, o);

            if (!(expected == got) || expected.size() != got.size())
            {
                ++nFuzzFail;
                Info<< "sparse-tail fuzz fail: trial=" << trial
                    << " baseSize=" << baseSize
                    << " activeLen=" << activeLen
                    << " find_last=" << b.find_last()
                    << " otherSize=" << otherSize
                    << " offset=" << offset << nl
                    << "  expected: (" << expected.size() << ") "
                    << flatOutput(expected) << nl
                    << "  got     : (" << got.size() << ") "
                    << flatOutput(got) << nl;
            }
        }

        Info<< "unset(offset,other) sparse-tail fuzz: " << nFuzz
            << " total, " << nFuzzFail << " failed" << nl;
    }

    // Deterministic boundary check: force thisExtent == offset+thatExtent
    // exactly (the tightest margin at which the block-wise path is still
    // considered safe), across both block-aligned and mis-aligned
    // 'lower' offsets, and several relative positions of the last bit
    // within its block (to exercise the shift/carry off-by-one boundary
    // specifically).
    {
        label nTest = 0;
        label nFail = 0;

        // otherSize/onLen chosen so thatExtent (find_last()+1) is exact
        for (label otherSize : {1, 5, 32, 33, 63, 64, 70})
        {
            // 'other' has its last 'on' bit exactly at otherSize-1,
            // so thatExtent == otherSize precisely.
            bitSet other(otherSize);
            forAll(other, i)
            {
                other[i] = ((i*5 + 2) % 7) < 4;
            }
            other.set(otherSize - 1);   // guarantee find_last() == otherSize-1

            for (label offset : {0, 1, 5, 31, 32, 33, 63, 64, 100})
            {
                if (offset == 0) continue;  // handled by the offset==0 shortcut

                const label thatExtent = otherSize;
                const label thisExtent = offset + thatExtent;  // exact boundary

                bitSet b(thisExtent + 5);   // size() a bit larger than extent
                forAll(b, i)
                {
                    b[i] = ((i*3 + 1) % 5) < 3;
                }
                b.unset(labelRange(thisExtent, b.size() - thisExtent));
                b.set(thisExtent - 1);   // guarantee find_last() == thisExtent-1

                ++nTest;

                bitSet expected(naiveUnset(b, offset, other));

                bitSet got(b);
                got.unset(offset, other);

                if (!(expected == got) || expected.size() != got.size())
                {
                    ++nFail;
                    Info<< "boundary fail: otherSize=" << otherSize
                        << " offset=" << offset
                        << " thisExtent=" << thisExtent << nl
                        << "  expected: (" << expected.size() << ") "
                        << flatOutput(expected) << nl
                        << "  got     : (" << got.size() << ") "
                        << flatOutput(got) << nl;
                }
            }
        }

        Info<< "unset(offset,other) exact-boundary tests: " << nTest
            << " total, " << nFail << " failed" << nl;
    }
}


void test_assign_offset()
{
    Info<< "\nassign(offset, source)\n";

    // Brute-force reference: overwrite bit-by-bit at 'offset', silently
    // ignoring/clipping the leading part of 'source' that maps to a
    // negative destination index (mirrors bitSet::set() semantics).
    //
    // NOTE: the real assign() explicitly resize()s to the full extent
    // (offset + source.size()) whenever source has any 'on' bit - even
    // if the trailing bits of source happen to be 'off'. A plain
    // bit-by-bit set() loop would only grow lazily on 'on' bits, so we
    // must replicate that explicit sizing policy here too, otherwise
    // the reference under-reports size() whenever source's last bit(s)
    // are 'off'.
    auto naiveAssign = [](bitSet dest, const label offset, const bitSet& source)
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

        // Always grow to the full extent, regardless of source content
        // (matches the documented assign() contract).
        if (dest.size() < upper)
        {
            dest.resize(upper);
        }

        forAll(source, i)
        {
            dest.set(i + offset, source.get(i));
        }
        return dest;
    };

    // Non-trivial destination pattern so both 0 and 1 bits already
    // present get properly overwritten (not just OR-ed).
    bitSet base(150);
    forAll(base, i)
    {
        base[i] = ((i*3 + 1) % 7) < 4;
    }

    auto makeSource = [](label n, label mod, label thresh)
    {
        bitSet b(n);
        forAll(b, i)
        {
            b[i] = ((i*5 + 2) % mod) < thresh;
        }
        return b;
    };

    label nTest = 0;
    label nFail = 0;

    for
    (
        label offset
      : {-40, -32, -5, -1, 0, 1, 5, 31, 32, 33, 64, 90, 120, 149, 200}
    )
    {
        for
        (
            const bitSet& source
          : {
                bitSet(),                      // empty
                bitSet(20, true),               // all-on
                bitSet(20),                     // all-off (none)
                makeSource(20, 7, 4),           // sparse, small
                makeSource(70, 11, 6),          // sparse, crosses many blocks
                makeSource(300, 3, 2)           // extends well beyond base
            }
        )
        {
            ++nTest;

            bitSet expected(naiveAssign(base, offset, source));

            bitSet got(base);
            got.assign(offset, source);

            if (!(expected == got) || expected.size() != got.size())
            {
                ++nFail;
                Info<< "fail: offset=" << offset
                    << " source.size=" << source.size() << nl
                    << "  expected: (" << expected.size() << ") "
                    << flatOutput(expected) << nl
                    << "  got     : (" << got.size() << ") "
                    << flatOutput(got) << nl;
            }
        }
    }

    Info<< "assign(offset,source) tests: " << nTest << " total, "
        << nFail << " failed" << nl;

    // Corner cases
    {
        bitSet a(base);
        a.assign(0, bitSet());   // source empty: no-op
        Info<< "assign with empty source is no-op: "
            << (a == base ? "pass" : "fail") << nl;
    }
    {
        // offset == 0: plain overwrite of the leading part of *this*
        bitSet a(base);
        const bitSet source(makeSource(70, 11, 6));
        a.assign(0, source);

        bitSet expect(naiveAssign(base, 0, source));
        Info<< "assign(0, source) matches bit-by-bit reference: "
            << (a == expect ? "pass" : "fail") << nl;
    }
    {
        // source extends past the end of *this*: must grow
        bitSet a(20);
        a.set(labelRange(0, 20));   // all-on, small
        const bitSet source(makeSource(300, 3, 2));

        bitSet expect(naiveAssign(a, 10, source));

        a.assign(10, source);

        Info<< "assign() growing past current size: "
            << (a == expect ? "pass" : "fail")
            << " size=" << a.size() << " (expect " << expect.size() << ")"
            << nl;
    }
    {
        // source is all-zero but extends past the end of *this*:
        // must still grow (per documented contract), filled with 0s.
        bitSet a(5);   // all-off
        const bitSet source(10);   // all-off, size 10

        a.assign(0, source);

        Info<< "assign() grows for all-zero source: "
            << (a.size() == 10 && a.none() ? "pass" : "fail")
            << " size=" << a.size() << " (expect 10)" << nl;
    }
    {
        // all-zero source, offset > 0, straddling several blocks
        bitSet a(40);
        forAll(a, i) { a[i] = ((i*7 + 3) % 5) < 2; }
        const bitSet before(a);
        const bitSet source(90);   // all-off, size 90

        bitSet expect(naiveAssign(a, 33, source));
        a.assign(33, source);

        Info<< "assign() all-zero source, offset=33, grows+clears: "
            << (a == expect ? "pass" : "fail")
            << " size=" << a.size() << " (expect " << expect.size() << ")"
            << nl;
    }

    // Randomised fuzz test
    {
        Random rndGen(20260904);

        label nFuzz = 0;
        label nFuzzFail = 0;

        for (label trial = 0; trial < 2000; ++trial)
        {
            const label baseSize = rndGen.position<label>(0, 300);
            const label sourceSize = rndGen.position<label>(0, 300);
            const label offset = rndGen.position<label>(-150, 350);

            bitSet b(baseSize);
            forAll(b, i)
            {
                b[i] = (rndGen.sample01<scalar>() < 0.5);
            }

            bitSet s(sourceSize);
            forAll(s, i)
            {
                s[i] = (rndGen.sample01<scalar>() < 0.35);
            }

            ++nFuzz;

            bitSet expected(naiveAssign(b, offset, s));

            bitSet got(b);
            got.assign(offset, s);

            if (!(expected == got) || expected.size() != got.size())
            {
                ++nFuzzFail;
                Info<< "fuzz fail: trial=" << trial
                    << " baseSize=" << baseSize
                    << " sourceSize=" << sourceSize
                    << " offset=" << offset << nl
                    << "  expected: (" << expected.size() << ") "
                    << flatOutput(expected) << nl
                    << "  got     : (" << got.size() << ") "
                    << flatOutput(got) << nl;
            }
        }

        Info<< "assign(offset,source) fuzz: " << nFuzz << " total, "
            << nFuzzFail << " failed" << nl;
    }
}


void test_toc()
{
    Info<< "\ntoc()/sortedToc() with offset\n";

    bitSet small;
    small.set(labelList{0, 2, 5});

    Info<< "toc(): " << flatOutput(small.toc()) << nl
        << "toc(100): " << flatOutput(small.toc(100)) << nl
        << "sortedToc(): " << flatOutput(small.sortedToc()) << nl
        << "sortedToc(-5): " << flatOutput(small.sortedToc(-5)) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    test_extraction();

    test_unset_offset();

    test_assign_offset();

    test_toc();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
