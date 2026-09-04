/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "PackedList.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<unsigned Width>
void Foam::PackedList<Width>::assign_at
(
    const label offset,
    const PackedList<Width>& source
)
{
    // Assumes: offset >= 0 and (offset + source.size()) <= this->size().
    // The caller is responsible for any growing/resizing and for the
    // trivial/negative-offset cases.

    const label len = source.size();

    if (len <= 0)
    {
        // Nothing to do
        return;
    }

    constexpr bool word_aligned = (bits_per_block % Width == 0);

    if constexpr (word_aligned)
    {
        // 1. Clear the target range first (bit-precise, block-wise).
        unset(labelRange(offset, len));

        // 2. Assign by OR-ing in the (possibly shifted) bits from 'source'.
        //
        // Process block-wise as much as possible.
        // Each output block is assembled from (up to) two adjacent
        // source blocks, shifted according to the sub-block element
        // offset of the destination start.

        // Number of input blocks (before trimming)
        const label nblocks = source.num_blocks();

        // The begin block/bit-shift (of the destination)
        const label bblock   = block_index(offset);
        const unsigned bshift = Width * sub_index(offset);

        // The input start
        auto iter = source.blocks_.cbegin();
        const auto last = (iter + nblocks);

        // The output start
        auto output = (blocks_.begin() + bblock);

        if (!bshift)
        {
            // Fast path: block-aligned start - direct OR, no shift required
            while (iter != last)
            {
                *output |= *iter;

                ++iter;
                ++output;
            }
        }
        else
        {
            // General case: each source block spans (up to) two
            // destination blocks.

            const auto endOutput = blocks_.end();
            unsigned int carry = 0;

            while (iter != last)
            {
                const unsigned int val = *iter;
                ++iter;

                // Low bits of val -> upper bits of the current output block
                *output |= (carry | (val << bshift));
                ++output;

                // High bits of val -> lower bits of the next output block
                carry = (val >> (bits_per_block - bshift));
            }

            if (carry && output != endOutput)
            {
                *output |= carry;
            }
        }
    }
    else
    {
        // Fallback: element-wise copy (blocks have unused padding bits)
        for (label i = 0; i < len; ++i)
        {
            setValue_unsafe(offset + i, source.get(i));
        }
    }
}


template<unsigned Width>
void Foam::PackedList<Width>::extract_copy
(
    const PackedList<Width>& source,
    PackedList<Width>& target,
    const labelRange& range
)
{
    if (&source == &target)
    {
        FatalErrorInFunction
            << "Cannot extract a copy into itself !"
            << abort(FatalError);
    }

    // With clamped start/size
    const labelRange clamped(range.subset0(source.size()));
    const label pos = clamped.start();
    const label len = clamped.size();

    if (len <= 0)
    {
        // nothing to do
        target.clear();
        return;
    }

    // Resize with minimal over-allocation
    target.resize_fill(len, 0u);

    // Block-wise copying is only valid when elements pack without any
    // unused padding bits within a block. If there are padding bits,
    // the two-block combine (below) would go awry.

    constexpr bool word_aligned = (bits_per_block % Width == 0);

    if constexpr (word_aligned)
    {
        // Process block-wise as much as possible.
        //
        // Each output block is assembled from (up to) two adjacent
        // source blocks, shifted according to the sub-block element
        // offset of the range start.

        // Number of output blocks (before trimming),
        // also a safe input range due to prior size clamping.
        const label nblocks = target.num_blocks();

        // The begin block/bit-shift
        const label bblock    = block_index(pos);
        const unsigned bshift = Width * sub_index(pos);

        // The input start/end
        auto iter = (source.blocks_.begin() + bblock);
        const auto last = (iter + nblocks);

        auto output = target.blocks_.begin();

        if (!bshift)
        {
            // Fast path: block-aligned start - direct copy, no shift required
            std::copy(iter, last, output);
        }
        else
        {
            // General case: combine two adjacent source blocks.

            const auto endIter = source.blocks_.end();
            while (iter < last)
            {
                unsigned int val = (*iter >> bshift);

                ++iter;

                if (iter < endIter)
                {
                    val |= (*iter << (bits_per_block - bshift));
                }

                *output = val;
                ++output;
            }
        }

        // Mask off any excess bits copied into the final (partial) block
        target.clear_trailing_bits();
    }
    else
    {
        // Fallback: element-wise copy (blocks have unused padding bits)
        label srcPos = pos;

        for (label i = 0; i < len; ++i)
        {
            target.setValue_unsafe(i, source.get(srcPos));
            ++srcPos;
        }
    }
}


template<unsigned Width>
void Foam::PackedList<Width>::extract_copy
(
    const PackedList<Width>& source,
    PackedList<Width>& target,
    const label pos,
    const label len
)
{
    // With start/len clamped to the source extents
    const labelRange range(clamp_pos_len(pos, len, source.size()));

    if (range.empty())
    {
        // Nothing to copy
        target.clear();
    }
    else
    {
        // Extract the range
        extract_copy(source, target, range);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned Width>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const labelUList& addr
)
:
    PackedList<Width>(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        setValue_unsafe(i, list.get(addr[i]));
    }
}


template<unsigned Width>
template<class Addr>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const IndirectListBase<label, Addr>& addr
)
:
    PackedList<Width>(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        setValue_unsafe(i, list.get(addr[i]));
    }
}


template<unsigned Width>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& source,
    const labelRange& range
)
:
    PackedList<Width>()
{
    extract_copy(source, *this, range);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width>
bool Foam::PackedList<Width>::uniform() const noexcept
{
    // Trivial cases
    if (empty())
    {
        return false;
    }
    else if (size() == 1)
    {
        return true;
    }

    // The value of the first element for testing
    const unsigned int val = this->get(0);

    if (!val)
    {
        // No bits set: just check there are no non-zero blocks
        // - like bitSet::none()
        return (-1 == find_first_block());
    }
    else if (val == PackedList<Width>::max_value)
    {
        // All bits set: just check there are no zero blocks
        // - like bitSet::all()
        return (-1 == find_first_not_block());
    }
    else
    {
        bool identical = true;

        const label nblocks = num_blocks(size());

        if (nblocks > 1)
        {
            // Check each complete block (nblocks-1)

            // Fill value for complete blocks
            const unsigned int blockval = repeated_value(val);

            const auto* block_ptr = blocks_.cdata();

            for (label blocki = 0; identical && blocki < (nblocks-1); ++blocki)
            {
                identical = (block_ptr[blocki] == blockval);
            }
        }

        // Partial block: check manually
        for
        (
            label elemi = elem_per_block*(nblocks-1);
            identical && elemi < size();
            ++elemi
        )
        {
            identical = (val == get(elemi));
        }

        return identical;
    }
}


template<unsigned Width>
bool Foam::PackedList<Width>::equal
(
    const PackedList<Width>& other
) const noexcept
{
    if (size() != other.size())
    {
        return false;
    }

    // Like
    // std::equal(cdata(), cdata()+num_blocks(), other.cdata())

    const label nblocks = this->num_blocks();

    const auto* block_ptr = this->blocks_.cdata();
    const auto* other_ptr = other.blocks_.cdata();

    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        if (block_ptr[blocki] != other_ptr[blocki])
        {
            return false;
        }
    }

    return true;
}


template<unsigned Width>
int Foam::PackedList<Width>::compare
(
    const PackedList<Width>& other
) const noexcept
{
    // Numeric-style ordering as most significant "digit", which matches
    // how elements are actually packed into blocks_
    // Returns -1/0/+1 for less/equal/greater.

    const label asize = this->size();
    const label bsize = other.size();

    if (!asize && !bsize)
    {
        // Both empty == equal
        return 0;
    }
    else if (!bsize)
    {
        // 'other' is empty
        return 1;
    }
    else if (!asize)
    {
        // 'this' is empty
        return -1;
    }
    else if (asize == bsize)
    {
        // Identical sizes - use block-wise comparison starting from the
        // most-significant (highest-index) block.
        // This is okay since trailing bits are always kept cleared

        const auto* a = this->blocks_.cdata();
        const auto* b = other.blocks_.cdata();

        for (label blocki = this->num_blocks()-1; blocki >= 0; --blocki)
        {
            if (a[blocki] < b[blocki]) return -1;
            if (a[blocki] > b[blocki]) return 1;
        }
        return 0;
    }
    else
    {
        // Different sizes - compare overlapping most-significant elements

        const label common = std::min(asize, bsize);

        label i = asize;
        label j = bsize;

        for (label k = 0; k < common; ++k)
        {
            const auto avalue = this->get(--i);
            const auto bvalue = other.get(--j);

            if (avalue < bvalue) return -1;
            if (avalue > bvalue) return 1;
        }

        // Identical common content, compare sizes
        return (asize < bsize) ? -1 : 1;
    }
}


template<unsigned Width>
Foam::labelList Foam::PackedList<Width>::values() const
{
    return this->unpack<label>();
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack() const
{
    static_assert
    (
        std::is_integral_v<IntType>,
        "Integral required for output."
    );
    static_assert
    (
        Width < std::numeric_limits<IntType>::digits,
        "Width of IntType is too small to hold result"
    );

    List<IntType> output(size());

    if (empty())
    {
        return output;
    }
    else if (uniform())
    {
        output = static_cast<IntType>(get(0));
        return output;
    }

    // NON-UNIFORM and len > 0

    label outi = 0;

    // Process n-1 complete blocks
    const label nblocks = num_blocks(size());

    for (label blocki=0; blocki < nblocks-1; ++blocki)
    {
        unsigned int blockval = blocks_[blocki];

        for (unsigned nget = elem_per_block; nget; --nget, ++outi)
        {
            output[outi] = IntType(blockval & PackedList<Width>::max_value);
            blockval >>= Width;
        }
    }

    // Any partial blocks
    for (/*nil*/; outi < size(); ++outi)
    {
        output[outi] = get(outi);
    }

    return output;
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack(const labelRange& range) const
{
    static_assert
    (
        std::is_integral_v<IntType>,
        "Integral required for unpack output."
    );
    static_assert
    (
        Width < std::numeric_limits<IntType>::digits,
        "Width of IntType is too small to hold unpack output."
    );


    // Could be more efficient but messier with block-wise access.
    // - automatically handles any invalid positions

    auto pos = range.start();

    List<IntType> output(range.size());

    for (IntType& out : output)
    {
        out = IntType(get(pos));
        ++pos;
    }

    return output;
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack(const labelUList& locations) const
{
    static_assert
    (
        std::is_integral_v<IntType>,
        "Integral required for unpack output."
    );
    static_assert
    (
        Width < std::numeric_limits<IntType>::digits,
        "Width of IntType is too small to hold unpack output."
    );


    label pos = 0;

    List<IntType> output(locations.size());

    for (IntType& out : output)
    {
        out = IntType(get(locations[pos]));
        ++pos;
    }

    return output;
}


template<unsigned Width>
void Foam::PackedList<Width>::set(const labelRange& range, unsigned int val)
{
    if (!val)
    {
        // Delegate to unset()
        unset(range);
        return;
    }

    labelRange slice(range);
    slice.adjust();  // No negative start, size adjusted accordingly

    // Range is invalid (zero-sized or entirely negative) - noop
    if (slice.empty())
    {
        return;
    }

    // Range finishes at or beyond the right side.
    // - zero fill any gaps that we might create.
    // - flood-fill the rest, which now corresponds to the full range.
    if (slice.end_value() >= size())
    {
        reserve(slice.end_value());
        resize(slice.begin_value(), 0u);
        resize(slice.end_value(), val);
        return;
    }

    fill_range(slice, val);
}


template<unsigned Width>
void Foam::PackedList<Width>::unset(const labelRange& range)
{
    // Require intersection with the current addressable range
    const labelRange slice(range.subset0(size()));

    // Range does not intersect (invalid, empty, or list is empty)
    if (slice.empty())
    {
        return;
    }

    // Range finishes at or beyond the right side.
    if (slice.end_value() >= size())
    {
        // The original size
        const label orig = size();

        resize(slice.begin_value(), 0u);
        resize(orig, 0u);
        return;
    }

    fill_range(slice, 0u);
}


template<unsigned Width>
void Foam::PackedList<Width>::assign
(
    const label offset,
    const PackedList<Width>& source
)
{
    if (source.empty())
    {
        // Nothing to do
        return;
    }

    // With clamped start/size
    const label lower = Foam::max(0, offset);
    const label upper = (offset + source.size());

    if (upper <= lower)
    {
        // No overlap
        return;
    }

    // Resize (growing only), also introducing trailing 0 bits since
    // this is meant to mimic a full bit copy - regardless of the
    // content of 'source'.
    if (this->size() < upper)
    {
        resize(upper);
    }

    // Uncommon case (negative offset clips the leading elements of 'source').
    // Handle simply/directly, element-by-element.
    if (offset < 0)
    {
        label srcPos = lower - offset;
        for (label i = lower; i < upper; ++i)
        {
            setValue_unsafe(i, source.get(srcPos));
            ++srcPos;
        }
        return;
    }

    // General case
    assign_at(offset, source);
}


// ************************************************************************* //
