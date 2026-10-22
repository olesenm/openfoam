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

\*---------------------------------------------------------------------------*/

#include "error.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<unsigned Width>
inline unsigned int Foam::PackedList<Width>::readValue(Istream& is)
{
    if constexpr (Width <= 1)
    {
        return readBool(is);
    }
    else
    {
        const unsigned int val = readLabel(is);

        if (val > max_value)
        {
            FatalIOErrorInFunction(is)
                << "Out-of-range value " << val << " for PackedList<" << Width
                << ">. Maximum permitted value is " << max_value << "."
                << exit(FatalIOError);
        }

        return val;
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::setPair(Istream& is)
{
    if constexpr (Width <= 1)
    {
        const label ind = readLabel(is);
        set(ind, true);
    }
    else
    {
        is.readBegin("Tuple2<label,uint32>");

        const label ind = readLabel(is);
        const unsigned int val = readLabel(is);

        is.readEnd("Tuple2<label,uint32>");

        if (val > max_value)
        {
            FatalIOErrorInFunction(is)
                << "Out-of-range value " << val << " for PackedList<" << Width
                << "> at index " << ind
                << ". Maximum permitted value is " << max_value << "."
                << exit(FatalIOError);
        }

        set(ind, val);

        is.check(FUNCTION_NAME);
    }
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::setValue_unsafe
(
    const label pos,
    const unsigned int val
) noexcept
{
    // Reference to the block content
    auto& ref_ = blocks_.data()[block_index(pos)];

    // The bit shift to access the given sub-portion
    const unsigned bshift = (Width * sub_index(pos));

    const unsigned int mask = (max_value << bshift);
    const unsigned int prev = ref_;

    if (val >= max_value)
    {
        ref_ |= mask;  // Overflow is max_value, so fill entirely
    }
    else
    {
        ref_ &= ~mask;
        ref_ |= mask & (val << bshift);
    }

    return (prev != ref_);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::clear_trailing_bits() noexcept
{
    // Mask off any partial rubbish in final block

    if (const unsigned int endOffset = sub_index(size()); endOffset)
    {
        unsigned int blk = block_index(size());

        blocks_.data()[blk] &= mask_lower(endOffset);
    }
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::trim(label minpos)
{
    if (empty())
    {
        return false;  // Trivial case
    }

    const label orig = size();
    if (orig < minpos)
    {
        minpos = orig;   // Don't allow allow accidental growth!
    }

    for (label blocki = this->num_blocks()-1; blocki >= 0; --blocki)
    {
        // Truncate to the block begin
        size_ = blocki * elem_per_block;

        if (unsigned int blockval = blocks_[blocki]; blockval)
        {
            // Some bits were found in the block, increment size again
            for (; blockval; ++size_)
            {
                blockval >>= Width;
            }
            break;
        }
        else if (size_ < minpos)
        {
            break;
        }
    }

    if (size_ < minpos)
    {
        size_ = minpos;
    }

    return (size() != orig);
}


template<unsigned Width>
inline Foam::label
Foam::PackedList<Width>::find_first_block() const noexcept
{
    if (size())
    {
        // Note: use cdata()[..] instead of blocks_[] to ensure noexcept
        // even when UList is compiled in full debug

        const label nblocks = this->num_blocks();
        const auto* block_ptr = blocks_.cdata();

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            if (block_ptr[blocki])  // ie, blocks_[blocki] but noexcept
            {
                return blocki;
            }
        }
    }

    return -1;
}


template<unsigned Width>
inline Foam::label
Foam::PackedList<Width>::find_first_not_block() const noexcept
{
    if (!size())
    {
        return -1;
    }

    // Check on complement (changes 0 <-> 1).
    // If any 1's now appear, there was a 0 bit before

    const label nblocks = this->num_blocks();
    const auto* block_ptr = blocks_.cdata();

    // Mask of the addressable bits within a fully packed block.
    constexpr unsigned int fullMask = mask_lower(elem_per_block);

    // All blocks preceding the final one are fully packed
    for (label blocki = 0; blocki < nblocks-1; ++blocki)
    {
        if (~(block_ptr[blocki]) & fullMask)
        {
            return blocki;
        }
    }

    // The final block: mask to the number of elements actually used,
    // or the full valid-bits mask if it happens to be completely filled
    const unsigned int endOffset = sub_index(size());
    const unsigned int emask = (endOffset ? mask_lower(endOffset) : fullMask);

    if (~(block_ptr[nblocks-1]) & emask)
    {
        return nblocks-1;
    }

    return -1;
}


template<unsigned Width>
inline Foam::label
Foam::PackedList<Width>::find_last_block() const noexcept
{
    if (size())
    {
        const auto* block_ptr = blocks_.cdata();

        for (label blocki = this->num_blocks()-1; blocki >= 0; --blocki)
        {
            if (block_ptr[blocki])  // ie, blocks_[blocki] but noexcept
            {
                return blocki;
            }
        }
    }

    return -1;
}


template<unsigned Width>
inline Foam::label
Foam::PackedList<Width>::find_last_not_block() const noexcept
{
    if (!size())
    {
        return -1;
    }

    // Check on complement (changes 0 <-> 1).
    // If any 1's now appear, there was a 0 bit before

    const label nblocks = this->num_blocks();
    const auto* block_ptr = blocks_.cdata();

    // Mask of the addressable bits within a fully packed block.
    constexpr unsigned int fullMask = mask_lower(elem_per_block);

    label blocki = (nblocks - 1);

    // The final block: mask to the number of elements actually used,
    // or the full valid-bits mask if it happens to be completely filled

    const unsigned int endOffset = sub_index(size());
    const unsigned int emask = (endOffset ? mask_lower(endOffset) : fullMask);

    {
        if (~(block_ptr[blocki]) & emask)
        {
            return blocki;
        }
        --blocki;
    }

    // Remaining (fully packed) blocks
    for (; blocki >= 0; --blocki)
    {
        if (~(block_ptr[blocki]) & fullMask)
        {
            return blocki;
        }
    }

    return -1;
}


template<unsigned Width>
inline void Foam::PackedList<Width>::fill_range
(
    const labelRange& slice,
    unsigned int val
)
{
    // Fill value for complete blocks
    const unsigned int blockval = repeated_value(val);

    // The more difficult case - everything in between.
    // 1. sequence may begin/end in the same block
    // 2. Cover more than one block
    //    a. with partial coverage in the first block
    //    b. with partial coverage in the end block

    // The begin block, bitmask for partial blocks (0 if unneeded)
    const unsigned int bblock = block_index(slice.begin_value());
    const unsigned int bmask  = mask_lower(sub_index(slice.begin_value()));

    // The end block, bitmask for partial blocks (0 if unneeded)
    const unsigned int eblock = block_index(slice.end_value());
    const unsigned int emask  = mask_lower(sub_index(slice.end_value()));

    if (bblock == eblock)
    {
        // Same block - fill between the begin/end bits.
        // Example:
        // bmask = 0000000000001111  (lower bits)
        // emask = 0000111111111111  (lower bits)
        // -> set  0000111111110000  (xor)

        const unsigned int m = (emask^bmask);
        blocks_[bblock] = ((blocks_[bblock] & ~m) | (blockval & m));
    }
    else
    {
        unsigned blocki = bblock;

        if (bmask)
        {
            // The first (partial) block
            // - only retain things below bmask, fill in the rest.
            blocks_[blocki] = ((blocks_[blocki] & bmask) | (blockval & ~bmask));
            ++blocki;
        }

        // Continue with block-wise fill
        for (; blocki < eblock; ++blocki)
        {
            blocks_[blocki] = blockval;
        }

        if (emask)
        {
            // The last (partial) block.
            // - only retain things above emask, fill in the rest.
            blocks_[eblock] = ((blocks_[eblock] & ~emask) | (blockval & emask));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned Width>
inline constexpr Foam::PackedList<Width>::PackedList() noexcept
:
    blocks_(),
    size_(0)
{}


template<unsigned Width>
inline Foam::PackedList<Width>::PackedList(const label numElem)
:
    blocks_(num_blocks(numElem), 0u),
    size_(numElem)
{}


template<unsigned Width>
inline Foam::PackedList<Width>::PackedList
(
    const label numElem,
    const unsigned int val
)
:
    blocks_(num_blocks(numElem), 0u),
    size_(numElem)
{
    if (val)
    {
        operator=(val);
    }
}


template<unsigned Width>
inline Foam::PackedList<Width>::PackedList(Istream& is)
:
    blocks_(),
    size_(0)
{
    readList(is);
}


template<unsigned Width>
inline Foam::PackedList<Width>::PackedList(const PackedList<Width>& list)
:
    blocks_(list.blocks_),
    size_(list.size_)
{}


template<unsigned Width>
inline Foam::PackedList<Width>::PackedList(PackedList<Width>&& list)
:
    blocks_(std::move(list.blocks_)),
    size_(list.size_)
{
    list.size_ = 0;
}


template<unsigned Width>
inline Foam::PackedList<Width>::PackedList(const labelUList& values)
:
    blocks_(num_blocks(values.size()), 0u),
    size_(values.size())
{
    const label len = values.size();

    // Could add more intelligent filling (blockwise),
    // but likely done fairly infrequently

    for (label i = 0; i < len; ++i)
    {
        const unsigned int val(values[i]);
        if (val) set(i, val);
    }
}


template<unsigned Width>
template<class Addr>
inline Foam::PackedList<Width>::PackedList
(
    const IndirectListBase<label, Addr>& values
)
:
    blocks_(num_blocks(values.size()), 0u),
    size_(values.size())
{
    const label len = values.size();

    // Could add more intelligent filling (blockwise),
    // but likely done fairly infrequently

    for (label i = 0; i < len; ++i)
    {
        const unsigned int val(values[i]);
        if (val) set(i, val);
    }
}


template<unsigned Width>
inline Foam::autoPtr<Foam::PackedList<Width>>
Foam::PackedList<Width>::clone() const
{
    return autoPtr<PackedList<Width>>::New(*this);
}


// * * * * * * * * * * * * * * * * References * * * * * * * * * * * * * * * * //

template<unsigned Width>
inline Foam::PackedList<Width>::reference::reference
(
    PackedList<Width>* parent,
    const label pos
) noexcept
:
    ref_(parent->blocks_.data()[block_index(pos)]),
    shift_(Width * sub_index(pos))
{}


template<unsigned Width>
inline unsigned int Foam::PackedList<Width>::reference::get() const noexcept
{
    return ((ref_ >> shift_) & max_value);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::reference::set
(
    const unsigned int val
) noexcept
{
    const unsigned int mask = (max_value << shift_);
    // const unsigned int prev = ref_;

    if (val >= max_value)
    {
        ref_ |= mask;  // Overflow is max_value, so fill entirely
    }
    else
    {
        ref_ &= ~mask;
        ref_ |= mask & (val << shift_);
    }

    // return (prev != ref_);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::reference::operator=
(
    const reference& other
) noexcept
{
    // Also safe for self-assignment
    this->set(other.get());
}


template<unsigned Width>
inline void Foam::PackedList<Width>::reference::operator=
(
    const unsigned int val
) noexcept
{
    this->set(val);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width>
inline void Foam::PackedList<Width>::checkIndex(const label i) const
{
    if (!size_)
    {
        FatalErrorInFunction
            << "attempt to access element " << i << " from zero sized list"
            << abort(FatalError);
    }
    else if (i < 0 || i >= size_)
    {
        FatalErrorInFunction
            << "index " << i << " out of range [0," << size_ << ")"
            << abort(FatalError);
    }
}


template<unsigned Width>
inline Foam::label Foam::PackedList<Width>::capacity() const noexcept
{
    return elem_per_block * blocks_.size();
}


template<unsigned Width>
inline void Foam::PackedList<Width>::resize_nocopy
(
    const label numElem
)
{
    this->resize(numElem);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::resize_fill
(
    const label numElem,
    unsigned int val
)
{
    if (numElem < 0)
    {
        // Bad input - clamp to zero silently or with warning
        #ifdef FULLDEBUG
        WarningInFunction
            << "clamping negative value for resize as 0" << endl;
        #endif

        blocks_ = 0u;  // Clear storage bits
        size_ = 0;
        return;
    }

    const label oldLen = num_blocks(size());
    const label newLen = num_blocks(numElem);

    // Extend and resize. Using reserve_exact() instead of reserve() to
    // minimize over-allocation
    reserve_exact(numElem);
    size_ = numElem;

    // If addressing space shrinks:
    // - clear old blocks so they remain with zero bits.

    for (label blocki = newLen; blocki < oldLen; ++blocki)
    {
        blocks_[blocki] = 0u;
    }

    // Fill addressable range
    fill(val);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::resize
(
    const label newSize,
    unsigned int val
)
{
    if (newSize < 0)
    {
        // Bad input - clamp to zero silently or with warning
        #ifdef FULLDEBUG
        WarningInFunction
            << "clamping negative value for resize as 0" << endl;
        #endif

        blocks_ = 0u;  // Clear storage bits
        size_ = 0;
        return;
    }

    const label oldSize = size();

    reserve(newSize);
    size_ = newSize;

    if (oldSize < newSize)
    {
        // Fill new elements or newly exposed elements
        if (val)
        {
            // Fill value for complete blocks
            const unsigned int blockval = repeated_value(val);

            // Fill complete blocks
            const label oldLen = num_blocks(oldSize);
            const label newLen = num_blocks(size());
            for (label blocki = oldLen; blocki < newLen; ++blocki)
            {
                blocks_[blocki] = blockval;
            }

            // Finish previous partial block, preserve existing value
            if (const unsigned int off = sub_index(oldSize); off)
            {
                const unsigned int blk = block_index(oldSize);
                const unsigned int mask = mask_lower(off);

                blocks_[blk] &= mask;
                blocks_[blk] |= ~mask & blockval;
            }

            clear_trailing_bits();
        }
    }
    else if (newSize < oldSize)
    {
        // The list is now shorter than before, so we zero assign the unused
        // blocks and any trailing junk. This costs slightly here, but make
        // things much simpler elsewhere.

        // Clear complete blocks
        const label oldLen = num_blocks(oldSize);
        const label newLen = num_blocks(size());
        for (label blocki = newLen; blocki < oldLen; ++blocki)
        {
            blocks_[blocki] = 0u;
        }

        clear_trailing_bits();
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::setCapacity(const label numElem)
{
    const label nblocks = num_blocks(numElem);

    blocks_.resize(nblocks, 0u);

    if (numElem < size())
    {
        size_ = numElem;
        clear_trailing_bits();
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::reserve(const label numElem)
{
    const label oldLen = blocks_.size();
    const label newLen = num_blocks(numElem);

    constexpr label minSize = 16;

    // Allocate more capacity if necessary
    if (oldLen < newLen)
    {
        blocks_.resize
        (
            // Allocation doubling with min-size
            Foam::max(minSize, Foam::max(newLen, 2*oldLen)),
            0u
        );
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::reserve_exact(const label numElem)
{
    const label oldLen = blocks_.size();
    const label newLen = num_blocks(numElem);

    constexpr label minSize = 16;

    // Allocate more capacity if necessary
    if (oldLen < newLen)
    {
        blocks_.resize
        (
            Foam::max(minSize, newLen),
            0u
        );
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::reset()
{
    blocks_ = 0u;
}


template<unsigned Width>
inline void Foam::PackedList<Width>::clear()
{
    reset();
    size_ = 0;
}


template<unsigned Width>
inline void Foam::PackedList<Width>::clearStorage()
{
    blocks_.clear();
    size_ = 0;
}


template<unsigned Width>
inline void Foam::PackedList<Width>::shrink_to_fit()
{
    // Any unneeded space allocated?
    if (const label nblocks = this->num_blocks(); nblocks < blocks_.size())
    {
        blocks_.resize(nblocks);
    }
}


template<unsigned Width>
inline Foam::label Foam::PackedList<Width>::num_blocks() const noexcept
{
    return num_blocks(size());
}


template<unsigned Width>
inline const char* Foam::PackedList<Width>::cdata_bytes() const noexcept
{
    return blocks_.cdata_bytes();
}


template<unsigned Width>
inline char* Foam::PackedList<Width>::data_bytes() noexcept
{
    return blocks_.data_bytes();
}


template<unsigned Width>
inline std::streamsize Foam::PackedList<Width>::size_data() const noexcept
{
    return num_blocks(size());
}


template<unsigned Width>
inline std::streamsize Foam::PackedList<Width>::size_bytes() const noexcept
{
    return num_blocks(size()) * sizeof(block_type);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::swap(PackedList<Width>& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-swap is a no-op
    }

    blocks_.swap(rhs.blocks_);
    std::swap(size_, rhs.size_);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::transfer(PackedList<Width>& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    blocks_.transfer(rhs.blocks_);
    size_ = rhs.size_;
    rhs.size_ = 0;
}


template<unsigned Width>
inline unsigned int Foam::PackedList<Width>::get(const label i) const noexcept
{
    if (i < 0 || i >= size())
    {
        return 0u;          // Out-of-bounds: return 0 (false)
    }
    else
    {
        // The corresponding block/offset
        const label blocki = block_index(i);
        const unsigned int off = sub_index(i);

        if constexpr (Width <= 1)
        {
            return
            (
                ((blocks_.cdata()[blocki]) >> (off)) & max_value
            );
        }
        else
        {
            return
            (
                ((blocks_.cdata()[blocki]) >> (Width * off)) & max_value
            );
        }
    }
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::set
(
    const label i,
    const unsigned int val
)
{
    if (i < 0)
    {
        #ifdef FULLDEBUG
        if (val)
        {
            WarningInFunction
                << "Ignoring attempt to set a negative index " << i
                << " range is [0," << size_ << ")"
                << endl;
        }
        #endif

        return false;       // Out-of-bounds: ignore
    }
    else if (i >= size())
    {
        if (val)
        {
            resize(i + 1);  // Lazy evaluation: adjust size for assign
        }
        else
        {
            // Ignore out-of-bounds for 'unset'
            return false;
        }
    }

    return this->setValue_unsafe(i, val);
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::unset(const label i)
{
    if (i < 0 || i >= size())
    {
        // Ignore out-of-bounds for 'unset'
        return false;
    }
    else
    {
        return this->setValue_unsafe(i, 0u);
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::push_back(const unsigned int val)
{
    const label idx = size();
    reserve(idx + 1);
    ++size_;

    if (val)
    {
        (void) this->setValue_unsafe(idx, val);
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::pop_back(label n)
{
    if (n >= size())
    {
        clear();
    }
    else if (n > 0)
    {
        resize(size() - n);
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::fill(const unsigned int val) noexcept
{
    if (empty())
    {
        return;  // Trivial case
    }

    const label nblocks = this->num_blocks();

    // Fill value for complete blocks
    const unsigned int blockval = (val ? repeated_value(val) : 0u);

    for (label blocki=0; blocki < nblocks; ++blocki)
    {
        blocks_.data()[blocki] = blockval;
    }

    if (val)
    {
        clear_trailing_bits();
    }
}

template<unsigned Width>
inline void Foam::PackedList<Width>::extract
(
    PackedList<Width>& target,
    label pos,
    label len
) const
{
    extract_copy(*this, target, pos, len);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::extract
(
    PackedList<Width>& target,
    const labelRange& range
) const
{
    extract_copy(*this, target, range);
}


template<unsigned Width>
inline Foam::PackedList<Width>
Foam::PackedList<Width>::extract(label pos, label len) const
{
    PackedList<Width> output;
    extract_copy(*this, output, pos, len);
    return output;
}


template<unsigned Width>
inline Foam::PackedList<Width>
Foam::PackedList<Width>::extract(const labelRange& range) const
{
    PackedList<Width> output;
    extract_copy(*this, output, range);
    return output;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<unsigned Width>
inline unsigned int
Foam::PackedList<Width>::operator[](const label i) const noexcept
{
    // Silently ignores out-of-bounds
    return this->get(i);
}


template<unsigned Width>
inline typename Foam::PackedList<Width>::reference
Foam::PackedList<Width>::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return reference(this, i);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::operator=(const PackedList<Width>& rhs)
{
    // Self-assignment silently ignored
    blocks_ = rhs.blocks_;
    size_ = rhs.size_;
}


template<unsigned Width>
inline void Foam::PackedList<Width>::operator=(PackedList<Width>&& rhs)
{
    transfer(rhs);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::operator=(const unsigned int val) noexcept
{
    fill(val);
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::operator==
(
    const PackedList<Width>& list
) const noexcept
{
    return this->equal(list);
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::operator!=
(
    const PackedList<Width>& list
) const noexcept
{
    return !(this->equal(list));
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::operator<
(
    const PackedList<Width>& other
) const noexcept
{
    return (this->compare(other) < 0);
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::operator>
(
    const PackedList<Width>& other
) const noexcept
{
    return (this->compare(other) > 0);
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::operator<=
(
    const PackedList<Width>& other
) const noexcept
{
    return (this->compare(other) <= 0);
}


template<unsigned Width>
inline bool Foam::PackedList<Width>::operator>=
(
    const PackedList<Width>& other
) const noexcept
{
    return (this->compare(other) >= 0);
}


// ************************************************************************* //
