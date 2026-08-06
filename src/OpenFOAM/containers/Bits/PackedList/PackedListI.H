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
    auto& ref_ = blocks_.data()[pos / elem_per_block];

    // The bit shift to access the given sub-portion
    const unsigned shift = (Width * (pos % elem_per_block));

    const unsigned int mask = (max_value << shift);
    const unsigned int prev = ref_;

    if (val >= max_value)
    {
        ref_ |= mask;  // Overflow is max_value, so fill entirely
    }
    else
    {
        ref_ &= ~mask;
        ref_ |= mask & (val << shift);
    }

    return (prev != ref_);
}


template<unsigned Width>
inline void Foam::PackedList<Width>::clear_trailing_bits() noexcept
{
    // Mask off any partial rubbish in final block

    if (const unsigned int off = size() % elem_per_block; off)
    {
        const unsigned int blk = size() / elem_per_block;

        blocks_.data()[blk] &= mask_lower(off);
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

    for (label blocki = num_blocks(size())-1; blocki >= 0; --blocki)
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
inline Foam::label Foam::PackedList<Width>::first_block() const noexcept
{
    if (size())
    {
        // Note: use cdata()[..] instead of blocks_[] to ensure noexcept
        // even when UList is compiled in full debug

        const label nblocks = num_blocks(size());
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
inline Foam::label Foam::PackedList<Width>::first_not_block() const noexcept
{
    if (!size())
    {
        return -1;
    }

    // Check on complement (changes 0 <-> 1).
    // If any 1's now appear, there was a 0 bit before

    const label nblocks = num_blocks(size());
    const auto* block_ptr = blocks_.cdata();

    // Extra bits in the final block?
    if (const unsigned int off = size() % elem_per_block; off)
    {
        for (label blocki = 0; blocki < nblocks-1; ++blocki)
        {
            if (~(block_ptr[blocki]))
            {
                return blocki;
            }
        }

        // The final block needs masking
        if (~(block_ptr[nblocks-1]) & mask_lower(off))
        {
            return nblocks-1;
        }
    }
    else
    {
        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            if (~(block_ptr[blocki]))
            {
                return blocki;
            }
        }
    }

    return -1;
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
    ref_(parent->blocks_.data()[pos / elem_per_block]),
    shift_(Width * (pos % elem_per_block))
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
inline void Foam::PackedList<Width>::resize
(
    const label newSize,
    const unsigned int val
)
{
    reserve(newSize);

    const label oldSize = size();
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
            if (const unsigned int off = oldSize % elem_per_block; off)
            {
                const unsigned int blk = oldSize / elem_per_block;
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

    // Allocate more capacity if necessary
    if (oldLen < newLen)
    {
        blocks_.resize
        (
            // SizeMin=16, allocation doubling
            Foam::max(16, Foam::max(newLen, 2*oldLen)),
            0u
        );
    }
}


template<unsigned Width>
inline void Foam::PackedList<Width>::reserve_exact(const label numElem)
{
    const label oldLen = blocks_.size();
    const label newLen = num_blocks(numElem);

    // Allocate more capacity if necessary
    if (oldLen < newLen)
    {
        blocks_.resize
        (
            // SizeMin=16
            Foam::max(16, newLen),
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
inline void Foam::PackedList<Width>::reset(const label newSize)
{
    // Clear all storage bits
    blocks_ = 0u;

    if (newSize >= 0)
    {
        reserve(newSize);  // Any extended storage filled with 0u
        size_ = newSize;   // ... can simply adjust addressed size
    }
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
    const label nblocks = num_blocks(size());
    if (nblocks < blocks_.size())
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
        const label blocki = i / elem_per_block;
        const unsigned int off = i % elem_per_block;

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
inline unsigned int Foam::PackedList<Width>::remove()
{
    // Location of last element and simultaneously the new size
    const label idx = size()-1;
    const unsigned int old = this->get(idx);

    if (idx < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    resize(idx);
    return old;
}


template<unsigned Width>
inline void Foam::PackedList<Width>::fill(const unsigned int val) noexcept
{
    if (empty())
    {
        return;  // Trivial case
    }

    const label nblocks = num_blocks(size());

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


// ************************************************************************* //
