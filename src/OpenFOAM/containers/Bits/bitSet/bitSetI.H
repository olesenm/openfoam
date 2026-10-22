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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::bitSet::bitSet() noexcept
:
    PackedList<1>()
{}


inline Foam::bitSet::bitSet(const label n)
:
    PackedList<1>(n)
{}


inline Foam::bitSet::bitSet(const label n, const bool val)
:
    bitSet(n)
{
    if (val) fill(val);
}


inline Foam::bitSet::bitSet(const bitSet& bitset)
:
    PackedList<1>(bitset)
{}


inline Foam::bitSet::bitSet(bitSet&& bitset)
:
    PackedList<1>(std::move(bitset))
{}


inline Foam::bitSet::bitSet(PackedList<1>&& list) noexcept
:
    PackedList<1>(std::move(list))
{}


inline Foam::bitSet::bitSet(const UList<bool>& bools)
:
    bitSet()
{
    assign(bools);
}


inline Foam::bitSet::bitSet(const label n, const labelUList& locations)
:
    bitSet(n)
{
    setMany(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::bitSet::bitSet
(
    const label n,
    const IndirectListBase<label, Addr>& locations
)
:
    bitSet(n)
{
    setMany(locations.begin(), locations.end());
}


inline Foam::bitSet::bitSet
(
    const label n,
    std::initializer_list<label> locations
)
:
    bitSet(n)
{
    setMany(locations.begin(), locations.end());
}


inline Foam::bitSet::bitSet(const labelUList& locations)
:
    bitSet()
{
    setMany(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::bitSet::bitSet
(
    const IndirectListBase<label, Addr>& locations
)
:
    bitSet()
{
    setMany(locations.begin(), locations.end());
}


inline Foam::autoPtr<Foam::bitSet> Foam::bitSet::clone() const
{
    return autoPtr<bitSet>::New(*this);
}


// * * * * * * * * * * * * * * * * References * * * * * * * * * * * * * * * * //

inline Foam::bitSet::reference::reference
(
    bitSet* parent,
    const label pos
) noexcept
:
    PackedList<1>::reference(parent, pos)
{}


inline void Foam::bitSet::reference::flip() noexcept
{
    const unsigned int mask = (max_value << shift_);
    ref_ ^= mask;
}


inline void Foam::bitSet::reference::operator=
(
    const reference& other
) noexcept
{
    // Also safe for self-assignment
    this->set(other.get());
}


inline void Foam::bitSet::reference::operator=
(
    const unsigned int val
) noexcept
{
    this->set(val);
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline constexpr Foam::bitSet::const_iterator::const_iterator() noexcept
:
    set_(nullptr),
    pos_(-1)
{}


inline Foam::bitSet::const_iterator::const_iterator
(
    const bitSet* parent
) noexcept
:
    set_(parent),
    pos_(set_->find_first())
{}


inline Foam::bitSet::const_iterator::const_iterator
(
    const bitSet* parent,
    label pos
) noexcept
:
    set_(parent),
    pos_(parent->find_next(pos-1))
{}


inline Foam::label Foam::bitSet::const_iterator::operator*() const noexcept
{
    return pos_;
}


inline Foam::bitSet::const_iterator&
Foam::bitSet::const_iterator::operator++() noexcept
{
    pos_ = set_->find_next(pos_);
    return *this;
}


inline bool Foam::bitSet::const_iterator::operator==
(
    const const_iterator& iter
) const noexcept
{
    return (iter.pos_ == pos_);
}


inline bool Foam::bitSet::const_iterator::operator!=
(
    const const_iterator& iter
) const noexcept
{
    return (iter.pos_ != pos_);
}


inline Foam::bitSet::const_iterator
Foam::bitSet::begin() const noexcept
{
    return const_iterator(this);
}


inline Foam::bitSet::const_iterator
Foam::bitSet::cbegin() const noexcept
{
    return const_iterator(this);
}


inline Foam::bitSet::const_iterator
Foam::bitSet::begin(label pos) const noexcept
{
    return const_iterator(this, pos);
}


inline Foam::bitSet::const_iterator
Foam::bitSet::cbegin(label pos) const noexcept
{
    return const_iterator(this, pos);
}


inline Foam::bitSet::const_iterator
Foam::bitSet::end() const noexcept
{
    return const_iterator();
}


inline Foam::bitSet::const_iterator
Foam::bitSet::cend() const noexcept
{
    return const_iterator();
}


inline Foam::label Foam::bitSet::find_first(label pos) const noexcept
{
    if (pos >= size())
    {
        // out-of-range
        return -1;
    }

    if (pos <= 0)
    {
        // Special case (search from the beginning): use block-wise search

        if (const label blocki = find_first_block(); blocki >= 0)
        {
            label index_found = (blocki * elem_per_block);

            // Detect first '1' bit within the block
            for
            (
                unsigned int blockval = blocks_.cdata()[blocki];
                blockval;
                blockval >>= 1u
            )
            {
                if (blockval & 1u)
                {
                    return index_found;
                }
                ++index_found;
            }
        }

        return -1;
    }

    const label nblocks = this->num_blocks();

    // The corresponding block/offset
    label blocki = block_index(pos);

    // Check the first (partial) block, if needed.
    if (unsigned int off = sub_index(pos); off > 0)
    {
        label index_found = pos;

        // No extra size() guard needed on index_found
        for
        (
            unsigned int blockval = (blocks_.cdata()[blocki] >> off);
            blockval;
            blockval >>= 1u
        )
        {
            if (blockval & 1u)
            {
                return index_found;
            }
            ++index_found;
        }

        // Continue with the next block
        ++blocki;
    }

    // Block-wise search
    for (; blocki < nblocks; ++blocki)
    {
        label index_found = (blocki * elem_per_block);

        // No extra size() guard needed on index_found
        for
        (
            unsigned int blockval = blocks_.cdata()[blocki];
            blockval;
            blockval >>= 1u
        )
        {
            if (blockval & 1u)
            {
                return index_found;
            }
            ++index_found;
        }
    }

    return -1;
}


inline Foam::label Foam::bitSet::find_next(label pos) const noexcept
{
    if (pos >= -1 && pos < size())
    {
        // Search one position beyond the specified one
        return find_first(pos+1);
    }
    else
    {
        // out-of-range
        return -1;
    }
}


inline Foam::label Foam::bitSet::find_first_not(label pos) const noexcept
{
    if (pos >= size())
    {
        // out-of-range
        return -1;
    }

    if (pos <= 0)
    {
        // Special case (search from the beginning): use block-wise search

        if (const label blocki = find_first_not_block(); blocki >= 0)
        {
            label index_found = (blocki * elem_per_block);

            // Detect first '0' bit within the block (check the complement)
            // No special masking - already done by find_first_not_block()

            for
            (
                unsigned int blockval = ~(blocks_.cdata()[blocki]);
                blockval;
                blockval >>= 1u
            )
            {
                if (blockval & 1u)
                {
                    return index_found;
                }
                ++index_found;
            }
        }

        return -1;
    }

    const label nblocks = this->num_blocks();

    // Mask of the valid (in-range) bits in the final block.
    // == 0 if no masking is required
    const unsigned int emask = mask_lower(sub_index(size()));

    // The corresponding block/offset
    label blocki = block_index(pos);

    // Check the first (partial) block, if needed.
    if (unsigned int off = sub_index(pos); off > 0)
    {
        label index_found = pos;
        unsigned int blockval = ~(blocks_.cdata()[blocki]);

        // Mask out any out-of-range bits for the final block
        if (emask && blocki == nblocks-1)
        {
            blockval &= emask;
        }

        for (blockval >>= off; blockval; blockval >>= 1u)
        {
            if (blockval & 1u)
            {
                return index_found;
            }
            ++index_found;
        }

        // Continue with the next block
        ++blocki;
    }

    // Block-wise search
    for (; blocki < nblocks; ++blocki)
    {
        label index_found = (blocki * elem_per_block);
        unsigned int blockval = ~(blocks_.cdata()[blocki]);

        // Mask out any out-of-range bits for the final block
        if (emask && blocki == nblocks-1)
        {
            blockval &= emask;
        }

        for (; blockval; blockval >>= 1u)
        {
            if (blockval & 1u)
            {
                return index_found;
            }
            ++index_found;
        }
    }

    return -1;
}


inline Foam::label Foam::bitSet::find_next_not(label pos) const noexcept
{
    if (pos >= -1 && pos < size())
    {
        // Search one position beyond the specified one
        return find_first_not(pos+1);
    }
    else
    {
        // out-of-range
        return -1;
    }
}


inline Foam::label Foam::bitSet::find_last() const noexcept
{
    if (const label blocki = find_last_block(); blocki >= 0)
    {
        if (unsigned int blockval = blocks_.cdata()[blocki]; blockval)
        {
            label pos = (blocki * elem_per_block) - 1;

            while (blockval)
            {
                blockval >>= 1u;
                ++pos;
            }

            return pos;
        }
    }

    return -1;
}


inline Foam::label Foam::bitSet::find_last_not() const noexcept
{
    if (const label blocki = find_last_not_block(); blocki >= 0)
    {
        unsigned int blockval = ~(blocks_.cdata()[blocki]);

        // Mask out any out-of-range bits for the final block
        if (blocki == this->num_blocks()-1)
        {
            if (unsigned int emask = mask_lower(sub_index(size())); emask)
            {
                blockval &= emask;
            }
        }

        label index_found = (blocki * elem_per_block) - 1;

        while (blockval)
        {
            blockval >>= 1u;
            ++index_found;
        }

        return index_found;
    }

    return -1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::bitSet::all() const noexcept
{
    if (empty()) return true;  // SIC. boost convention
    return (-1 == find_first_not_block());
}


inline bool Foam::bitSet::any() const noexcept
{
    return (-1 != find_first_block());
}


inline bool Foam::bitSet::none() const noexcept
{
    return (-1 == find_first_block());
}


inline bool Foam::bitSet::empty_or_test(label i) const noexcept
{
    return (empty() || test(i));
}


inline unsigned int Foam::bitSet::count(const bool on) const noexcept
{
    unsigned int total = 0;

    const label nblocks = num_blocks(size());

    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        total += BitOps::bit_count(blocks_.cdata()[blocki]);
    }

    if (!on)
    {
        // Return the number of bits that are OFF.
        return (unsigned(size()) - total);
    }

    return total;
}


inline Foam::labelList Foam::bitSet::sortedToc(const label offset) const
{
    return toc(offset);
}


inline void Foam::bitSet::resize_last()
{
    if (const label pos = find_last(); pos >= 0)
    {
        resize(pos+1);
    }
    else
    {
        clear();
    }
}


inline void Foam::bitSet::swap(bitSet& bitset)
{
    PackedList<1>::swap(bitset);
}


inline void Foam::bitSet::transfer(bitSet& bitset)
{
    PackedList<1>::transfer(bitset);
}


inline void Foam::bitSet::fill(const bool val)
{
    if (empty())
    {
        return;
    }
    else if (val)
    {
        std::fill_n(blocks_.data(), num_blocks(), ~0u);
        clear_trailing_bits();
    }
    else
    {
        // or: PackedList<1>::reset();
        // or: blocks_ = 0u;
        std::fill_n(blocks_.data(), num_blocks(), 0u);
    }
}


inline bool Foam::bitSet::test_set(const label i, const bool val)
{
    bool old = this->test(i);
    PackedList<1>::set(i, val);
    return old;
}


inline void Foam::bitSet::set(const bitSet& bitset)
{
    orEq(bitset);
}


inline Foam::label Foam::bitSet::set(const labelUList& locations)
{
    return setMany(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::label Foam::bitSet::set
(
    const IndirectListBase<label, Addr>& locations
)
{
    return setMany(locations.begin(), locations.end());
}


inline Foam::label Foam::bitSet::unset(const labelUList& locations)
{
    return unset(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::label Foam::bitSet::unset
(
    const IndirectListBase<label, Addr>& locations
)
{
    return unset(locations.begin(), locations.end());
}


inline void Foam::bitSet::unset(const bitSet& other)
{
    minusEq(other);
}


inline void Foam::bitSet::flip()
{
    if (size())
    {
        const label nblocks = num_blocks(size());
        const auto* end_block = (blocks_.cdata() + nblocks);

        for (auto* blk = blocks_.data(); blk != end_block; ++blk)
        {
            *blk = ~(*blk);
        }
        clear_trailing_bits();
    }
}


inline void Foam::bitSet::flip(const label i)
{
    if (i >= 0 && i < size())
    {
        reference(this, i).flip();
    }
}


inline void Foam::bitSet::extract(bitSet& target, label pos, label len) const
{
    extract_copy(*this, target, pos, len);
}


inline void Foam::bitSet::extract(bitSet& target, const labelRange& range) const
{
    extract_copy(*this, target, range);
}


inline Foam::bitSet Foam::bitSet::extract(label pos, label len) const
{
    bitSet output;
    extract_copy(*this, output, pos, len);
    return output;
}


inline Foam::bitSet Foam::bitSet::extract(const labelRange& range) const
{
    bitSet output;
    extract_copy(*this, output, range);
    return output;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::bitSet::operator()(const label pos) const noexcept
{
    return this->test(pos);
}


inline unsigned int Foam::bitSet::operator[](const label i) const noexcept
{
    return this->get(i);
}


inline Foam::bitSet::reference Foam::bitSet::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return reference(this, i);
}


inline Foam::bitSet& Foam::bitSet::operator=(const bitSet& bitset)
{
    PackedList<1>::operator=(bitset);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator=(bitSet&& bitset)
{
    transfer(bitset);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator=(const bool val)
{
    fill(val);
    return *this;
}


inline Foam::bitSet Foam::bitSet::operator~() const
{
    bitSet result(*this);
    result.flip();
    return result;
}


inline Foam::bitSet& Foam::bitSet::operator&=(const bitSet& other)
{
    andEq(other);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator|=(const bitSet& other)
{
    orEq(other);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator^=(const bitSet& other)
{
    xorEq(other);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator-=(const bitSet& other)
{
    minusEq(other);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator<<=(label n)
{
    left_shiftEq(n);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator>>=(label n)
{
    right_shiftEq(n);
    return *this;
}


inline Foam::bitSet Foam::bitSet::operator<<(label n) const
{
    bitSet result(*this);
    result <<= n;
    return result;
}


inline Foam::bitSet Foam::bitSet::operator>>(label n) const
{
    bitSet result(*this);
    result >>= n;
    return result;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline Foam::bitSet Foam::operator&(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    result &= b;

    result.resize_last();
    return result;
}


inline Foam::bitSet Foam::operator|(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    result |= b;

    result.resize_last();
    return result;
}


inline Foam::bitSet Foam::operator^(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    result ^= b;

    result.resize_last();
    return result;
}


inline Foam::bitSet Foam::operator-(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    result |= b;

    result.resize_last();
    return result;
}


// ************************************************************************* //
