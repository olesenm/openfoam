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

#include "bitSet.H"
#include "IOstreams.H"
#include "UPstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bitSet, 0);

    // TBD: add IO support of compound type?
    // defineNamedCompoundTypeName(bitSet, List<1>);
    // addNamedCompoundToRunTimeSelectionTable(bitSet, bitSet, List<1>);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::bitSet::minusEq(const bitSet& other)
{
    if (&other == this)
    {
        // Self '-=' : clears all bits
        if (debug & 2)
        {
            InfoInFunction
                << "Perform -= on self: clears all bits" << nl;
        }

        reset();
        return;
    }
    else if (none() || other.none())
    {
        // no-op: nothing can change
        return;
    }


    label nblocks(0);

    // Determine the overlap
    {
        const label thisExtent = (this->find_last()+1);
        const label thatExtent = (other.find_last()+1);

        // min overlap
        if (label common = std::min(thisExtent, thatExtent); common > 0)
        {
            nblocks = num_blocks(common);
        }
    }

    // The operation (on overlapping blocks)
    {
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] &= ~rhs[blocki];
        }
    }
}


void Foam::bitSet::andEq(const bitSet& other)
{
    if (FOAM_UNLIKELY(&other == this))
    {
        // Self '&=' : no-op

        if (debug & 2)
        {
            InfoInFunction
                << "Perform &= on self: ignore" << nl;
        }

        return;
    }
    else if (none())
    {
        // no-op: nothing is set - no intersection possible
        return;
    }
    else if (other.none())
    {
        // no-op: other has nothing set - no intersection possible
        reset();
        return;
    }


    label nblocks(0);

    // Determine the overlap
    {
        const label thisExtent = (this->find_last()+1);
        const label thatExtent = (other.find_last()+1);

        // min overlap
        if (label common = std::min(thisExtent, thatExtent); common > 0)
        {
            nblocks = num_blocks(common);
        }

        if (thatExtent < thisExtent)
        {
            // Clear bits (and blocks) that do not overlap at all
            const auto origSize = size();
            resize(thatExtent);
            resize(origSize);
        }
    }

    // The operation (on overlapping blocks)
    {
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] &= rhs[blocki];
        }
    }
}


void Foam::bitSet::orEq(const bitSet& other)
{
    if (&other == this)
    {
        // Self '|=' : no-op

        if (debug & 2)
        {
            InfoInFunction
                << "Perform |= on self: ignore" << nl;
        }

        return;
    }
    else if (other.none())
    {
        // no-op: nothing can change
        return;
    }


    label nblocks(0);

    // Determine the overlap
    {
        const label thisExtent = (this->find_last()+1);
        const label thatExtent = (other.find_last()+1);

        // max overlap
        if (label common = std::max(thisExtent, thatExtent); common > 0)
        {
            nblocks = num_blocks(common);
        }

        if (size() < thatExtent)
        {
            // Accommodate any extra bits from 'other'
            resize(thatExtent);
        }
    }

    // The operation (on overlapping blocks)
    {
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] |= rhs[blocki];
        }
    }
}


void Foam::bitSet::xorEq(const bitSet& other)
{
    if (&other == this)
    {
        // Self '^=' : clears all bits

        if (debug & 2)
        {
            InfoInFunction
                << "Perform ^= on self: clears all bits" << nl;
        }

        reset();
        return;
    }
    else if (other.none())
    {
        // no-op: nothing can change
        return;
    }


    label nblocks(0);

    // Determine the overlap
    {
        const label thisExtent = (this->find_last()+1);
        const label thatExtent = (other.find_last()+1);

        // max overlap
        if (label common = std::max(thisExtent, thatExtent); common > 0)
        {
            nblocks = num_blocks(common);
        }

        if (size() < thatExtent)
        {
            // Accommodate any extra bits from 'other'
            resize(thatExtent);
        }
    }

    // The operation (on overlapping blocks)
    {
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] ^= rhs[blocki];
        }
    }
}


void Foam::bitSet::left_shiftEq(label n)
{
    // Left-shift (moves content towards higher indices),
    // vacated low bits are cleared.
    // Bits shifted past the end are discarded and the size is unchanged.

    if (n >= size())
    {
        reset();
    }
    else if (n > 0)
    {
        // Left shift by N

        // Last addressable block, num_blocks() is >= 1
        const label last = num_blocks() - 1;

        // block_index(n)
        const label blkOffset = (n / elem_per_block);

        // sub_index(n)
        const auto bitShift = static_cast<unsigned int>(n % elem_per_block);

        if (bitShift)
        {
            const unsigned int invShift = (bits_per_block - bitShift);

            // Copy from lower to higher
            for (label i = last - blkOffset; i > 0; --i)
            {
                blocks_[i + blkOffset] =
                (
                    (blocks_[i] << bitShift) | (blocks_[i-1] >> invShift)
                );
            }
            blocks_[blkOffset] = (blocks_[0] << bitShift);
        }
        else
        {
            // Block-aligned: straight block-wise copy (high to low)
            for (label i = last - blkOffset; i > 0; --i)
            {
                blocks_[i + blkOffset] = blocks_[i];
            }

            blocks_[blkOffset] = blocks_[0];
        }

        // Vacated lower blocks become all-zero
        for (label i = 0; i < blkOffset; ++i)
        {
            blocks_[i] = 0u;
        }

        // Remove any bits that shifted in beyond the current size
        clear_trailing_bits();
    }
}


void Foam::bitSet::right_shiftEq(label n)
{
    // Right-shift (moves content from upper to lower indices),
    // vacated high bits are cleared.
    // Bits shifted past the end are discarded and the size is unchanged.

    if (n >= size())
    {
        reset();
    }
    else if (n > 0)
    {
        // Right shift by N

        // Last addressable block, num_blocks() is >= 1
        const label last = num_blocks() - 1;

        // block_index(n)
        const label blkOffset = (n / elem_per_block);

        // sub_index(n)
        const auto bitShift = static_cast<unsigned int>(n % elem_per_block);

        if (bitShift)
        {
            const unsigned int invShift = (bits_per_block - bitShift);

            // Copy from higher to lower
            for (label i = blkOffset; i < last; ++i)
            {
                blocks_[i - blkOffset] =
                (
                    (blocks_[i] >> bitShift) | (blocks_[i+1] << invShift)
                );
            }
            blocks_[last-blkOffset] = (blocks_[last] >> bitShift);
        }
        else
        {
            // Block-aligned: straight block-wise copy (low to high)
            for (label i = blkOffset; i <= last; ++i)
            {
                blocks_[i - blkOffset] = blocks_[i];
            }
        }

        // Vacated upper blocks become all-zero
        for (label i = (last - blkOffset + 1); i <= last; ++i)
        {
            blocks_[i] = 0u;
        }
    }
}


void Foam::bitSet::extract_copy
(
    const bitSet& source,
    bitSet& target,
    const labelRange& range
)
{
    if (&source == &target)
    {
        FatalErrorInFunction
            << "Cannot extract a copy into itself !"
            << abort(FatalError);
    }

    // With clamped range
    const label len = range.subset0(source.size()).size();

    if (len <= 0)
    {
        // Nothing to extract
        target.clear();
    }
    // Trivial cases - avoid any block-wise copying/shifting
    else if (source.empty() || source.none())
    {
        // All bits false
        target.resize_fill(len, false);
    }
    else if (source.all())
    {
        // All bits true
        target.resize_fill(len, true);
    }
    else
    {
        // General case: block-wise extraction with sub-block bit shifting
        PackedList<1>::extract_copy(source, target, range);
    }
}


void Foam::bitSet::extract_copy
(
    const bitSet& source,
    bitSet& target,
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

Foam::bitSet::bitSet(Istream& is)
:
    PackedList<1>()
{
    is  >> *this;
}


Foam::bitSet::bitSet(const bitSet& bitset, const labelUList& addr)
:
    bitSet(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, bitset.get(addr[i]));
    }
}


Foam::bitSet::bitSet(const bitSet& source, const labelRange& range)
:
    bitSet()
{
    extract_copy(source, *this, range);
}


Foam::bitSet::bitSet(const label minSize, const bitSet& bitset)
:
    bitSet(bitset)
{
    if (bitset.size() < minSize)
    {
        // Need to extend the size
        resize(minSize);
    }
    else if (bitset.size() > minSize)
    {
        // Avoid over allocated size
        const label extent = (this->find_last()+1);

        // Shrink to the actually used length or the expected length
        resize(std::max(extent, minSize));
    }
}


Foam::bitSet::bitSet(const label n, const labelRange& range)
:
    bitSet(n)
{
    this->set(range);
}


Foam::bitSet::bitSet(const labelRange& range)
{
    this->set(range);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::bitSet::assign(const UList<bool>& bools)
{
    resize_fill(bools.size(), 0u);

    unsigned bitIdx = 0u;
    auto* packed = blocks_.data();

    // Set according to indices that are true
    for (const auto b : bools)
    {
        if (b)
        {
            *packed |= (1u << bitIdx);
        }

        if (++bitIdx >= elem_per_block)
        {
            bitIdx = 0u;
            ++packed;
        }
    }
}


bool Foam::bitSet::intersects(const bitSet& other) const
{
    if (size() && other.size())
    {
        const label nblocks = num_blocks(std::min(size(), other.size()));
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            if (bool(blocks_[blocki] & rhs[blocki]))
            {
                return true;
            }
        }
    }

    return false;
}


void Foam::bitSet::unset(const label offset, const bitSet& other)
{
    if (other.none() || this->none())
    {
        // no-op: nothing to clear, or nothing can change
        return;
    }
    else if (offset == 0)
    {
        unset(other);
        return;
    }

    // one-past last active positon
    const label thisExtent = (this->find_last()+1);
    const label thatExtent = (other.find_last()+1);

    // With clamped start/size.
    // - clearing bits beyond the current size is a no-op anyhow.
    const label lower = Foam::max(0, offset);
    const label upper = Foam::min(thisExtent, (offset + thatExtent));

    if (upper <= lower)
    {
        // No overlap
        return;
    }

    if (other.all())
    {
        unset(labelRange(lower, (upper-lower)));
        return;
    }

    // Uncommon case (negative offset clips the leading bits of 'other')
    //
    // Handle simply/directly, bit-by-bit.
    if (offset < 0)
    {
        for (label i : other)
        {
            unset(i + offset);  // silently handles out-of-range
        }

        return;
    }

    // ---------------

    // Process block-wise as much as possible, which may need
    // to carry values between blocks.
    // Note that clearing (AND NOT) does not disturb bits lying
    // outside of the target range but within the same block(s).

    // The actually usable source length:
    const label nblocks =
        num_blocks(Foam::min(thatExtent, (this->size() - lower)));

    const label bblock   = block_index(lower);
    const unsigned bmask = sub_index(lower);

    // The input start/end
    auto iter = other.blocks_.cbegin();
    const auto last = (iter + nblocks);

    auto output = (blocks_.begin() + bblock);

    if (!bmask)
    {
        // Fast path: block-aligned start - direct AND NOT, no shift
        while (iter != last)
        {
            *output &= ~(*iter);

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

            *output &= ~(carry | (val << bmask));
            ++output;

            carry = (val >> (elem_per_block - bmask));
        }

        if (carry && output != endOutput)
        {
            *output &= ~carry;
        }
    }
}


void Foam::bitSet::assign(const label offset, const bitSet& source)
{
    if (source.empty())
    {
        // Nothing to do
        return;
    }

    // With clamped start/size
    const label lower = Foam::max(0, offset);
    const label upper = (offset + source.size());

    // TBD: avoid self assignment?

    if (upper <= lower)
    {
        // No overlap
        return;
    }

    // Resize (growing only), also introducing trailing 0 bits since
    // this is meant to mimic a full bit copy - regardless of the
    // content of 'source' (ie, even if all zero).
    if (this->size() < upper)
    {
        resize(upper);
    }

    if (source.all())
    {
        // Trivial case
        set(labelRange(lower, (upper-lower)));
        return;
    }

    // Clear the target range first and assign by OR-ing  bits from 'source'
    unset(labelRange(lower, (upper-lower)));

    if (source.none())
    {
        // Already cleared the bits just above
        return;
    }

    // Uncommon case (negative offset clips the leading bits of 'source').
    // Handle simply/directly, bit-by-bit.
    if (offset < 0)
    {
        for (label i : source)
        {
            if (auto x = i + offset; x >= 0)
            {
                set(x);
            }
        }
        return;
    }

    // General case: reuse the block-wise merge from the base.
    // The target range was already cleared above, so a plain OR-merge
    // (without re-clearing) would suffice, but assign_at() re-clearing
    // the same (already-zero) range is harmless and keeps a single,
    // well-tested implementation.
    PackedList<1>::assign_at(offset, source);
}


Foam::labelList Foam::bitSet::toc(const label offset) const
{
    // Number of used (set) entries
    const label total = any() ? count() : 0;

    if (!total)
    {
        return labelList();
    }

    labelList output(total);
    label nItem = 0;

    // Process block-wise, detecting any '1' bits

    const label nblocks = num_blocks(size());
    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        if (unsigned int blockval = blocks_[blocki]; blockval)
        {
            for (label pos = (blocki * elem_per_block); blockval; ++pos)
            {
                if (blockval & 1u)
                {
                    output[nItem] = (pos + offset);
                    ++nItem;
                }
                blockval >>= 1u;
            }
            if (nItem == total) break;  // Terminate early
        }
    }

    return output;
}


Foam::List<bool> Foam::bitSet::values() const
{
    List<bool> output(size(), false);

    // Process block-wise, detecting any '1' bits

    const label nblocks = num_blocks(size());
    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        label pos = (blocki * elem_per_block);

        for
        (
            unsigned int blockval = blocks_[blocki];
            blockval;
            blockval >>= 1u
        )
        {
            if (blockval & 1u)
            {
                output[pos] = true;
            }
            ++pos;
        }
    }

    return output;
}


// * * * * * * * * * * * * * *  Parallel Functions * * * * * * * * * * * * * //

// Special purpose broadcast for bitSet which is more efficient than
// either a regular broadcast (with serialization) or the usual
// broadcast for lists.
//
// The initial broadcast sends the extent (last bit on) and the length.
// The receive clears out its bits and resizes (ie, all zeros and the proper
// length).
// The final broadcast can be skipped (if there is no non-zero content),
// or simply minimized the amount of data broadcast.

void Foam::bitSet::broadcast(int communicator, int root)
{
    if (communicator < 0)
    {
        communicator = UPstream::worldComm;
    }

    if (!UPstream::is_parallel(communicator))
    {
        return;
    }

    // The number of data blocks for the operation
    label nblocks(0);

    // Determine the extent/sizing
    {
        int64_t sizing[2] = { 0, 0 };

        if (root == UPstream::myProcNo(communicator))
        {
            // Sender: extent of content and size
            sizing[0] = static_cast<int64_t>(find_last()+1);
            sizing[1] = static_cast<int64_t>(size());

            UPstream::broadcast(sizing, 2, communicator, root);
        }
        else
        {
            // Receiver: use a clean bitset with the updated size
            UPstream::broadcast(sizing, 2, communicator, root);

            clear();  // Clear old contents
            resize(sizing[1]);
        }

        nblocks = num_blocks(sizing[0]);
    }

    if (nblocks > 0)
    {
        UPstream::broadcast(this->data(), nblocks, communicator, root);
    }
}


void Foam::bitSet::reduceAnd(int communicator, bool syncSizes)
{
    if (communicator < 0)
    {
        communicator = UPstream::worldComm;
    }

    if (!UPstream::is_parallel(communicator))
    {
        return;
    }

    // The number of data blocks for the operation
    label nblocks(0);

    // Operation is an intersection.
    // - can restrict to a smaller region than the original size
    if (syncSizes)
    {
        int64_t common(find_last()+1);  // Size to include last bit

        UPstream::mpiAllReduce_min(&common, 1, communicator);

        nblocks = num_blocks(common);

        if (common > 0)
        {
            // Clear bits (and blocks) that do not overlap at all
            const auto origSize = size();
            resize(common);
            resize(origSize);
        }
        else
        {
            // No intersections
            reset();
        }
    }
    else
    {
        nblocks = this->num_blocks();
    }

    if (nblocks > 0)
    {
        UPstream::mpiAllReduce<UPstream::opCodes::op_bit_and>
        (
            this->data(),
            nblocks,
            communicator
        );
    }
}


void Foam::bitSet::reduceOr(int communicator, bool syncSizes)
{
    if (communicator < 0)
    {
        communicator = UPstream::worldComm;
    }

    if (!UPstream::is_parallel(communicator))
    {
        return;
    }

    // The number of data blocks needed for the operation
    label nblocks(0);

    // Operation can increase the addressed size
    if (syncSizes)
    {
        // Operational sizes
        int64_t sizing[2] =
        {
            static_cast<int64_t>(find_last()+1),  // Size to include last bit
            static_cast<int64_t>(size())          // The overall size
        };

        UPstream::mpiAllReduce_max(sizing, 2, communicator);

        nblocks = num_blocks(sizing[0]);

        // Extend local size to the max size encountered
        // This is greedy, but produces consistent sizing
        if (size() < sizing[1])
        {
            resize(sizing[1]);
        }

        // Alternative:
        // Extend local size to include any 'on' bits.
        // The resulting bitsets will not have identical sizes on all ranks.
        // if (size() < sizing[0])
        // {
        //     resize(sizing[0]);
        // }
    }
    else
    {
        nblocks = this->num_blocks();
    }

    if (nblocks > 0)
    {
        UPstream::mpiAllReduce<UPstream::opCodes::op_bit_or>
        (
            this->data(),
            nblocks,
            communicator
        );
    }
}


Foam::bitSet Foam::bitSet::gatherValues(bool localValue, int communicator)
{
    if (communicator < 0)
    {
        communicator = UPstream::worldComm;
    }

    bitSet allValues;

    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(communicator) as well?
        allValues.resize(1);
        allValues.set(0, localValue);
    }
    else
    {
        List<bool> bools;
        if (UPstream::master(communicator))
        {
            bools.resize(UPstream::nProcs(communicator), false);
        }

        UPstream::mpiGather
        (
            &localValue,    // Send
            bools.data(),   // Recv
            1,              // Num send/recv data per rank
            communicator
        );

        // Transcribe to bitSet (on master)
        allValues.assign(bools);
    }

    return allValues;
}


// Note that for allGather()
// - MPI_Gather of individual bool values and broadcast the packed result
// - this avoids bit_or on 32bit values everywhere, since we know a priori
//   that each rank only contributes 1bit of info

Foam::bitSet Foam::bitSet::allGather(bool localValue, int communicator)
{
    if (communicator < 0)
    {
        communicator = UPstream::worldComm;
    }

    bitSet allValues(bitSet::gatherValues(localValue, communicator));

    if (UPstream::is_parallel(communicator))
    {
        // Identical size on all ranks
        allValues.resize(UPstream::nProcs(communicator));

        UPstream::broadcast
        (
            allValues.data(),
            allValues.num_blocks(),
            communicator
            // root = 0
        );
    }

    return allValues;
}


// ************************************************************************* //
