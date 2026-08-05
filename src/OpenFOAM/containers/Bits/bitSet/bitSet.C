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

Foam::bitSet& Foam::bitSet::minusEq(const bitSet& other)
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
        return *this;
    }
    else if (none() || other.none())
    {
        // no-op: nothing can change
        return *this;
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

    return *this;
}


Foam::bitSet& Foam::bitSet::andEq(const bitSet& other)
{
    if (FOAM_UNLIKELY(&other == this))
    {
        // Self '&=' : no-op

        if (debug & 2)
        {
            InfoInFunction
                << "Perform &= on self: ignore" << nl;
        }

        return *this;
    }
    else if (none())
    {
        // no-op: nothing is set - no intersection possible
        return *this;
    }
    else if (other.none())
    {
        // no-op: other has nothing set - no intersection possible
        reset();
        return *this;
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

    return *this;
}


Foam::bitSet& Foam::bitSet::orEq(const bitSet& other)
{
    if (&other == this)
    {
        // Self '|=' : no-op

        if (debug & 2)
        {
            InfoInFunction
                << "Perform |= on self: ignore" << nl;
        }

        return *this;
    }
    else if (other.none())
    {
        // no-op: nothing can change
        return *this;
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

    return *this;
}


Foam::bitSet& Foam::bitSet::xorEq(const bitSet& other)
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
        return *this;
    }
    else if (other.none())
    {
        // no-op: nothing can change
        return *this;
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

    return *this;
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


Foam::bitSet::bitSet(const bitSet& bitset, const labelRange& range)
:
    bitSet(range.size())
{
    label pos = range.start();
    const label len = range.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, bitset.get(pos));
        ++pos;
    }
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
    fill(false);
    resize(bools.size());

    unsigned bitIdx = 0u;
    auto* packed = blocks_.data();

    // Set according to indices that are true
    for (const auto b : bools)
    {
        if (b)
        {
            *packed |= (1u << bitIdx);
        }

        if (++bitIdx >= PackedList<1>::elem_per_block)
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


void Foam::bitSet::set(const labelRange& range)
{
    labelRange slice(range);
    slice.adjust();  // No negative start, size adjusted accordingly

    // Range is invalid (zero-sized or entirely negative) - noop
    if (slice.empty())
    {
        return;
    }

    // Range finishes at or beyond the right side.
    // - zero fill any gaps that we might create.
    // - flood-fill the reset, which now corresponds to the full range.
    if (slice.end_value() >= size())
    {
        reserve(slice.end_value());
        resize(slice.begin_value(), false);
        resize(slice.end_value(), true);
        return;
    }

    // The more difficult case - everything in between.
    // 1. sequence may begin/end in the same block
    // 2. Cover more than one block
    //    a. with partial coverage in the first block
    //    b. with partial coverage in the end block

    // The begin block/offset
    unsigned int bblock = slice.begin_value() / elem_per_block;
    unsigned int bmask  = slice.begin_value() % elem_per_block;

    // The end block/offset
    unsigned int eblock = slice.end_value() / elem_per_block;
    unsigned int emask  = slice.end_value() % elem_per_block;

    // Transform offsets to lower bit masks
    if (bmask) bmask = mask_lower(bmask);
    if (emask) emask = mask_lower(emask);

    if (bblock == eblock)
    {
        // Same block - flll between the begin/end bits.
        // Example:
        // bmask = 0000000000001111  (lower bits)
        // emask = 0000111111111111  (lower bits)
        // -> set  0000111111110000  (xor)

        blocks_[bblock] |= (emask^bmask);
    }
    else
    {
        if (bmask)
        {
            // The first (partial) block
            // - set everything above the bmask.
            blocks_[bblock] |= (~bmask);
            ++bblock;
        }

        // Fill these blocks
        for (unsigned blocki = bblock; blocki < eblock; ++blocki)
        {
            blocks_[blocki] = (~0u);
        }

        if (emask)
        {
            // The last (partial) block.
            // - set everything below emask.
            blocks_[eblock] |= (emask);
        }
    }
}


void Foam::bitSet::unset(const labelRange& range)
{
    // Require intersection with the current bitset
    const labelRange slice = range.subset0(size());

    // Range does not intersect (invalid, empty, bitset is empty)
    if (slice.empty())
    {
        return;
    }

    // Range finishes at or beyond the right side.
    if (slice.end_value() >= size())
    {
        // The original size
        const label orig = size();

        resize(slice.begin_value(), false);
        resize(orig, false);
        return;
    }


    // The more difficult case - everything in between.
    // 1. sequence may begin/end in the same block
    // 2. Cover more than one block
    //    a. with partial coverage in the first block
    //    b. with partial coverage in the end block

    // The begin block/offset
    unsigned int bblock = slice.begin_value() / elem_per_block;
    unsigned int bmask  = slice.begin_value() % elem_per_block;

    // The end block/offset
    unsigned int eblock = slice.end_value() / elem_per_block;
    unsigned int emask  = slice.end_value() % elem_per_block;

    // Transform offsets to lower bit masks
    if (bmask) bmask = mask_lower(bmask);
    if (emask) emask = mask_lower(emask);

    if (bblock == eblock)
    {
        // Same block - flll between the begin/end bits.
        // Example:
        // bmask = 0000000000001111  (lower bits)
        // emask = 0000111111111111  (lower bits)
        // -> set  0000111111110000  (xor)
        // -> ~    1111000000001111

        blocks_[bblock] &= (~(emask^bmask));
    }
    else
    {
        if (bmask)
        {
            // The first (partial) block
            // - only retain things below bmask.
            blocks_[bblock] &= (bmask);
            ++bblock;
        }

        // Clear these blocks
        for (unsigned blocki = bblock; blocki < eblock; ++blocki)
        {
            blocks_[blocki] = (0u);
        }

        if (emask)
        {
            // The last (partial) block.
            // - only retain things above bmask.
            blocks_[eblock] &= (~emask);
        }
    }
}


Foam::labelList Foam::bitSet::toc() const
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
        unsigned int blockval = blocks_[blocki];

        if (blockval)
        {
            for (label pos = (blocki * elem_per_block); blockval; ++pos)
            {
                if (blockval & 1u)
                {
                    output[nItem] = pos;
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
