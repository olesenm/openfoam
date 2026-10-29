/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

template<class IntType>
inline constexpr Foam::CSRAddressing<IntType>::CSRAddressing() noexcept
:
    offsets_()
{}


template<class IntType>
inline Foam::CSRAddressing<IntType>::CSRAddressing(label numRows)
:
    offsets_()
{
    if (numRows > 0)
    {
        // Sanity
        offsets_.resize(numRows+1, 0);
    }
}


template<class IntType>
inline Foam::CSRAddressing<IntType>::CSRAddressing
(
    Foam::direct_construct_t,
    const UList<IntType>& rowOffsets
)
:
    offsets_()
{
    if (rowOffsets.size() > 1)  // Enforce sizing sanity
    {
        offsets_ = rowOffsets;
    }
}


template<class IntType>
inline Foam::CSRAddressing<IntType>::CSRAddressing
(
    Foam::direct_construct_t,
    List<IntType>&& rowOffsets
)
:
    offsets_(std::move(rowOffsets))
{
    if (offsets_.size() == 1)  // Enforce sizing sanity
    {
        offsets_.clear();
    }
}


template<class IntType>
inline Foam::CSRAddressing<IntType>::CSRAddressing
(
    const labelUList& rowSizes,
    bool checkOverflow
)
:
    offsets_()
{
    reset(rowSizes, checkOverflow);
}


template<class IntType>
inline Foam::CSRAddressing<IntType>::CSRAddressing
(
    CSRAddressing<IntType>& addr,
    bool reuse
)
:
    offsets_(addr.offsets_, reuse)
{}


template<class IntType>
inline Foam::CSRAddressing<IntType>::CSRAddressing
(
    Istream& is
)
:
    offsets_(is)
{}


//
// Constructor-like
//

template<class IntType>
inline void Foam::CSRAddressing<IntType>::single(IntType len)
{
    offsets_.resize_fill(2, 0);
    if (len > 0)
    {
        offsets_[1] = len;
    }
}


template<class IntType>
inline void Foam::CSRAddressing<IntType>::single
(
    IntType beg,
    IntType len,
    IntType total
)
{
    const IntType ending = (beg + len);

    if (total <= ending)
    {
        offsets_.resize_nocopy(2);
        offsets_[0] = beg;
        offsets_[1] = ending;
    }
    else
    {
        offsets_.resize_nocopy(3);
        offsets_[0] = beg;
        offsets_[1] = ending;
        offsets_[2] = total;
    }
}


template<class IntType>
inline void Foam::CSRAddressing<IntType>::sentinel(IntType total)
{
    // Note this creates a degenerate offset table with only 1 entry
    // instead of the usual minimum of 2 but does mean that we can
    // evaluate totalSize() without per-rank information at all.
    offsets_.resize_fill(1, std::max(IntType(0), total));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
inline bool Foam::CSRAddressing<IntType>::empty() const noexcept
{
    // Filters out size=1 offsets
    const label len = (offsets_.size() - 1);
    return (len < 1) || (*(offsets_.cdata() + len) == 0);
}


template<class IntType>
inline bool Foam::CSRAddressing<IntType>::single() const noexcept
{
    return (offsets_.size() == 2);
}


template<class IntType>
inline Foam::label Foam::CSRAddressing<IntType>::length() const noexcept
{
    // Filters out size=1 offsets
    const label len = (offsets_.size() - 1);
    return (len < 1) ? 0 : len;
}


template<class IntType>
inline Foam::label Foam::CSRAddressing<IntType>::size() const noexcept
{
    // Filters out size=1 offsets
    const label len = (offsets_.size() - 1);
    return (len < 1) ? 0 : len;
}


template<class IntType>
inline IntType Foam::CSRAddressing<IntType>::begin_value() const noexcept
{
    // No size=1 filtering
    return (offsets_.empty() ? 0 : *(offsets_.cdata()));
}


template<class IntType>
inline IntType Foam::CSRAddressing<IntType>::end_value() const noexcept
{
    // No size=1 filtering : same as totalSize()
    return (offsets_.empty() ? 0 : *(offsets_.cdata() + (offsets_.size()-1)));
}


template<class IntType>
inline IntType
Foam::CSRAddressing<IntType>::totalSize() const noexcept
{
    // No size=1 filtering : same as end_value()
    return (offsets_.empty() ? 0 : *(offsets_.cdata() + (offsets_.size()-1)));
}


template<class IntType>
inline bool Foam::CSRAddressing<IntType>::contains(IntType i) const noexcept
{
    return (i >= begin_value() && i < end_value());
}


template<class IntType>
inline bool Foam::CSRAddressing<IntType>::contains
(
    label row,
    IntType i
) const noexcept
{
    const auto* data_ptr = offsets_.cdata();

    return
    (
        (row >= 0 && row < offsets_.size()-1)
     && (data_ptr[row] <= i && i < data_ptr[row+1])
    );
}


template<class IntType>
inline Foam::labelList
Foam::CSRAddressing<IntType>::sizes() const
{
    return localSizes();
}


template<class IntType>
inline const Foam::UList<IntType>
Foam::CSRAddressing<IntType>::localStarts() const
{
    const auto len = (offsets_.size() - 1);

    // Filters out size=1 offsets
    if (len < 1) return UList<IntType>::null();

    return offsets_.slice(0, len);
}


template<class IntType>
inline IntType Foam::CSRAddressing<IntType>::localStart(label row) const
{
    // No bounds checking
    return offsets_[row];
}


template<class IntType>
inline IntType Foam::CSRAddressing<IntType>::localSize(label row) const
{
    // No bounds checking
    return (offsets_[row+1] - offsets_[row]);
}


template<class IntType>
inline IntType Foam::CSRAddressing<IntType>::localEnd(label row) const
{
    // No bounds checking
    return offsets_[row+1];
}


template<class IntType>
inline Foam::label Foam::CSRAddressing<IntType>::maxSize() const noexcept
{
    return maxNonLocalSize(-1);
}


template<class IntType>
inline Foam::IntRange<IntType>
Foam::CSRAddressing<IntType>::range(label row) const noexcept
{
    // Filters out size=1 offsets
    if (row >= 0 && row < (offsets_.size() - 1))
    {
        const auto* data_ptr = offsets_.cdata();

        return IntRange<IntType>
        (
            data_ptr[row],
            (data_ptr[row+1] - data_ptr[row])
        );
    }
    else
    {
        // Out of bounds
        return IntRange<IntType>();
    }
}


template<class IntType>
inline IntType Foam::CSRAddressing<IntType>::toGlobal
(
    label row,
    IntType i
) const
{
    // No bounds checking
    return i + offsets_[row];
}


template<class IntType>
inline IntType Foam::CSRAddressing<IntType>::toLocal
(
    label row,
    IntType i
) const
{
    if (!contains(row, i))
    {
        FatalErrorInFunction
            << "Index " << i << " does not belong on rank/row " << row
            << "\nOffsets:" << offsets_
            << Foam::abort(FatalError);
    }

    // Already bounds checked
    return (i - offsets_[row]);
}


template<class IntType>
inline Foam::label
Foam::CSRAddressing<IntType>::findRow(IntType i) const
{
    if (i < begin_value() || i >= end_value())  //<- !contains(i)
    {
        return -1;
    }

    return static_cast<label>
    (
        std::upper_bound(offsets_.cbegin(), offsets_.cend(), i)
      - offsets_.cbegin() - 1  //< uses (-1) for start instead of end index
    );
}


template<class IntType>
inline Foam::label
Foam::CSRAddressing<IntType>::findRow(label row, IntType i) const
{
    if (row < 0 || row >= offsets_.size()-1)
    {
        // Invalid row, use regular routine
        return findRow(i);
    }

    if (i < begin_value() || i >= end_value())  //<- !contains(i)
    {
        // Out of range
        return -1;
    }

    // Handle (frequent) queries within the given row,
    // or make it the split point for restricting binary searches

    if (offsets_[row] <= i && i < offsets_[row+1])
    {
        // is contained within row
        return row;
    }

    if (i < offsets_[row])
    {
        // Restrict search to rows below
        const auto first = offsets_.cbegin();

        return static_cast<label>
        (
            std::upper_bound(first, first+row, i)
          - first - 1  //< uses (-1) for start instead of end index
        );
    }
    else
    {
        // Restrict search to rows above
        const auto first = (offsets_.cbegin() + row+1);

        return
        (
            (row+1)
          + static_cast<label>
            (
                std::upper_bound(first, offsets_.end(), i)
              - first - 1  //< uses (-1) for start instead of end index
            )
        );
    }
}


template<class IntType>
inline Foam::label
Foam::CSRAddressing<IntType>::findRowBelow(label row, IntType i) const
{
    // Simple checks first
    if
    (
        (row <= 0)                  // Nothing below first row
     || (row >= offsets_.size())    // Invalid row
     || (i >= offsets_[row])        // Exclude current row and above
     || (i < begin_value())         // Beyond the begin
    )
    {
        return -1;
    }

    {
        // Restrict search to rows below
        const auto first = offsets_.cbegin();

        return static_cast<label>
        (
            std::upper_bound(first, first+row, i)
          - first - 1  //< uses (-1) for start instead of end index
        );
    }
}


template<class IntType>
inline Foam::label
Foam::CSRAddressing<IntType>::findRowAbove(label row, IntType i) const
{
    // Simple checks first
    if
    (
        (row < 0)                   // Invalid row
     || (row >= offsets_.size()-1)  // Nothing above the last row
     || (i < offsets_[row+1])       // Exclude current row and below
     || (i >= end_value())          // Beyond the end
    )
    {
        return -1;
    }

    {
        // Restrict search to rows above
        const auto first = (offsets_.cbegin() + row+1);

        return
        (
            (row+1)
          + static_cast<label>
            (
                std::upper_bound(first, offsets_.end(), i)
              - first - 1  //< uses (-1) for start instead of end index
            )
        );
    }
}


template<class IntType>
inline void Foam::CSRAddressing<IntType>::clear()
{
    offsets_.clear();
}


template<class IntType>
inline void Foam::CSRAddressing<IntType>::resize(label numRows)
{
    if (numRows <= 0)
    {
        // Clear
        offsets_.clear();
    }
    else if (numRows < size())
    {
        // Shrink
        offsets_.resize(numRows+1);
    }
    else if (numRows > size())
    {
        // Extend number of rows, each with local size of 0
        IntType endOffset(0);

        if (!offsets_.empty())
        {
            endOffset = offsets_.back();
        }

        offsets_.resize(numRows+1, endOffset);
    }
}


// ************************************************************************* //
