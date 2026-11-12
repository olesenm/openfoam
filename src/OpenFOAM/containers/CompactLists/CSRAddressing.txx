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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class IntType>
void Foam::CSRAddressing<IntType>::reportOverflowAndExit
(
    label idx,
    int64_t prevOffset,
    int64_t count
)
{
    if (idx < 0)
    {
        // No overflow tagged
        return;
    }

    FatalErrorInFunction
        << "Overflow : sum of sizes exceeds max for "
        << (8*sizeof(IntType)) << "-bit integral after index " << idx;

    if (prevOffset >= 0 && count >= 0)
    {
        FatalError
            << " while trying to add (" << count
            << ") to offset (" << prevOffset << ")";
    }

    FatalError
        << "\nPlease consider a larger integral datatype.\n"
        << Foam::exit(FatalError);
}


template<class IntType>
template<class ListListType>
Foam::CSRAddressing<IntType>
Foam::CSRAddressing<IntType>::pack_impl
(
    const ListListType& lists,
    bool checkOverflow
)
{
    CSRAddressing<IntType> addr;

    if (const label len = lists.size(); len > 0)
    {
        auto& starts = addr.offsets();

        starts.resize(len+1);

        IntType total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = static_cast<IntType>(lists[i].size());

            // No contribution to sizes
            if (count < 1)
            {
                starts[i] = total;
                continue;
            }

            // No overflow check for int64...
            if constexpr (sizeof(IntType) < sizeof(int64_t))
            {
                if (checkOverflow && isOverflow(total, count))
                {
                    reportOverflowAndExit(i, total, count);
                }
            }

            starts[i] = total;
            total += count;
        }
        starts[len] = total;
    }

    return addr;
}


template<class IntType>
template<class SubListType>
Foam::CSRAddressing<IntType>
Foam::CSRAddressing<IntType>::pack
(
    const UList<SubListType>& lists,
    bool checkOverflow
)
{
    return
    (
        CSRAddressing<IntType>::
        pack_impl<UList<SubListType>>
        (
            lists,
            checkOverflow
        )
    );
}


template<class IntType>
template<class SubListType, class Addr>
Foam::CSRAddressing<IntType>
Foam::CSRAddressing<IntType>::pack
(
    const IndirectListBase<SubListType, Addr>& lists,
    bool checkOverflow
)
{
    return
    (
        CSRAddressing<IntType>::
        pack_impl<IndirectListBase<SubListType, Addr>>
        (
            lists,
            checkOverflow
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IntType>
template<class Addr>
Foam::CSRAddressing<IntType>::CSRAddressing
(
    const IndirectListBase<label, Addr>& rowSizes,
    bool checkOverflow
)
:
    offsets_()
{
    if (const label len = rowSizes.size(); len > 0)
    {
        offsets_.resize(len+1);
        auto& starts = offsets_;

        IntType total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = static_cast<IntType>(rowSizes[i]);

            // Ignore negative values as possible placeholders,
            // and 0 has no effect on total anyhow.
            if (count < 1)
            {
                starts[i] = total;
                continue;
            }

            if constexpr (sizeof(IntType) < sizeof(int64_t))
            {
                if (checkOverflow && isOverflow(total, count))
                {
                    reportOverflowAndExit(i, total, count);
                }
            }

            starts[i] = total;
            total += count;
        }
        starts[len] = total;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
Foam::List<Foam::IntRange<IntType>>
Foam::CSRAddressing<IntType>::ranges() const
{
    List<IntRange<IntType>> values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label i = 0; i < len; ++i)
    {
        values[i].reset(offsets_[i], offsets_[i+1] - offsets_[i]);
    }

    return values;
}


template<class IntType>
void Foam::CSRAddressing<IntType>::reset
(
    const labelUList& rowSizes,
    [[maybe_unused]] bool checkOverflow
)
{
    if (const label len = rowSizes.size(); len > 0)
    {
        offsets_.resize_nocopy(len+1);
        auto& starts = offsets_;

        IntType total = 0;
        for (label i = 0; i < len; ++i)
        {
            const auto count = static_cast<IntType>(rowSizes[i]);

            // Ignore negative values as possible placeholders,
            // and 0 has no effect on total anyhow.
            if (count < 1)
            {
                starts[i] = total;
                continue;
            }

            // No overflow check for int64...
            if constexpr (sizeof(IntType) < sizeof(int64_t))
            {
                if (checkOverflow && isOverflow(total, count))
                {
                    reportOverflowAndExit(i, total, count);
                }
            }

            starts[i] = total;
            total += count;
        }
        starts[len] = total;
    }
    else
    {
        offsets_.clear();
    }
}


template<class IntType>
void Foam::CSRAddressing<IntType>::setLocalSize
(
    label row,
    label len
)
{
    if
    (
        const label numRows = this->length();
        (row < 0 || numRows <= row) || (len < 0)
    )
    {
        FatalErrorInFunction
            << "Attempting to set rank/row:" << row << " of " << numRows
            << " to length=" << len
            << Foam::abort(FatalError);
    }
    else  // (len >= 0)
    {
        const IntType delta = (len - (offsets_[row+1] - offsets_[row]));

        // TBD: additional overflow check
        if (delta)
        {
            for (label i = row+1; i < offsets_.size(); ++i)
            {
                offsets_[i] += delta;
            }
        }
    }
}


template<class IntType>
Foam::labelList Foam::CSRAddressing<IntType>::localSizes() const
{
    labelList values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label i = 0; i < len; ++i)
    {
        values[i] = static_cast<label>(offsets_[i+1] - offsets_[i]);
    }

    return values;
}


template<class IntType>
Foam::label Foam::CSRAddressing<IntType>::maxNonLocalSize
(
    label row
) const noexcept
{
    IntType maxLen = 0;

    const auto* data_ptr = offsets_.cdata();

    // Filters out size=1 offsets
    const label len = (offsets_.size() - 1);

    for (label i = 0; i < len; ++i)
    {
        if (i != row)
        {
            IntType local = (data_ptr[i+1] - data_ptr[i]);
            maxLen = std::max(maxLen, local);
        }
    }

    return static_cast<label>(maxLen);
}


// ************************************************************************* //
