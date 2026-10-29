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

Class
    Foam::CSRAddressing

Description
    An offsets addressing table for compact sparse row (CSR) addressing.
    CSR-style addressing is used in compact list-of-list data,
    global mesh addressing etc.  The CSRAddressing is templated
    on the integral type used for the offsets. In many places using
    an \c int64 for the offsets may be necessary even when the rest of
    the code uses \c int32 for its sizes.

    The offset table is the size of the number of rows+1
    whose elements are the
    accumulated sizes of the rows, i.e.
      - offset[i] gives the index of first element of row i
      - offset[i+1] - offset[i] is the number of elements in row i

    Note that an empty CSRAddressing should have empty offsets.

    Regardless of the integral type used for the offsets themselves,
    the row indices always use \c Foam::label since these correspond
    to \em sub-lists etc which are restricted to the normal addressing
    constraints.

SourceFiles
    CSRAddressingI.H
    CSRAddressing.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_CSRAddressing_H
#define Foam_CSRAddressing_H

#include "List.H"
#include "labelRange.H"
#include <algorithm>
#include <limits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Int> class CSRAddressing;
template<class Int> Istream& operator>>(Istream&, CSRAddressing<Int>&);
template<class Int> Ostream& operator<<(Ostream&, const CSRAddressing<Int>&);

/*---------------------------------------------------------------------------*\
                        Class CSRAddressing Declaration
\*---------------------------------------------------------------------------*/

template<class IntType>
class CSRAddressing
{
    // Private Data

        //- Offset (addressing) table
        List<IntType> offsets_;


    // Private Member Functions

        //- Construct from the sizes of list of lists
        template<class ListListType>
        static CSRAddressing<IntType> pack_impl
        (
            const ListListType& lists,
            const bool checkOverflow = false
        );


protected:

    // Protected Member Functions

        //- Report overflow at specified (non-negative) index
        static void reportOverflowAndExit
        (
            label idx,
            int64_t prevOffset = -1,  // The last valid offset value
            int64_t count = 0         // The count to add to prevOffset
        );

        //- Check if adding count to prevOffset would exceed the maximum
        //- value of T (tested before the addition, to avoid signed overflow)
        template<class T>
        static bool isOverflow(T prevOffset, T count) noexcept
        {
            return
            (
                (count > 0)
             && (prevOffset > (std::numeric_limits<T>::max() - count))
            );
        }

public:

    // Public Data Types

        //- The integer type used for offsets
        typedef IntType value_type;


    // Static Member Functions

        //- Return a null CSRAddressing (reference to a nullObject).
        static const CSRAddressing<IntType>& null() noexcept
        {
            return NullObjectRef<CSRAddressing<IntType>>();
        }

        //- Construct from the sizes of list of lists
        template<class SubListType>
        static CSRAddressing<IntType> pack
        (
            const UList<SubListType>& lists,
            bool checkOverflow = false
        );

        //- Construct from the sizes of an indirect list of lists
        template<class SubListType, class Addr>
        static CSRAddressing<IntType> pack
        (
            const IndirectListBase<SubListType, Addr>& lists,
            bool checkOverflow = false
        );


    // Generated Methods

        //- Copy construct
        CSRAddressing(const CSRAddressing&) = default;

        //- Move construct
        CSRAddressing(CSRAddressing&&) = default;

        //- Copy assignment
        CSRAddressing& operator=(const CSRAddressing&) = default;

        //- Move assignment
        CSRAddressing& operator=(CSRAddressing&&) = default;


    // Constructors

        //- Default construct (empty)
        inline constexpr CSRAddressing() noexcept;

        //- Construct from number of rows
        explicit inline CSRAddressing(label numRows);

        //- Copy construct from offsets.
        //- Discards/ignores degenerate offsets (with size less than 2).
        explicit inline CSRAddressing
        (
            Foam::direct_construct_t,
            const UList<IntType>& rowOffsets
        );

        //- Move construct from offsets.
        //- Discards/ignores degenerate offsets (with size less than 2).
        explicit inline CSRAddressing
        (
            Foam::direct_construct_t,
            List<IntType>&& rowOffsets
        );

        //- Construct from row-sizes, with optional overflow check
        explicit inline CSRAddressing
        (
            const labelUList& rowSizes,
            bool checkOverflow = false
        );

        //- Construct from an indirect list of row lists,
        //- with optional overflow check
        template<class Addr>
        explicit CSRAddressing
        (
            const IndirectListBase<label, Addr>& rowSizes,
            bool checkOverflow = false
        );

        //- Construct as copy or re-use as specified.
        inline CSRAddressing(CSRAddressing& addr, bool reuse);

        //- Construct from Istream. Reads the offsets
        explicit inline CSRAddressing(Istream& is);


    // Member Functions

    // Access

        //- True if the number of rows is zero, or has totalSize == 0
        inline bool empty() const noexcept;

        //- The primary size (the number of rows)
        inline label length() const noexcept;

        //- The primary size (the number of rows)
        inline label size() const noexcept;

        //- True if single-row content
        inline bool single() const noexcept;

        //- The total addressed size
        inline IntType totalSize() const noexcept;

        //- The value corresponding to the first offset
        inline IntType begin_value() const noexcept;

        //- The value corresponding to the final offset,
        //- which is 1 beyond the end of the span
        inline IntType end_value() const noexcept;

        //- The offset table (1 larger than length)
        const List<IntType>& offsets() const noexcept { return offsets_; }

        //- Return non-const access to the offset table
        List<IntType>& offsets() noexcept { return offsets_; }


    // Queries

        //- True if the addressing span includes the value \p i.
        inline bool contains(IntType i) const noexcept;

        //- True if the addressing for \p row includes the value \p i.
        inline bool contains(label row, IntType i) const noexcept;

        //- The local row sizes. Same as localSizes
        inline labelList sizes() const;

        //- The local row starts (a slice into offsets)
        inline const UList<IntType> localStarts() const;

        //- The local row sizes
        labelList localSizes() const;

        //- Start offset for given row
        inline IntType localStart(label row) const;

        //- End offset for given row
        inline IntType localEnd(label row) const;

        //- Size of given row
        inline IntType localSize(label row) const;

        //- The max of local sizes, excluding the specified row
        label maxNonLocalSize(label row) const noexcept;

        //- Return start/size ranges
        List<IntRange<IntType>> ranges() const;

        //- Return start/size range for given row, or an empty
        //- range when out-of-bounds.
        inline IntRange<IntType> range(label row) const noexcept;

        //- The max row length used
        inline label maxSize() const noexcept;


    // Edit

        //- Clear addressing
        inline void clear();

        //- Set offsets as a single entry content (eg, local-only content)
        //- with a beginning offset of 0.
        inline void single(IntType len);

        //- Set offsets for a single entry with specified beginning offset,
        //- optionally with a total size (for trailing information).
        inline void single(IntType beg, IntType len, IntType total = 0);

        //- Set as special purpose \em sentinel value corresponding to
        //- a total-size. The corresponding offsets are \em "degenerate"
        //- since there will only be a single offset instead of an offset
        //- pair. However, the totalSize() method can still be used.
        inline void sentinel(IntType total);

        //- Change length of CSRAddressing
        inline void resize(label numRows);

        //- Reset row-sizes, with optional overflow check
        void reset(const labelUList& rowSizes, bool checkOverflow = false);

        //- Alter local addressing size for given row
        void setLocalSize(label row, label len);


    // Other

        //- From row-local index to the global index (linear addressing)
        inline IntType toGlobal(label row, IntType i) const;

        //- From global (linear addressing) to local index on \p row
        inline IntType toLocal(label row, IntType i) const;

        //- Find the row corresponding to the linear location \p i,
        //- checking the \p row first, followed by a binary search.
        //  \return the row number or -1 (not found)
        inline label findRow(label row, IntType i) const;

        //- Find the row above \p row corresponding to linear location \p i.
        //  \return the row number or -1 (not found)
        inline label findRowAbove(label row, IntType i) const;

        //- Find the row below \p row  corresponding to linear location \p i.
        //  \return the row number or -1 (not found)
        inline label findRowBelow(label row, IntType i) const;

        //- Find the row corresponding to the addressing for \p i.
        //  \return the row number or -1 (not found)
        inline label findRow(IntType idx) const;


    // Member Operators

        //- Compare for equality - uses the offsets
        bool operator==(const CSRAddressing& rhs) const
        {
            return (this->offsets() == rhs.offsets());
        }

        //- Compare for inequality - uses the offsets
        bool operator!=(const CSRAddressing& rhs) const
        {
            return !(*this == rhs);
        }

        //- Compare for less-than - uses the offsets
        bool operator<(const CSRAddressing& rhs) const
        {
            return (this->offsets() < rhs.offsets());
        }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read offsets addressing from Istream
template<class IntType>
Istream& operator>>(Istream& is, CSRAddressing<IntType>& addr)
{
    return (is >> addr.offsets());
}

//- Write offsets addressing to Ostream
template<class IntType>
Ostream& operator<<(Ostream& os, const CSRAddressing<IntType>& addr)
{
    return (os << addr.offsets());
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CSRAddressingI.H"
#include "CSRAddressing.txx"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
