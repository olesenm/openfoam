/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
    Copyright (C) 2026 Keysight Technologies
    Copyright (C) 2026 Mark Olesen
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
    Foam::globalIndex

Description
    Calculates a non-overlapping list of offsets based on an input size
    (eg, number of cells) from different MPI ranks.
    These offsets can be used to define a unique (global) indexing
    number.

    For example,
    \verbatim
        globalIndex globalFaces(mesh.nFaces());
        label globalFacei = globalFaces.toGlobal(facei);
    \endverbatim

Note
    The globalIndex is generally used when the ranks need information
    about each other's position, whereas GlobalOffset and OffsetRange
    are used when only the local addressing information (and total Size)
    is needed without additional information about other ranks.

Warning
    The unique integer is a Foam::label, so may be subject to overflow
    problems (2G max) with a 32-bit label compilation.

SourceFiles
    globalIndex.cxx
    globalIndex.txx
    globalIndexI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_globalIndex_H
#define Foam_globalIndex_H

#include "Pstream.H"
#include "CSRAddressing.H"
#include "CompactListList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class IntType> class GlobalIndex;

/*---------------------------------------------------------------------------*\
                       Class globalIndexBase Declaration
\*---------------------------------------------------------------------------*/

//- Template invariant parts for globalIndex etc.
class globalIndexBase
{
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

        //- My processor number (on worldComm)
        static int myProcNumber()
        {
            return UPstream::myProcNo(UPstream::worldComm);
        }

public:

    // Public Data Types

        //- Dispatch tag:
        //- Construct with a single (local size) entry, no communication.
        struct gatherNone{};

        //- Dispatch tag:
        //- Construct 'one-sided' by gathering the local sizes.
        //  The master will have all content, the sub-ranks will have none.
        struct gatherOnly{};

        //- Dispatch tag:
        //- An augmented version of gatherOnly.
        //- Construct 'one-sided' by gathering the local sizes but also
        //- distribute the total-size to each rank.
        //  The master will have all content, sub-ranks will have
        //  \em single-value content (offsets containing only the total size).
        struct gatherTotal{};

        //- Dispatch tag:
        //- An augmented version of gatherOnly.
        //- Construct 'one-sided' by gathering the local sizes but also
        //- distribute the start/total to each rank.
        //  The master will have all content, sub-ranks will have
        //  \em "localized" offsets of either three values (beg, end, total)
        //  or two values (begin, end) if there is no trailing content.
        struct gatherOffset{};


    // Properties

        //- Default communication type for operations (eg, gather/scatter)
        static constexpr UPstream::commsTypes defaultCommsType() noexcept
        {
            return UPstream::commsTypes::nonBlocking;
        }


    // Helper Functions

        //- Based on the local input size(s),
        //- calculate a globally-consistent local start offset.
        //  No overflow checks.
        //  \note Restricted to integral types.
        template<class IntType>
        static IntType calcOffset
        (
            IntType localSize,
            int communicator = UPstream::worldComm  //!< communicator
        )
        {
            return UPstream::find_offset(localSize, communicator);
        }

        //- Based on the local input size(s),
        //- calculate a globally-consistent local start offset.
        //  With optional overflow checking.
        //  \note Restricted to integral types.
        template<class IntType>
        static IntType calcOffset
        (
            IntType localSize,
            int communicator,  //!< communicator
            [[maybe_unused]] bool checkOverflow
        )
        {
            IntType start = UPstream::find_offset(localSize, communicator);

            // Assume no overflow for large integrals
            if constexpr (sizeof(IntType) < sizeof(int64_t))
            {
                if
                (
                    checkOverflow && UPstream::is_parallel(communicator)
                 && globalIndexBase::isOverflow(start, localSize)
                )
                {
                    reportOverflowAndExit
                    (
                        UPstream::myProcNo(communicator), start, localSize
                    );
                }
            }
            return start;
        }

        //- Based on the local input size(s),
        //- calculate a globally-consistent range (offset/size).
        //  No overflow checks.
        //  \note Restricted to integral types.
        template<class IntType>
        std::enable_if_t<std::is_integral_v<IntType>, IntRange<IntType>>
        static calcRange
        (
            IntType localSize,
            int communicator = UPstream::worldComm  //!< communicator
        )
        {
            IntType start = UPstream::find_offset(localSize, communicator);
            return IntRange<IntType>(start, localSize);
        }

        //- Based on the local input size(s),
        //- calculate a globally-consistent range (offset/size).
        //- With optional overflow checking.
        //  \note Restricted to integral types.
        template<class IntType>
        std::enable_if_t<std::is_integral_v<IntType>, IntRange<IntType>>
        static calcRange
        (
            IntType localSize,
            int communicator,  //!< communicator
            bool checkOverflow
        )
        {
            IntType start =
                globalIndexBase::calcOffset
                (
                    localSize,
                    communicator,
                    checkOverflow
                );

            return IntRange<IntType>(start, localSize);
        }

        //- Get the receive sizes.
        //- On master the length == nProcs and element [0] corresponds
        //- to the max of all other receive sizes.
        //- On non-master, the list is zero-sized.
        //  \note Restricted to integral types.
        template<class IntType>
        static List<IntType> calcRecvSizes
        (
            IntType localSize,
            int communicator = UPstream::worldComm  //!< communicator
        );

        //- Calculate offsets from a list of local sizes,
        //- with optional check for int32 (label) overflow
        //  The OffsetType defaults to the input IntType (when void),
        //  e.g. calcOffsets<int64_t>(const UList<int32_t>&)
        template<class OffsetType = void, class IntType>
        static List
        <
            std::conditional_t
            <std::is_void_v<OffsetType>, IntType, OffsetType>
        >
        calcOffsets
        (
            const UList<IntType>& counts,
            bool checkOverflow = false
        );

        //- Calculate offsets from an indirect list of local sizes,
        //- with optional check for int32 (label) overflow.
        //  The OffsetType defaults to the input IntType (when void)
        template<class OffsetType = void, class IntType, class Addr>
        static List
        <
            std::conditional_t
            <std::is_void_v<OffsetType>, IntType, OffsetType>
        >
        calcOffsets
        (
            const IndirectListBase<IntType, Addr>& counts,
            bool checkOverflow = false
        );

        //- Calculate offsets from list of lists,
        //- with optional check for int32 (label) overflow.
        //  The OffsetType defaults to the input IntType (when void)
        template<class OffsetType = void, class SubListType>
        static List
        <
            std::conditional_t
            <std::is_void_v<OffsetType>, Foam::label, OffsetType>
        >
        calcListOffsets
        (
            const UList<SubListType>& lists,
            bool checkOverflow = false
        );
};


/*---------------------------------------------------------------------------*\
                         Class GlobalIndex Declaration
\*---------------------------------------------------------------------------*/

template<class IntType>
class GlobalIndex
:
    public globalIndexBase,
    public CSRAddressing<IntType>
{
    // Private Data Types

        //- The template instance for the addressing
        typedef CSRAddressing<IntType> addr_type;


    // Private Member Functions

        //- Reset offsets for a non-parallel communicator:
        //- localSize for the current rank, zero elsewhere.
        //- Empty if the communicator has no ranks.
        void reset_nonParallel(label localSize, int communicator);


protected:

    // Protected Member Functions

        // Multiple inheritance...
        using globalIndexBase::reportOverflowAndExit;

        // Multiple inheritance...
        using globalIndexBase::isOverflow;


public:

    //- Require integral type.
    //  This should normally be at least as large as <label>, but
    //  \c static_assert(sizeof(IntType) >= sizeof(Foam::label), ...)
    //  would prevent generalized templating.
    static_assert(std::is_integral_v<IntType>, "Integral required");

    // STL type definitions

        //- Type of values contained
        typedef IntType value_type;

        //- The type that can represent the size
        typedef IntType size_type;


    // Generated Methods

        //- Default construct (empty)
        GlobalIndex() noexcept = default;

        //- Copy construct
        GlobalIndex(const GlobalIndex&) = default;

        //- Move construct
        GlobalIndex(GlobalIndex&&) = default;

        //- Copy assignment
        GlobalIndex& operator=(const GlobalIndex&) = default;

        //- Move assignment
        GlobalIndex& operator=(GlobalIndex&&) = default;


    // Constructors

        //- Copy construct from a list of offsets.
        //- No communication required
        inline explicit GlobalIndex(const UList<IntType>& listOffsets);

        //- Move construct from a list of offsets.
        //- No communication required
        inline explicit GlobalIndex(List<IntType>&& listOffsets);

        //- Construct from a list of sizes and calculate the offsets.
        //- No communication required
        inline GlobalIndex
        (
            globalIndexBase::gatherNone,
            const labelUList& localSizes
        );

        //- Construct from local size, using gather/broadcast
        //- with default/specified communicator if parallel.
        inline explicit GlobalIndex
        (
            label localSize,
            int communicator = UPstream::worldComm,
            bool parallel = UPstream::parRun()  //!< use parallel comms
        );

        //- Construct with a single (local size) entry, no communication
        inline GlobalIndex
        (
            globalIndexBase::gatherNone,
            label localSize,
            int communicator = -1   //!< (unused) no communicator needed
        );

        //- Construct 'one-sided' from local sizes.
        //- The master will have all content, sub-ranks no content.
        //
        //  \note can be used when UPstream::parRun() is false.
        inline GlobalIndex
        (
            globalIndexBase::gatherOnly,
            label localSize,
            int communicator = UPstream::worldComm,  //!< communicator
            bool parallel = UPstream::parRun()       //!< use parallel comms
        );

        //- Construct an augmented 'one-sided' from local sizes
        //- and distribute the total-size to each rank.
        //- The master will have all content, sub-ranks will have
        //- \em single-value content (offsets containing only the total size).
        //
        //  \note can be used when UPstream::parRun() is false.
        inline GlobalIndex
        (
            globalIndexBase::gatherTotal,
            label localSize,
            int communicator = UPstream::worldComm,  //!< communicator
            bool parallel = UPstream::parRun()       //!< use parallel comms
        );

        //- Construct an augmented 'one-sided' from local sizes
        //- and distribute the start/total to each rank.
        //- The master will have all content, sub-ranks will have
        //- \em "localized" offsets of either three values (beg, end, total)
        //- or two values (begin, end) if there is no trailing content.
        //
        //  \note can be used when UPstream::parRun() is false.
        inline GlobalIndex
        (
            globalIndexBase::gatherOffset,
            label localSize,
            int communicator = UPstream::worldComm,  //!< communicator
            bool parallel = UPstream::parRun()       //!< use parallel comms
        );

        //- Construct from Istream.
        //- No communication required
        explicit GlobalIndex(Istream& is) : addr_type(is) {}


    // Member Functions

        //- Check for default constructed or total-size == 0
        bool empty() const noexcept { return addr_type::empty(); }

        //- The number of items covered by the offsets
        label length() const noexcept { return addr_type::length(); }

        //- Global sum of localSizes. Same as totalSize()
        FOAM_DEPRECATED_STRICT(2022-10, "totalSize() - unambiguous")
        IntType size() const noexcept { return addr_type::totalSize(); }

        //- The total addressed size, which corresponds to the
        //- end offset and also the sum of all localSizes.
        IntType totalSize() const noexcept { return addr_type::totalSize(); }

        //- The local sizes. Same as localSizes()
        labelList sizes() const { return addr_type::sizes(); }

        //- The local starts
        const UList<IntType> localStarts() const
        {
            return addr_type::localStarts();
        }

        //- The local sizes
        labelList localSizes() const { return addr_type::localSizes(); }

        //- Global max of localSizes
        IntType maxSize() const noexcept { return addr_type::maxSize(); }


    // Access

        //- Direct access to the offsets
        using addr_type::offsets;


    // Dimensions

        //- True if local-only content (ie, nProcs == 1).
        //- Such content is often created with gatherNone.
        bool single() const noexcept { return addr_type::single(); }

        //- The value corresponding to the first offset
        IntType begin_value() const noexcept
        {
            return addr_type::begin_value();
        }

        //- The value corresponding to the last offset (end offset),
        //- which is 1 beyond the end of the range.
        IntType end_value() const noexcept
        {
            return addr_type::end_value();
        }


    // Edit

        //- Reset to be empty (no offsets)
        void clear() { addr_type::clear(); }

        //- Set offsets as a single entry content (eg, local-only content)
        //- with a beginning offset of 0.
        void single(IntType len) { addr_type::single(len); }

        //- Set offsets for a single entry with specified beginning offset,
        //- optionally with a total size (for trailing information).
        void single(IntType beg, IntType len, IntType total = 0)
        {
            addr_type::single(beg, len, total);
        }

        //- Set special purpose \em sentinel for total-size only.
        //- The offsets are \em "degenerate" since there will only be one
        //- offset instead of at least two or more offsets.
        //- However, the totalSize() method is still valid.
        void sentinel(IntType total) { addr_type::sentinel(total); }

        //- Change the number of entries (nProcs) in the offsets table.
        //- Extending will fill with empty local sizes.
        void resize(label n) { addr_type::resize(n); }

        //- Alter local size for given processor
        void setLocalSize(label proci, IntType len)
        {
            addr_type::setLocalSize(proci, len);
        }

        //- Reset offsets from a list of local sizes,
        //- with optional check for overflow.
        //- No communication required
        void reset(const labelUList& counts, bool checkOverflow = false)
        {
            addr_type::reset(counts, checkOverflow);
        }

        //- Reset from local size, using gather/broadcast
        //- with default/specified communicator if parallel.
        void reset
        (
            label localSize,
            int communicator = UPstream::worldComm,  //!< communicator
            bool parallel = UPstream::parRun()       //!< use parallel comms
        );

        //- Reset to a single (local size) entry, no communication
        inline void reset
        (
            globalIndexBase::gatherNone,
            label localSize,
            int communicator = -1   //!< (unused) no communicator needed
        );

        //- Reset as 'one-sided' from local sizes.
        //- The master will have all content, sub-ranks no content.
        //
        //  \note can be used when UPstream::parRun() is false.
        void reset
        (
            globalIndexBase::gatherOnly,
            label localSize,
            int communicator = UPstream::worldComm,  //!< communicator
            bool parallel = UPstream::parRun()       //!< use parallel comms
        );

        //- Reset as augmented 'one-sided' from local sizes.
        //- The master will have all content, sub-ranks will have
        //- \em single-value content (offsets containing only the total size).
        //
        //  \note can be used when UPstream::parRun() is false.
        void reset
        (
            globalIndexBase::gatherTotal,
            label localSize,
            int communicator = UPstream::worldComm,  //!< communicator
            bool parallel = UPstream::parRun()       //!< use parallel comms
        );

        //- Reset as augmented 'one-sided' from local sizes.
        //- The master will have all content, sub-ranks will have
        //- \em "localized" offsets of either three values (beg, end, total)
        //- or two values (begin, end) if there is no trailing content.
        //
        //  \note can be used when UPstream::parRun() is false.
        void reset
        (
            globalIndexBase::gatherOffset,
            label localSize,
            int communicator = UPstream::worldComm,  //!< communicator
            bool parallel = UPstream::parRun()       //!< use parallel comms
        );

        //- Reset the globalIndex. Same as copy assignment.
        void reset(const GlobalIndex& rhs) { addr_type::operator=(rhs); }


    // General searching etc

        //- True if contained within the offsets range
        bool contains(IntType i) const noexcept
        {
            return addr_type::contains(i);
        }

        //- The max of localSizes, excluding the specified processor
        IntType maxNonLocalSize(label proci) const noexcept
        {
            return addr_type::maxNonLocalSize(proci);
        }

        //- Is on processor proci
        bool isLocal(label proci, IntType i) const
        {
            return addr_type::contains(proci, i);
        }

        //- Find processor with specified global id.
        //- Check proci first, followed by binary search.
        //  \return the processor number or -1 (not found)
        label findProc(label proci, IntType i) const
        {
            return addr_type::findRow(proci, i);
        }

        //- Find processor above proci with specified global id - binary search.
        //  \return the processor above proci or -1 otherwise
        //  (including for out-of-range indices)
        label findProcAbove(label proci, IntType i) const
        {
            return addr_type::findRowAbove(proci, i);
        }

        //- Find processor below proci with specified global id - binary search.
        //  Binary search.
        //  \return the processor below proci or -1 otherwise
        //  (including for out-of-range indices)
        label findProcBelow(label proci, IntType i) const
        {
            return addr_type::findRowBelow(proci, i);
        }

        //- Which processor does global id come from?
        //- Checks proci first (assumed to occur reasonably frequently)
        //- followed by a binary search.
        //- Fatal for out-of-range indices
        inline label whichProcID(label proci, IntType i) const;


    // Searching etc related to myProcNo
    // (Caution: uses UPstream::worldComm)

        //- The max of localSizes, excluding current (myProcNo) rank
        IntType maxNonLocalSize() const
        {
            return maxNonLocalSize(myProcNumber());
        }

        //- Is on local processor
        bool isLocal(IntType i) const
        {
            return addr_type::contains(myProcNumber(), i);
        }

        //- Find processor with specified global id.
        //- Check myProcNo first, followed by binary search.
        //  \return the processor number or -1 (not found)
        label findProc(IntType i) const
        {
            return addr_type::findRow(myProcNumber(), i);
        }

        //- Which processor does global id come from?
        //- Checks myProcNo first, followed by a binary search.
        //- Fatal for out-of-range indices
        label whichProcID(IntType i) const
        {
            return whichProcID(myProcNumber(), i);
        }


    // Offsets, sizing

        //- Start of proci data
        IntType localStart(label proci) const
        {
            return addr_type::localStart(proci);
        }

        //- Local start on myProcNo()
        IntType localStart() const { return localStart(myProcNumber()); }

        //- End of proci data
        IntType localEnd(label proci) const
        {
            return addr_type::localEnd(proci);
        }

        //- Local end on myProcNo()
        IntType localEnd() const { return localEnd(myProcNumber()); }

        //- Size of proci data
        IntType localSize(label proci) const
        {
            return addr_type::localSize(proci);
        }

        //- Local size on myProcNo()
        IntType localSize() const { return localSize(myProcNumber()); }

        //- Return start/size range of proci data
        IntRange<IntType> range(label proci) const noexcept
        {
            return addr_type::range(proci);
        }

        //- Return start/size range of local (myProcNo) data
        IntRange<IntType> range() const { return range(myProcNumber()); }

        //- Return start/size/total slab addressing for proci data
        inline OffsetRange<IntType> slice(label proci) const noexcept;

        //- Return start/size/total slab addressing for local (myProcNo) data
        OffsetRange<IntType> slice() const { return slice(myProcNumber()); }


    // Value transformations

        //- From local to global on proci
        IntType toGlobal(label proci, IntType i) const
        {
            return addr_type::toGlobal(proci, i);
        }

        //- From local to global index on myProcNo
        IntType toGlobal(IntType i) const
        {
            return addr_type::toGlobal(myProcNumber(), i);
        }

        //- From global to local on proci
        //  FatalError if not on given processor.
        IntType toLocal(label proci, IntType i) const
        {
            return addr_type::toLocal(proci, i);
        }

        //- From global to local on current processor.
        //  FatalError if not on local processor.
        IntType toLocal(IntType i) const
        {
            return addr_type::toLocal(myProcNumber(), i);
        }
};


/*---------------------------------------------------------------------------*\
                         Class globalIndex Declaration
\*---------------------------------------------------------------------------*/

class globalIndex
:
    public GlobalIndex<label>
{
    // Private Member Functions

        //- Sort and bin. validBins contains bins with non-zero size.
        static CompactListList<label>
        bin
        (
            const labelUList& offsets,
            const labelUList& globalIds,
            labelList& order,
            DynamicList<label>& validBins
        );

public:

    // STL type definitions

        //- Type of values contained
        typedef label value_type;

        //- The type that can represent the size
        typedef label size_type;


    // Static Member Functions

        //- Return a null globalIndex (reference to a nullObject).
        //- Behaves like an empty globalIndex
        static const globalIndex& null() noexcept
        {
            return NullObjectRef<globalIndex>();
        }


    // Generated Methods

        //- Default construct (empty)
        globalIndex() noexcept = default;


    // Constructors

        // Inherit all constructors
        using GlobalIndex<label>::GlobalIndex;

        //- Construct from Istream.
        //- No communication required
        explicit globalIndex(Istream& is) : GlobalIndex<label>(is) {}


    // Member Functions

        //- The span size covered by the offsets, zero if empty
        label span() const noexcept
        {
            return (end_value() - begin_value());
        }

        //- Return start/size ranges for all data
        List<labelRange> ranges() const;


    // Dimensions

        //- The number of processors covered by the offsets,
        //- same as the primary length()
        label nProcs() const noexcept { return length(); }

        //- Range of process indices for all addressed offsets (processes)
        inline labelRange allProcs() const noexcept;

        //- Range of process indices for addressed sub-offsets (processes)
        inline labelRange subProcs() const noexcept;


    // Edit

        // Reset offsets
        using GlobalIndex<label>::reset;


    // Queries and renumbering

        // Note: range() traditionally returns labelRange instead of
        // IntRange<label>, but slice is newer and is just
        // GlobalOffset<label>

        //- Return start/size range of proci data
        labelRange range(label proci) const noexcept
        {
            return GlobalIndex<label>::range(proci);
        }

        //- Return start/size range of local (myProcNo) data
        labelRange range() const
        {
            return GlobalIndex<label>::range(myProcNumber());
        }


    // Value transformations

        //- From local to global on proci, or myProcNo
        using GlobalIndex<label>::toGlobal;

        //- From local to global on proci
        inline labelList toGlobal
        (
            label proci,
            const labelUList& labels
        ) const;

        //- From local to global index on myProcNo
        labelList toGlobal(const labelUList& labels) const
        {
            return toGlobal(myProcNumber(), labels);
        }

        //- From local to global index on proci (inplace)
        inline void inplaceToGlobal
        (
            label proci,
            labelUList& labels
        ) const;

        //- From local to global index on myProcNo (inplace)
        void inplaceToGlobal(labelUList& labels) const
        {
            inplaceToGlobal(myProcNumber(), labels);
        }


    // Helper Functions

        //- Calculate ranges (offset/size) from a list of local sizes,
        //- with optional check for label overflow
        static List<labelRange> calcRanges
        (
            const labelUList& counts,
            bool checkOverflow = false
        );

        //- Split the top-level offsets into inter-node and local-node
        //- components suitable to a two-stage hierarchy.
        bool splitNodeOffsets
        (
            //! [out] Offsets between nodes (only non-empty on node leaders)
            labelList& interNodeOffsets,
            //! [out] Offsets within a node (only non-empty on node leaders)
            labelList& localNodeOffsets,
            //! The communicator. Must resolve to const world-comm
            int communicator = UPstream::worldComm,
            //! Retain absolute values for the localNode offsets
            const bool absoluteLocalNodeOffsets = false
        ) const;


    // Misc low-level gather routines

        //- Collect single values in processor order on master (== procIDs[0]).
        //  Handles contiguous/non-contiguous data.
        //  non-zero output field (master only)
        template<class ProcIDsContainer, class Type>
        [[nodiscard]]
        static List<Type> listGatherValues
        (
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            const Type& localValue,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        );

        //- Collect data in processor order on master (== procIDs[0]).
        //  Handles contiguous/non-contiguous data, skips empty fields.
        template<class ProcIDsContainer, class Type>
        static void gather
        (
            const labelUList& offsets,  //!< offsets (master only)
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            const UList<Type>& fld,     //!< [in] all ranks
            //! [out] result (master only). Must be adequately sized!
            UList<Type>& allFld,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        );

        //- Collect indirect data in processor order on master
        //  Handles contiguous/non-contiguous data, skips empty fields.
        template<class ProcIDsContainer, class Type, class Addr>
        static void gather
        (
            const labelUList& offsets,  //!< offsets (master only)
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            const IndirectListBase<Type, Addr>& fld,  //!< [in] all ranks
            //! [out] result (master only). Must be adequately sized!
            UList<Type>& allFld,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        );

        //- Inplace collect in processor order on master (== procIDs[0]).
        template<class ProcIDsContainer, class Type>
        static void gatherInplace
        (
            const labelUList& offsets,  //!< offsets (master only)
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            List<Type>& fld,            //!< [in,out]
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        );

        //- Collect data in processor order on master (== procIDs[0]).
        //  \note the globalIndex offsets needed on master only.
        template<class ProcIDsContainer, class Type>
        void gather
        (
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            const UList<Type>& fld,     //!< [in] input field
            //! [out] resized to have results on master, empty elsewhere.
            List<Type>& allFld,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        ) const;

        //- Inplace collect in processor order on master (== procIDs[0]).
        //  \note the globalIndex offsets needed on master only.
        template<class ProcIDsContainer, class Type>
        void gatherInplace
        (
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            List<Type>& fld,            //!< [in,out]
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        ) const
        {
            gatherInplace(offsets(), comm, procIDs, fld, tag, commsType);
        }


    // Gather

        //- Collect data in processor order on master
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        template<class Type>
        void gather
        (
            //! [in] input on all ranks
            const UList<Type>& sendData,
            //! [out] resized to have results on master, empty elsewhere.
            List<Type>& allData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        ) const;

        //- Collect data indirectly in processor order on master
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        template<class Type, class Addr>
        void gather
        (
            //! [in] input on all ranks
            const IndirectListBase<Type, Addr>& sendData,
            //! [out] resized to have results on master, empty elsewhere.
            List<Type>& allData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        ) const;

        //- Collect data in processor order on master
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        //
        //  \return output (master), zero-sized on non-master
        template<class Type, class OutputContainer = List<Type>>
        [[nodiscard]]
        OutputContainer gather
        (
            const UList<Type>& sendData,  //!< [in] all ranks
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        ) const;

        //- Collect data indirectly in processor order on master.
        //  Communication with default/specified communicator, message tag.
        //
        //  \return output (master), zero-sized on non-master
        template<class Type, class Addr, class OutputContainer = List<Type>>
        [[nodiscard]]
        OutputContainer gather
        (
            const IndirectListBase<Type, Addr>& sendData, //!< [in] all ranks
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        ) const;

        //- Inplace collect data in processor order on master
        //- (in serial: a no-op).
        //  Communication with default/specified communicator, message tag.
        //  After the gather, the field is zero-sized on non-master.
        template<class Type>
        void gatherInplace
        (
            List<Type>& fld,  //!< [in,out]
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        ) const;

        //- Use MPI_Gatherv call for contiguous data when possible
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator.
        //  \attention The nProcs for globalIndex and communicator
        //      must match!!
        template<class Type, class OutputContainer = List<Type>>
        void mpiGather
        (
            const UList<Type>& sendData,
            //! [out] output on master, zero-sized on non-master
            OutputContainer& allData,
            const label comm = UPstream::worldComm,  //!< communicator

            // For fallback routines:
            UPstream::commsTypes commsType = defaultCommsType(),
            const int tag = UPstream::msgType()
        ) const;

        //- Use MPI_Gatherv call for contiguous data when possible
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator.
        //  \attention The nProcs for globalIndex and communicator
        //      must match!!
        //
        //  \return output (master), zero-sized on non-master
        template<class Type, class OutputContainer = List<Type>>
        [[nodiscard]]
        OutputContainer mpiGather
        (
            const UList<Type>& sendData,
            const label comm = UPstream::worldComm,  //!< communicator

            // For fallback routines:
            UPstream::commsTypes commsType = defaultCommsType(),
            const int tag = UPstream::msgType()
        ) const;

        //- Use MPI_Gatherv call to inplace collect contiguous data
        //- when possible.
        //- (in serial: a no-op).
        //  Communication with default/specified communicator.
        //  \attention The nProcs for globalIndex and communicator
        //      must match!!
        //
        //  After the gather, the field is zero-sized on non-master.
        template<class Type>
        void mpiGatherInplace
        (
            //! [in,out]
            List<Type>& fld,
            const label comm = UPstream::worldComm,  //!< communicator

            // For fallback routines:
            UPstream::commsTypes commsType = defaultCommsType(),
            const int tag = UPstream::msgType()
        ) const;


    // Gather Operations

        //- Use MPI_Gatherv call to collect contiguous data when possible
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator.
        //
        //  The allData is output (master), zero-sized on non-master
        template<class Type, class OutputContainer = List<Type>>
        static void mpiGatherOp
        (
            const UList<Type>& sendData,
            OutputContainer& allData,
            const label comm = UPstream::worldComm,  //!< communicator

            // For fallback routines:
            UPstream::commsTypes commsType = defaultCommsType(),
            const int tag = UPstream::msgType()
        );

        //- Use MPI_Gatherv call to collect contiguous data when possible
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator.
        //
        //  \return output (master), zero-sized on non-master
        template<class Type, class OutputContainer = List<Type>>
        [[nodiscard]]
        static OutputContainer mpiGatherOp
        (
            const UList<Type>& sendData,
            const label comm = UPstream::worldComm,  //!< communicator

            // For fallback routines:
            UPstream::commsTypes commsType = defaultCommsType(),
            const int tag = UPstream::msgType()
        );

        //- Use MPI_Gatherv call to inplace collect contiguous data
        //- when possible.
        //- (in serial: a no-op).
        //  Communication with default/specified communicator.
        //
        //  After the gather, the field is zero-sized on non-master.
        template<class Type>
        static void mpiGatherInplaceOp
        (
            //! [in,out]
            List<Type>& fld,
            const label comm = UPstream::worldComm,  //!< communicator

            // For fallback routines:
            UPstream::commsTypes commsType = defaultCommsType(),
            const int tag = UPstream::msgType()
        );

        //- Collect data in processor order on master
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        template<class Type>
        static void gatherOp
        (
            const UList<Type>& sendData,
            //! [out] output on master, zero-sized on non-master
            List<Type>& allData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Collect data in processor order on master
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        template<class Type, class Addr>
        static void gatherOp
        (
            const IndirectListBase<Type, Addr>& sendData,
            //! [out] output on master, zero-sized on non-master
            List<Type>& allData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Collect and return data in processor order on master
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        //
        //  \return output (master), zero-sized on non-master
        template<class Type, class OutputContainer = List<Type>>
        [[nodiscard]]
        static OutputContainer gatherOp
        (
            const UList<Type>& sendData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Collect and return data in processor order on master
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        //
        //  \return output (master), zero-sized on non-master
        template<class Type, class Addr, class OutputContainer = List<Type>>
        [[nodiscard]]
        static OutputContainer gatherOp
        (
            const IndirectListBase<Type, Addr>& sendData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Inplace collect data in processor order on master
        //- (in serial: a no-op).
        //  Communication with default/specified communicator, message tag.
        //
        //  After the gather, the field is zero-sized on non-master.
        template<class Type>
        static void gatherInplaceOp
        (
            //! [in,out]
            List<Type>& fld,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        );


    // Scatter

        //- Distribute data in processor order.
        //  Requires fld to be correctly sized!
        //  Communication with default/specified communicator, message tag.
        template<class ProcIDsContainer, class Type>
        static void scatter
        (
            const labelUList& offsets,  //!< offsets (master only)
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            const UList<Type>& allFld,
            UList<Type>& fld,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        );

        //- Distribute data in processor order.
        //  Requires fld to be correctly sized!
        //  Communication with default/specified communicator, message tag.
        //  \note the globalIndex offsets needed on master only.
        template<class ProcIDsContainer, class Type>
        void scatter
        (
            const label comm,           //!< communicator
            const ProcIDsContainer& procIDs,
            const UList<Type>& allFld,
            UList<Type>& fld,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType()
        ) const
        {
            scatter(offsets(), comm, procIDs, allFld, fld, tag, commsType);
        }

        //- Distribute data in processor order.
        //  Requires fld to be correctly sized!
        //  Communication with default/specified communicator, message tag.
        //  \note the globalIndex offsets needed on master only.
        template<class Type>
        void scatter
        (
            const UList<Type>& allData,
            UList<Type>& localData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        ) const;

        //- Distribute data in processor order
        //- (in serial: performs a simple copy).
        //  Communication with default/specified communicator, message tag.
        //  \note the globalIndex offsets needed on master only.
        template<class Type, class OutputContainer = List<Type>>
        [[nodiscard]]
        OutputContainer scatter
        (
            const UList<Type>& allData,
            const int tag = UPstream::msgType(),
            UPstream::commsTypes commsType = defaultCommsType(),
            const label comm = UPstream::worldComm  //!< communicator
        ) const;


    // Scatter

        //- Get (potentially remote) data.
        //- Elements required given as global indices
        //  Communication with default/specified communicator, message tag.
        template<class Type, class CombineOp>
        void get
        (
            List<Type>& allFld,
            const labelUList& globalIds,
            const CombineOp& cop,
            const label comm = UPstream::worldComm,  //!< communicator
            const int tag = UPstream::msgType()
        ) const;


    // Housekeeping

        //- Dispatch tag:
        //- Construct 'one-sided' from the non-master local sizes
        //- using gather but no broadcast.
        struct gatherNonLocal{};

        //- Construct 'one-sided' from the non-master local sizes
        //- (ie, master size is treated as 0).
        //- The master will have all content, sub-ranks no content.
        //
        //  \note can be used when UPstream::parRun() is false.
        globalIndex
        (
            globalIndex::gatherNonLocal,
            label localSize,
            int communicator = UPstream::worldComm  //!< communicator
        )
        {
            // non-master sizes only
            if (UPstream::master(communicator)) localSize = 0;
            reset(globalIndex::gatherOnly{}, localSize, communicator);
        }

        //- Reset as 'one-sided' from the non-master local sizes
        //- (ie, master size is treated as 0).
        //- The master will have all content, sub-ranks no content.
        //
        //  \note can be used when UPstream::parRun() is false.
        void reset
        (
            globalIndex::gatherNonLocal,
            label localSize,
            int communicator = UPstream::worldComm  //!< communicator
        )
        {
            // non-master sizes only
            if (UPstream::master(communicator)) localSize = 0;
            reset(globalIndex::gatherOnly{}, localSize, communicator);
        }

        //- Construct from local size, using gather/broadcast
        //- with default/specified communicator if parallel.
        FOAM_DEPRECATED_FOR(2022-03, "construct without message tag")
        globalIndex
        (
            label localSize,
            const int tag,          // message tag (unused)
            const label comm,       // communicator
            const bool parallel     // use parallel comms
        )
        {
            reset(localSize, comm, parallel);
        }

        //- Reset from local size, using gather/broadcast
        //- with default/specified communicator if parallel.
        FOAM_DEPRECATED_FOR(2022-03, "reset without message tag")
        void reset
        (
            label localSize,
            const int tag,          // message tag (unused)
            const label comm,       // communicator
            const bool parallel     // use parallel comms
        )
        {
            reset(localSize, comm, parallel);
        }

        //- Prefer localStart() to avoid confusing with offsets()
        FOAM_DEPRECATED_FOR(2022-02, "use localStart()")
        label offset(const label proci) const { return localStart(proci); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "globalIndexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "globalIndex.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
