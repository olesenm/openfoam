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
#include "CompactListList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class globalIndex;

Istream& operator>>(Istream& is, globalIndex& gi);
Ostream& operator<<(Ostream& os, const globalIndex& gi);

/*---------------------------------------------------------------------------*\
                         Class globalIndex Declaration
\*---------------------------------------------------------------------------*/

class globalIndex
{
    // Private Data

        //- Offset (addressing) table - like CompactListList
        labelList offsets_;


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

        //- Report overflow at specified (non-negative) index
        static void reportOverflowAndExit
        (
            const label idx,
            const label prevOffset = -1,  // The last valid offset value
            const label count = 0         // The count to add to prevOffset
        );

        //- Check for overflow and calls reportOverflowAndExit...
        static inline void checkOverflowAndExit
        (
            const label idx,
            const label prevOffset,     // The last valid offset value
            const label count           // The count to add to prevOffset
        );

        //- Return start/size range addressing for proci data
        inline labelRange range_impl(label proci) const noexcept;

        //- Return start/size/total slab addressing for proci data
        inline OffsetRange<label> slice_impl(label proci) const noexcept;


public:

    // STL type definitions

        //- Type of values contained
        typedef label value_type;

        //- The type that can represent the size
        typedef label size_type;


    // Public Data Types

        //- Dispatch tag:
        //- Construct with a single (local size) entry, no communication.
        struct gatherNone{};

        //- Dispatch tag:
        //- Construct 'one-sided' from local sizes,
        //- using gather but no broadcast.
        struct gatherOnly{};

        //- Dispatch tag:
        //- Construct 'one-sided' from the non-master local sizes
        //- using gather but no broadcast.
        struct gatherNonLocal{};


    // Static Member Functions

        //- Return a null globalIndex (reference to a nullObject).
        //- Behaves like an empty globalIndex
        static const globalIndex& null() noexcept
        {
            return NullObjectRef<globalIndex>();
        }


    // Properties

        //- Default communication type for gather/scatter
        static constexpr UPstream::commsTypes defaultCommsType() noexcept
        {
            return UPstream::commsTypes::nonBlocking;
        }


    // Constructors

        //- Default construct (empty)
        globalIndex() noexcept = default;

        //- Copy construct from a list of offsets.
        //- No communication required
        inline explicit globalIndex(const labelUList& listOffsets);

        //- Move construct from a list of offsets.
        //- No communication required
        inline explicit globalIndex(labelList&& listOffsets);

        //- Construct from a list of sizes and calculate the offsets.
        //- No communication required
        inline globalIndex
        (
            globalIndex::gatherNone,
            const labelUList& localSizes
        );

        //- Construct from local size, using gather/broadcast
        //- with default/specified communicator if parallel.
        inline explicit globalIndex
        (
            label localSize,
            const label comm = UPstream::worldComm,   //!< communicator
            const bool parallel = UPstream::parRun()  //!< use parallel comms
        );

        //- Construct with a single (local size) entry, no communication
        inline globalIndex
        (
            globalIndex::gatherNone,
            label localSize,
            int communicator = -1   //!< (unused) no communicator needed
        );

        //- Construct 'one-sided' from local sizes.
        //- Uses UPstream::listGatherValues, but no broadcast.
        //- Will be empty on non-master processes.
        //
        //  \note can be used when Pstream::parRun() is false.
        inline globalIndex
        (
            globalIndex::gatherOnly,
            label localSize,
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Construct 'one-sided' from the non-master local sizes
        //- (ie, master size is treated as 0).
        //- Uses UPstream::listGatherValues, but no broadcast.
        //- Will be empty on non-master processes.
        //
        //  \note can be used when Pstream::parRun() is false.
        inline globalIndex
        (
            globalIndex::gatherNonLocal,
            label localSize,
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Construct from Istream.
        //- No communication required
        explicit globalIndex(Istream& is);


    // Member Functions

        //- Check for default constructed or total-size == 0
        inline bool empty() const noexcept;

        //- The number of items covered by the offsets
        inline label length() const noexcept;

        //- Global sum of localSizes. Same as totalSize()
        FOAM_DEPRECATED_STRICT(2022-10, "totalSize() - unambiguous")
        label size() const noexcept { return totalSize(); }

        //- The span size covered by the offsets, zero if empty
        inline label span() const noexcept;

        //- The total addressed size, which corresponds to the
        //- end offset and also the sum of all localSizes.
        inline label totalSize() const noexcept;

        //- The local sizes. Same as localSizes()
        inline labelList sizes() const;

        //- The local starts
        inline const labelUList localStarts() const;

        //- The local sizes
        labelList localSizes() const;

        //- Return start/size ranges for all data
        List<labelRange> ranges() const;

        //- Global max of localSizes
        inline label maxSize() const;


    // Access

        //- Const-access to the offsets
        inline const labelList& offsets() const noexcept;

        //- Write-access to the offsets, for changing after construction
        inline labelList& offsets() noexcept;


    // Dimensions

        //- True if local-only content (ie, nProcs == 1).
        //- Such content is often created with gatherNone.
        inline bool single() const noexcept;

        //- The number of processors covered by the offsets,
        //- same as the primary length()
        inline label nProcs() const noexcept;

        //- Range of process indices for all addressed offsets (processes)
        inline labelRange allProcs() const noexcept;

        //- Range of process indices for addressed sub-offsets (processes)
        inline labelRange subProcs() const noexcept;

        //- The value corresponding to the first offset
        inline label begin_value() const noexcept;

        //- The value corresponding to the last offset (end offset),
        //- which is 1 beyond the end of the range.
        inline label end_value() const noexcept;

        //- The first offset range. It is (0,0) if globalIndex is empty
        labelRange front() const;

        //- The last offset range. It is (0,0) if globalIndex is empty
        labelRange back() const;


    // Edit

        //- Reset to be empty (no offsets)
        inline void clear();

        //- Change the number of entries (nProcs) in the offsets table.
        //- Extending will fill with empty local sizes.
        void resize(const label n);

        //- Reset from local size, using gather/broadcast
        //- with default/specified communicator if parallel.
        void reset
        (
            label localSize,
            const label comm = UPstream::worldComm,   //!< communicator
            const bool parallel = UPstream::parRun()  //!< use parallel comms
        );

        //- Reset offsets from a list of local sizes,
        //- with optional check for label overflow.
        //- No communication required
        void reset
        (
            const labelUList& counts,
            const bool checkOverflow = false
        );

        //- Reset to a single (local size) entry, no communication
        inline void reset
        (
            globalIndex::gatherNone,
            label localSize,
            int communicator = -1   //!< (unused) no communicator needed
        );

        //- Reset as 'one-sided' from local sizes.
        //- Uses UPstream::listGatherValues, but no broadcast.
        //- Will be empty on non-master processes.
        //
        //  \note can be used when Pstream::parRun() is false.
        inline void reset
        (
            globalIndex::gatherOnly,
            label localSize,
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Reset as 'one-sided' from the non-master local sizes
        //- (ie, master size is treated as 0).
        //- Uses UPstream::listGatherValues, but no broadcast.
        //- Will be empty on non-master processes.
        //
        //  \note can be used when Pstream::parRun() is false.
        inline void reset
        (
            globalIndex::gatherNonLocal,
            label localSize,
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Reset the globalIndex. Same as copy assignment.
        inline void reset(const globalIndex& gi);

        //- Alter local size for given processor
        void setLocalSize(const label proci, const label len);


    // Queries and renumbering

        //- True if contained within the offsets range
        inline bool contains(const label i) const noexcept;

        //- Start of proci data
        inline label localStart(const label proci) const;

        //- End of proci data
        inline label localEnd(const label proci) const;

        //- Size of proci data
        inline label localSize(const label proci) const;

        //- The max of localSizes, excluding the specified processor
        label maxNonLocalSize(const label proci) const;

        //- Return start/size range of proci data
        inline labelRange range(label proci) const noexcept;

        //- Return start/size/total slab addressing for proci data
        inline OffsetRange<label> slice(label proci) const noexcept;

        //- Is on processor proci
        inline bool isLocal(const label proci, const label i) const;

        //- From local to global on proci
        inline label toGlobal(const label proci, const label i) const;

        //- From local to global on proci
        inline labelList toGlobal
        (
            const label proci,
            const labelUList& labels
        ) const;

        //- From local to global index on proci (inplace)
        inline void inplaceToGlobal
        (
            const label proci,
            labelUList& labels
        ) const;


        //- From global to local on proci
        inline label toLocal(const label proci, const label i) const;

        //- Find processor with specified global id.
        //- Check proci first, followed by binary search.
        //  \return the processor number or -1 (not found)
        inline label findProc(const label proci, const label i) const;

        //- Find processor above proci with specified global id - binary search.
        //  \return the processor above proci or -1 otherwise
        //  (including for out-of-range indices)
        inline label findProcAbove(const label proci, const label i) const;

        //- Find processor below proci with specified global id - binary search.
        //  Binary search.
        //  \return the processor below proci or -1 otherwise
        //  (including for out-of-range indices)
        inline label findProcBelow(const label proci, const label i) const;

        //- Which processor does global id come from?
        //- Checks proci first (assumed to occur reasonably frequently)
        //- followed by a binary search.
        //- Fatal for out-of-range indices
        inline label whichProcID(const label proci, const label i) const;


    // Queries relating to myProcNo (Caution: uses UPstream::worldComm)

        //- Local start on myProcNo()
        inline label localStart() const;

        //- Local end on myProcNo()
        inline label localEnd() const;

        //- Local size on myProcNo()
        inline label localSize() const;

        //- The max of localSizes, excluding current (myProcNo) rank
        inline label maxNonLocalSize() const;

        //- Return start/size range of local (myProcNo) data
        inline labelRange range() const;

        //- Return start/size/total slab addressing for local (myProcNo) data
        inline OffsetRange<label> slice() const;

        //- Is on local processor
        inline bool isLocal(const label i) const;

        //- From local to global index
        inline label toGlobal(const label i) const;

        //- From local to global index
        inline labelList toGlobal(const labelUList& labels) const;

        //- From local to global index (inplace)
        inline void inplaceToGlobal(labelUList& labels) const;

        //- From global to local on current processor.
        //  FatalError if not on local processor.
        inline label toLocal(const label i) const;

        //- Which processor does global id come from?
        //  Uses myProcNo for the initial local check.
        inline label whichProcID(const label i) const;


    // Iteration

        //- Forward input iterator with const access that is used to
        //- iterate across the globalIndex offsets() table.
        //  The de-referenced value is the range() with (start, size),
        //  but it also provides separate index, start, size information.
        class const_iterator
        {
            //- The parent class for which this is an iterator
            const globalIndex* parent_;

            //- The index into the parent
            label index_;

        public:

        // Constructors

            //- Construct from globalIndex list at given index
            explicit const_iterator
            (
                const globalIndex* globalIdx,
                const label i = 0
            ) noexcept;


        // Member Operators

            //- The index into the arrays
            inline label index() const noexcept;

            //- The local start
            inline label start() const;

            //- The local size
            inline label size() const;

            //- The local range
            inline labelRange range() const;

            //- The local range
            inline labelRange operator*() const;

            inline const_iterator& operator++();
            inline const_iterator operator++(int);

            inline const_iterator& operator--();
            inline const_iterator operator--(int);

            inline bool operator==(const const_iterator& iter) const noexcept;
            inline bool operator!=(const const_iterator& iter) const noexcept;
        };


        //- A const_iterator set to the beginning
        inline const_iterator cbegin() const noexcept;

        //- A const_iterator set to beyond the end
        inline const const_iterator cend() const noexcept;

        //- A const_iterator set to the beginning
        inline const_iterator begin() const noexcept;

        //- A const_iterator set to beyond the end
        inline const const_iterator end() const noexcept;

        //- Return const_iterator at offset proci from begin,
        //- clamped to [0,nProcs] range
        inline const_iterator cbegin(const label proci) const noexcept;

        //- Return const_iterator at offset proci from begin,
        //- clamped to [0,nProcs] range
        inline const_iterator begin(const label proci) const noexcept;

        // FUTURE?
        // //- Return start/size range of proci data
        // labelRange operator[](const label proci) const
        // {
        //     return this->range(proci);
        // }


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
            [[maybe_unused]] const bool checkOverflow
        )
        {
            IntType start = UPstream::find_offset(localSize, communicator);

            // Assume no overflow for large integrals
            if constexpr (sizeof(IntType) < sizeof(int64_t))
            {
                if (checkOverflow && UPstream::is_parallel(communicator))
                {
                    checkOverflowAndExit
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
            const label comm,  //!< communicator
            const bool checkOverflow
        )
        {
            IntType start =
                globalIndex::calcOffset(localSize, comm, checkOverflow);
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
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Calculate offsets from a list of local sizes,
        //- with optional check for label overflow
        static labelList calcOffsets
        (
            const labelUList& counts,
            const bool checkOverflow = false
        );

        //- Calculate offsets from an indirect list of local sizes,
        //- with optional check for label overflow
        template<class Addr>
        static labelList calcOffsets
        (
            const IndirectListBase<label, Addr>& counts,
            const bool checkOverflow = false
        );

        //- Calculate offsets from list of lists,
        //- with optional check for label overflow
        template<class SubListType>
        static labelList calcListOffsets
        (
            const List<SubListType>& lists,
            const bool checkOverflow = false
        );

        //- Calculate ranges (offset/size) from a list of local sizes,
        //- with optional check for label overflow
        static List<labelRange> calcRanges
        (
            const labelUList& counts,
            const bool checkOverflow = false
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
            const label communicator = UPstream::worldComm,
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
            gatherInplace(offsets_, comm, procIDs, fld, tag, commsType);
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
            scatter(offsets_, comm, procIDs, allFld, fld, tag, commsType);
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


    // Member Operators

        //- Compare for equality - uses the offsets
        bool operator==(const globalIndex& rhs) const
        {
            return (this->offsets() == rhs.offsets());
        }

        //- Compare for inequality - uses the offsets
        bool operator!=(const globalIndex& rhs) const
        {
            return !(*this == rhs);
        }

        //- Compare for less-than - uses the offsets
        bool operator<(const globalIndex& rhs) const
        {
            return (this->offsets() < rhs.offsets());
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, globalIndex& gi);
        friend Ostream& operator<<(Ostream& os, const globalIndex& gi);


    // Housekeeping

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
