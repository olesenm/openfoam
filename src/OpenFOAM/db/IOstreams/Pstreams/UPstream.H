/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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
    Foam::UPstream

Description
    Inter-processor communications stream

SourceFiles
    UPstream.C
    UPstream.txx
    UPstreamCommsStruct.C
    UPstreamReduceOffsets.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstream_H
#define Foam_UPstream_H

#include "wordList.H"
#include "labelList.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "Map.H"
#include "Enum.H"
#include "ListOps.H"
#include "OffsetRange.H"
#include "OffsetRangeIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Implementation details for UPstream/Pstream/MPI etc.
namespace PstreamDetail {}

//- Interface handling for UPstream/Pstream/MPI etc.
namespace PstreamUtils {}

/*---------------------------------------------------------------------------*\
                          Class UPstream Declaration
\*---------------------------------------------------------------------------*/

class UPstream
{
public:

    //- Int ranges are used for MPI ranks (processes)
    typedef IntRange<int> rangeType;

    //- Communications types
    enum class commsTypes : char
    {
        buffered,       //!< "buffered"                 : (MPI_Bsend, MPI_Recv)
        scheduled,      //!< "scheduled" (MPI standard) : (MPI_Send, MPI_Recv)
        nonBlocking,    //!< "nonBlocking" (immediate)  : (MPI_Isend, MPI_Irecv)
        // Aliases
        blocking = buffered  //!< compatibility name for buffered
    };

    //- Enumerated names for the communication types
    static const Enum<commsTypes> commsTypeNames;

    //- Different MPI-send modes (ignored for commsTypes::buffered)
    enum class sendModes : char
    {
        normal,         //!< (MPI_Send, MPI_Isend)
        sync            //!< (MPI_Ssend, MPI_Issend)
    };

    //- Mapping of some fundamental and aggregate types to MPI data types
    enum class dataTypes : char
    {
        // NOTE: changes here require adjustment in
        // PstreamGlobals, UPstreamTraits

        // Fundamental Types [10]:
        Basic_begin,    //!< (internal use) begin marker [basic/all types]
        type_byte = Basic_begin,  //!< byte, char, unsigned char, ...
        type_int16,
        type_int32,
        type_int64,
        type_uint16,
        type_uint32,
        type_uint64,
        type_float,
        type_double,
        type_long_double,

        //! (internal use) end marker [basic types]
        Basic_end, User_begin = Basic_end,
        //!< (internal use) begin marker [user types]

        // User Types [6]:
        type_3float = User_begin, //!< 3*float (eg, floatVector)
        type_3double,   //!< 3*double (eg, doubleVector)
        type_6float,    //!< 6*float (eg, floatSymmTensor, complex vector)
        type_6double,   //!< 6*double (eg, doubleSymmTensor, complex vector)
        type_9float,    //!< 9*float (eg, floatTensor)
        type_9double,   //!< 9*double (eg, doubleTensor)

        // Internal markers [1]
        invalid,        //!< invalid type (NULL)
        //! (internal use) end marker [user types]
        User_end = invalid, DataTypes_end = invalid
        //!< (internal use) end marker [all types]
    };

    //- Mapping of some MPI op codes.
    //  Currently excluding min/max location until they are needed
    enum class opCodes : char
    {
        // NOTE: changes here require adjustment in
        // PstreamGlobals, UPstreamTraits

        Basic_begin,  //!< (internal use) begin marker [reduce/window types]

        // Reduce or window operations [10]
        op_min = Basic_begin, //!< min(x,y)
        op_max,         //!< max(x,y)
        op_sum,         //!< (x + y)
        op_prod,        //!< (x * y)
        op_bool_and,    //!< Logical \c and
        op_bool_or,     //!< Logical \c or
        op_bool_xor,    //!< Logical \c xor
        op_bit_and,     //!< Bit-wise \c and for (unsigned) integral types
        op_bit_or,      //!< Bit-wise \c or for (unsigned) integral types
        op_bit_xor,     //!< Bit-wise \c xor for (unsigned) integral types

        //! (internal use) end marker [reduce types]
        Basic_end, Extra_begin = Basic_end,
        //!< (internal use) begin marker [window types]

        // Window-only operations [2]
        op_replace = Extra_begin,  //!< Replace (window only)
        op_no_op,       //!< No-op (window only)

        // Internal markers [1]
        invalid,        //!< invalid op (NULL)
        //! (internal use) end marker [window types]
        Extra_end = invalid, OpCodes_end = invalid
        //!< (internal use) end marker [all types]
    };

    //- Some bit masks corresponding to topology controls
    //  These selectively enable topology-aware handling
    enum class topoControls : int
    {
        broadcast = 1,      //!< broadcast [MPI]
        reduce = 2,         //!< reduce/all-reduce [MPI]
        gather = 4,         //!< gather/all-gather [MPI]
        combine = 16,       //!< combine/gather (reduction) [manual algorithm]
        mapGather = 32,     //!< mapGather (reduction) [manual algorithm]
        gatherList = 64,    //!< gatherList [manual algorithm]
    };


    // Public Classes

        //- Wrapper for OpenFOAM internal communicator index
        class communicator;  // Forward Declaration

        //- Wrapper for MPI_Comm
        class Communicator;  // Forward Declaration

        //- Wrapper for MPI_File routines
        class File;  // Forward Declaration

        //- Wrapper for MPI_Request
        class Request;  // Forward Declaration

        //- Wrapper for MPI_Win
        class Window;  // Forward Declaration

        //- Structure for communicating between processors
        class commsStruct
        {
            // Private Data

                //- The rank id \em directly above (-1 for master)
                int above_;

                //- The rank ids \em directly below (ie, leaders)
                List<int> below_;

                //- The ids of \em all ranks below (not just directly below).
                //  Corresponds to a depth-first order of nodes.
                List<int> allBelow_;

                //- The procIDs of all processors not below myProcNo
                //- (inverse of allBelow_ without myProcNo)
                List<int> allNotBelow_;


        public:

            // Constructors

                //- Default construct with above == -1
                commsStruct() noexcept : above_(-1) {}

                //- Move construct from components
                commsStruct
                (
                    const int above,
                    List<int>&& below,
                    List<int>&& allBelow,
                    List<int>&& allNotBelow
                );

                //- Copy construct from below, allBelow components
                commsStruct
                (
                    const int numProcs,
                    const int myProcID,
                    const int above,
                    const UList<int>& below,
                    const UList<int>& allBelow
                );


            // Member Functions

                //- Root of the tree (ie, nothing above)
                bool is_root() const noexcept { return (above_ < 0); }

                //- Internal node in the tree
                bool is_inode() const noexcept
                {
                    return (above_ >= 0 && !below_.empty());
                }

                //- Leaf/terminal in the tree (ie, nothing below)
                bool is_leaf() const noexcept { return below_.empty(); }


            // Access

                //- The rank id \em directly above (-1 for master)
                int above() const noexcept { return above_; }

                //- The rank ids \em directly below (ie, leaders)
                const List<int>& below() const noexcept { return below_; }

                //- The ids of \em all ranks below (not just directly below).
                //  Corresponds to a depth-first order of nodes.
                const List<int>& allBelow() const noexcept { return allBelow_; }

                //- The procIDs of all processors not below myProcNo.
                //- The inverse set of allBelow without myProcNo.
                //  No particular guarantee on the ordering.
                const List<int>& allNotBelow() const noexcept
                {
                    return allNotBelow_;
                }

                //- The number of processors addressed by the structure
                int nProcs() const noexcept;


            // Edit

                //- Reset to default constructed state
                void reset();

                //- Reset to linear (flat) communication
                void reset_linear
                (
                    const int myProci,
                    const int numProcs
                );

                //- Reset to tree selection (linear when very small),
                //- possibly with communicator-specific adjustments
                void reset
                (
                    const int myProci,
                    const int numProcs,
                    const int communicator
                );


            // Member / Friend Operators

                bool operator==(const commsStruct&) const;
                bool operator!=(const commsStruct&) const;

                friend Ostream& operator<<(Ostream&, const commsStruct&);
        };

        //- Collection of communication structures
        class commsStructList
        {
            // Private Data

                //- The communicator index
                int comm_;

                //- Linear (flat) communication instead of tree communication
                bool flat_;

                //- The communication tree
                List<commsStruct> tree_;

        public:

            // Constructors

                //- Construct empty with invalid communicator
                commsStructList() noexcept
                :
                    comm_(-1),
                    flat_(false)
                {}

                //- Construct empty with given communicator
                explicit commsStructList(int comm, bool flat=false) noexcept
                :
                    comm_(comm),
                    flat_(flat)
                {}


            // Static Functions

                //- An empty structure. Used for placeholders etc.
                static const commsStructList& null();


            // Member Functions

                //- True if communicator is non-negative (ie, was assigned)
                bool good() const noexcept { return (comm_ >= 0); }

                //- The communicator internal index
                int comm() const noexcept { return comm_; }

                //- Linear (flat) communication instead of tree communication
                bool linear() const noexcept { return flat_; }

                //- Change communication type (linear vs tree communication).
                //- Resets existing structure.
                void linear(bool on);

                //- Clear the list
                void clear() { return tree_.clear(); }

                //- True if the list is empty
                bool empty() const noexcept { return tree_.empty(); }

                //- The number of entries
                label size() const noexcept { return tree_.size(); }

                //- Reset communicator index, fill tree with empty entries
                void init(int communicator);

                //- Reset communicator index, fill tree with empty entries.
                //- Specify flat vs tree communication.
                void init(int communicator, bool flat);

                //- Reset communicator index, clear tree entries
                void reset(int communicator);

                //- Reset communicator index, clear tree entries.
                //- Specify flat vs tree communication.
                void reset(int communicator, bool flat);

                //- Get existing or create (demand-driven) entry
                const UPstream::commsStruct& get(int proci) const;

                //- Get existing or create (demand-driven) entry
                const UPstream::commsStruct& operator[](int proci) const
                {
                    return get(proci);
                }

                //- Print un-directed graph in graphviz dot format
                void printGraph(Ostream& os, int proci = 0) const;

                //- Write the communication tree
                Ostream& writeList(Ostream& os) const
                {
                    return tree_.writeList(os);
                }


            // Member / Friend Operators

                //- Write the communication tree
                friend Ostream& operator<<
                (
                    Ostream& os,
                    const commsStructList& rhs
                )
                {
                    return rhs.writeList(os);
                }
        };


private:

    // Private Data

        //- Communications type of this stream
        commsTypes commsType_;


    // Private Static Data

        //- By default this is not a parallel run
        static bool parRun_;

        //- Have support for threads?
        static bool haveThreads_;

        //- Initial MPI_Comm_dup(MPI_COMM_WORLD,...) disabled? (default: false)
        static bool noInitialCommDup_;

        //- Standard transfer message type
        static int msgType_;

        //- Index to the world-communicator as defined at startup
        //- (after any multi-world definitions).
        //- Is unaffected by any later changes to worldComm.
        static int constWorldComm_;

        //- Index to the inter-node communicator (between nodes),
        //- defined based on constWorldComm_
        static int commInterNode_;

        //- Index to the intra-host communicator (within a node),
        //- defined based on constWorldComm_
        static int commLocalNode_;

        //- The number of shared/host nodes in the (const) world communicator.
        static int numNodes_;

        //- Names of all worlds
        static wordList allWorlds_;

        //- Per processor the world index (into allWorlds_)
        static labelList worldIDs_;


    // Communicator specific data

        //- My processor number
        static DynamicList<int> myProcNo_;

        //- List of process IDs
        static DynamicList<List<int>> procIDs_;

        //- Parent communicator
        static DynamicList<label> parentComm_;

        //- Free communicators
        static DynamicList<label> freeComms_;

        //- Linear communication schedule
        static DynamicList<commsStructList> linearCommunication_;

        //- Multi level communication schedule
        static DynamicList<commsStructList> treeCommunication_;


    // Private Member Functions

        //- Set data for parallel running
        static void setParRun(const int nProcs, const bool haveThreads);

        //- Initialise entries for new communicator.
        //
        //  Resets corresponding entry in myProcNo_, procIDs_,
        //  linearCommunication_, treeCommunication_
        //  \return the communicator index
        static label getAvailableCommIndex(const label parentIndex);

        //- Define inter-host/intra-host communicators (uses commConstWorld).
        //  Optionally specify a given number per node.
        static bool setHostCommunicators(const int numPerNode = 0);

        //- Define inter-host/intra-host communicators based on
        //- shared-memory information. Uses comm-world.
        static bool setSharedMemoryCommunicators();

        //- Allocate MPI components of communicator with given index.
        //  This represents a "top-down" approach, creating a communicator
        //  based on the procIDs_ groupings.
        //
        //  Modifies myProcNo_, reads and modifies procIDs_
        static void allocateCommunicatorComponents
        (
            const label parentIndex,
            const label index
        );

        //- Allocate MPI components as duplicate of the parent communicator
        //
        //  Modifies myProcNo_, procIDs_
        static void dupCommunicatorComponents
        (
            const label parentIndex,
            const label index
        );

        //- Allocate MPI components for the given index by splitting
        //- the parent communicator on the given \em colour.
        //  This represents a "bottom-up" approach, when the individual ranks
        //  only know which group they should belong to, but don't yet know
        //  which other ranks will be in their group.
        //
        //  Modifies myProcNo_, procIDs_
        static void splitCommunicatorComponents
        (
            const label parentIndex,
            const label index,
            const int colour,
            //! Use MPI_Allgather+MPI_Comm_create_group vs MPI_Comm_split
            const bool two_step = true
        );

        //- Free MPI components of communicator.
        //  Does not touch the first two communicators (SELF, WORLD)
        static void freeCommunicatorComponents(const label index);


    // Private implementation helpers

        //- Test for communicator equality.
        //  True if they have the same index or address the same ranks
        static bool sameProcs_impl(int comm1, int comm2)
        {
            return
            (
                (comm1 == comm2) ||
                (
                    // With guard against bad index
                    (comm1 >= 0 && comm1 <= procIDs_.size())
                 && (comm2 >= 0 && comm2 <= procIDs_.size())
                 && (procIDs_[comm1] == procIDs_[comm2])
                )
            );
        }

        //- Test equality of communicator procs with the given list of ranks
        template<typename IntType>
        static bool sameProcs_impl
        (
            int communicator,
            const UList<IntType>& procs
        )
        {
            return
            (
                // With guard against bad index
                (communicator >= 0 && communicator <= procIDs_.size())
             && ListOps::equal(procIDs_[communicator], procs)
            );
        }

        //- Test the equality of two lists of ranks
        template<typename Type1, typename Type2>
        static bool sameProcs_impl
        (
            const UList<Type1>& procs1,
            const UList<Type2>& procs2
        )
        {
            return ListOps::equal(procs1, procs2);
        }


protected:

    // Protected Member Functions

    // Static Functions

        //- Broadcast buffer contents to all ranks (default: from rank=0).
        //- The sizes must match on all processes.
        //  For \b non-parallel : do nothing.
        //
        //  \note The method uses a \c void pointer and the required data type
        //  (as per MPI). This means it should almost never be called directly
        //  but always via a compile-time checked caller.
        //  \return True on success
        static bool mpi_broadcast
        (
            void* buf,
            std::streamsize count,
            const UPstream::dataTypes dataTypeId,
            const int communicator,
            //! The broadcast root (0 == master, -1 = last rank)
            int root = 0
        );

        //- In-place reduction of \c values with result on rank 0.
        //  Includes internal parallel guard and checks on data types, opcode.
        //
        //  \note The method uses a \c void pointer and the required data type
        //  (as per MPI). This means it should almost never be called directly
        //  but always via a compile-time checked caller.
        static void mpi_reduce
        (
            void* values,
            int count,
            const UPstream::dataTypes dataTypeId,
            const UPstream::opCodes opCodeId,
            const int communicator,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        );

        //- In-place reduction of \c values with same result on all ranks.
        //  Includes internal parallel guard and checks on data types, opcode.
        //
        //  \note The method uses a \c void pointer and the required data type
        //  (as per MPI). This means it should almost never be called directly
        //  but always via a compile-time checked caller.
        static void mpi_allreduce
        (
            void* values,
            int count,
            const UPstream::dataTypes dataTypeId,
            const UPstream::opCodes opCodeId,
            const int communicator,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        );

        //- In-place scan/exscan reduction of \c values.
        //  Includes internal parallel guard and checks on data types, opcode.
        //
        //  \note The method uses a \c void pointer and the required data type
        //  (as per MPI). This means it should almost never be called directly
        //  but always via a compile-time checked caller.
        static void mpi_scan_reduce
        (
            void* values,
            int count,
            const UPstream::dataTypes dataTypeId,
            const UPstream::opCodes opCodeId,
            const int communicator,
            //! Use exclusive scan
            const bool exclusive
        );

        //- Send buffer contents of specified data type to given processor.
        //
        //  \note The method uses a \c void pointer and the required data type
        //  (as per MPI). This means it should almost never be called directly
        //  but always via a compile-time checked caller.
        //  \return True on success (or Fatal)
        static bool mpi_send
        (
            const UPstream::commsTypes commsType,
            const void* buf,
            std::streamsize count,
            const UPstream::dataTypes dataTypeId,
            const int toProcNo,
            const int tag,              // eg, UPstream::msgType()
            const int communicator,     // eg, UPstream::worldComm
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr,
            //! optional send mode (normal | sync)
            const UPstream::sendModes sendMode = UPstream::sendModes::normal
        );

        //- Receive buffer contents of specified data type from given processor.
        //
        //  \note The method uses a \c void pointer and the required data type
        //  (as per MPI). This means it should almost never be called directly
        //  but always via a compile-time checked caller.
        //  The commsType will be ignored if UPstream::Request is specified.
        //
        //  \return number of elements read. May change in the future
        static std::streamsize mpi_receive
        (
            const UPstream::commsTypes commsType,
            void* buf,
            std::streamsize count,
            const UPstream::dataTypes dataTypeId,
            const int fromProcNo,
            const int tag,              // eg, UPstream::msgType()
            const int communicator,     // eg, UPstream::worldComm
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        );

        //- Receive identically-sized (contiguous) data from all ranks,
        //- placing the result on rank 0.
        //  Includes internal parallel guard.
        //  For non-parallel, does not copy any data.
        //  If needed, this must be done by the caller.
        static void mpi_gather
        (
            //! All ranks: location of individual value to send
            const void* sendData,
            //! Master: receive buffer with all values.
            //! Other ranks: ignored
            void* recvData,
            //! Number of send/recv data per rank. Globally consistent!
            int count,
            const UPstream::dataTypes dataTypeId,
            const int communicator,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        );

        //- Send identically-sized (contiguous) data from rank 0
        //- to all other ranks.
        //  Includes internal parallel guard.
        static void mpi_scatter
        (
            //! Master: send buffer with all values.
            //! Other ranks: ignored
            const void* sendData,
            //! All ranks: location to receive individual value
            void* recvData,
            //! Number of send/recv data per rank. Globally consistent!
            int count,
            const UPstream::dataTypes dataTypeId,
            const int communicator,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        );

        //- Gather/scatter identically-sized data
        //  Send data from proc slot, receive into all slots
        static void mpi_allgather
        (
            //! All ranks: the base of the data locations
            void* allData,
            //! Number of send/recv data per rank. Globally consistent!
            int count,
            const UPstream::dataTypes dataTypeId,
            const int communicator,
            //! [out] request information (for non-blocking)
            UPstream::Request* req = nullptr
        );

        //- Receive variable length data from all ranks,
        //- placing the result on rank 0.
        //- (caution: known to scale poorly)
        static void mpi_gatherv
        (
            const void* sendData,
            int sendCount,    //!< Ignored on master if recvCount[0] == 0
            void* recvData,                 //!< Ignored on non-root rank
            const UList<int>& recvCounts,   //!< Ignored on non-root rank
            const UList<int>& recvOffsets,  //!< Ignored on non-root rank
            const UPstream::dataTypes dataTypeId,
            const int communicator
        );

        //- Send variable length data from rank 0 to all ranks.
        //- (caution: known to scale poorly)
        static void mpi_scatterv
        (
            const void* sendData,           //!< Ignored on non-root rank
            const UList<int>& sendCounts,   //!< Ignored on non-root rank
            const UList<int>& sendOffsets,  //!< Ignored on non-root rank
            void* recvData,
            int recvCount,
            const UPstream::dataTypes dataTypeId,
            const int communicator
        );

public:

    //- Declare name of the class and its debug switch
    ClassName("UPstream");


    // Static Data

        //- Use of host/node topology-aware routines
        //  0: disabled
        //  1: split by hostname [default]
        //  2: split by shared
        //  >=4: (debug) split with given number per node
        static int nodeCommsControl_;

        //- Minimum number of nodes before topology-aware routines are enabled
        //  <= 2 : always
        //  >= 3 : when there are more than N nodes
        static int nodeCommsMin_;

        //- Selection of topology-aware routines as a bitmask combination
        //- of the topoControls enumerations
        static int topologyControl_;

        //- Test for selection of given topology-aware routine
        static bool usingTopoControl(UPstream::topoControls ctrl) noexcept
        {
            return static_cast<bool>(topologyControl_ & int(ctrl));
        }

        //- Should compact transfer be used in which floats replace doubles
        //- reducing the bandwidth requirement at the expense of some loss
        //- in accuracy
        static bool floatTransfer;

        //- Number of processors to change from linear to tree communication
        static int nProcsSimpleSum;

        //- Number of processors to change to nonBlocking consensual
        //- exchange (NBX). Ignored for zero or negative values.
        static int nProcsNonblockingExchange;

        //- Number of polling cycles in processor updates
        static int nPollProcInterfaces;

        //- Default commsType
        static commsTypes defaultCommsType;

        //- Optional maximum message size (bytes)
        static int maxCommsSize;

        //- Tuning parameters for non-blocking exchange (NBX)
        static int tuning_NBX_;

        //- MPI buffer-size (bytes)
        static const int mpiBufferSize;


    // Standard Communicators

        //- Communicator for all ranks.
        //- May differ from commGlobal() if local worlds are in use
        static label worldComm;

        //- Debugging: warn for use of any communicator differing from warnComm
        static label warnComm;

        //- Communicator for all ranks, irrespective of any local worlds.
        //  This value \em never changes during a simulation.
        static constexpr int commGlobal() noexcept { return 0; }

        //- Communicator within the current rank only
        //  This value \em never changes during a simulation.
        static constexpr int commSelf() noexcept { return 1; }

        //- Communicator for all ranks (respecting any local worlds).
        //  This value \em never changes after startup. Unlike the commWorld()
        //  which can be temporarily overriden.
        static int commConstWorld() noexcept { return constWorldComm_; }

        //- Communicator for all ranks (respecting any local worlds)
        static label commWorld() noexcept { return worldComm; }

        //- Set world communicator. Negative values are a no-op.
        //  \returns old world communicator index
        static label commWorld(label communicator) noexcept
        {
            label old(worldComm);
            if (communicator >= 0) worldComm = communicator;
            return old;
        }

        //- The current communicator debugging setting
        static label commWarn() noexcept { return warnComm; }

        //- Alter communicator debugging setting.
        //- Warns for use of any communicator differing from specified.
        //- Negative values disable.
        //  \returns the previous warn index
        static label commWarn(label communicator) noexcept
        {
            label old(warnComm);
            warnComm = communicator;
            return old;
        }

        //- Number of currently defined communicators
        static label nComms() noexcept { return parentComm_.size(); }

        //- Debugging: print the communication tree
        static void printCommTree(int communicator, bool linear = false);


    // Host Communicators

        //- Communicator between nodes/hosts (respects any local worlds)
        static int commInterNode() noexcept
        {
            return (parRun_ ? commInterNode_ : constWorldComm_);
        }

        //- Communicator within the node/host (respects any local worlds)
        static int commLocalNode() noexcept
        {
            return (parRun_ ? commLocalNode_ : constWorldComm_);
        }

        //- Both inter-node and local-node communicators have been created
        static bool hasNodeCommunicators() noexcept
        {
            return
            (
                (commInterNode_ > constWorldComm_)
             && (commLocalNode_ > constWorldComm_)
            );
        }

        //- True if node topology-aware routines have been enabled,
        //- it is running in parallel, the starting point is the
        //- world-communicator and it is not an odd corner case
        //- (ie, all processes on one node, all processes on different nodes)
        static bool usingNodeComms(const int communicator);


    // Constructors

        //- Construct for given communication type
        explicit UPstream(const commsTypes commsType) noexcept
        :
            commsType_(commsType)
        {}


    // Member Functions

        //- Create new communicator with sub-ranks on the parent communicator
        [[nodiscard]]
        static label newCommunicator
        (
            //! The parent communicator
            const label parent,

            //! The contiguous sub-ranks of parent to use
            const labelRange& subRanks,

            //! Call allocateCommunicatorComponents()
            const bool withComponents = true
        );

        //- Create new communicator with sub-ranks on the parent communicator
        [[nodiscard]]
        static label newCommunicator
        (
            //! The parent communicator
            const label parent,

            //! The sub-ranks of parent to use (ignore negative values)
            const labelUList& subRanks,

            //! Call allocateCommunicatorComponents()
            const bool withComponents = true
        );

        //- Duplicate the parent communicator
        //
        //  Always calls dupCommunicatorComponents() internally
        [[nodiscard]]
        static label dupCommunicator
        (
            //! The parent communicator
            const label parent
        );

        //- Allocate a new communicator by splitting the parent communicator
        //- on the given \em colour.
        //  Always calls splitCommunicatorComponents() internally
        [[nodiscard]]
        static label splitCommunicator
        (
            //! The parent communicator
            const label parent,

            //! The colouring to select which ranks to include.
            //! Negative values correspond to 'ignore'
            const int colour,

            //! Use MPI_Allgather+MPI_Comm_create_group vs MPI_Comm_split
            const bool two_step = true
        );

        //- Free a previously allocated communicator.
        //  Ignores placeholder (negative) communicators.
        static void freeCommunicator
        (
            const label communicator,
            const bool withComponents = true
        );

        //- Return physical processor number (i.e. processor number in
        //- worldComm) given communicator and processor
        static int baseProcNo(label comm, int procID);

        //- Return processor number in communicator (given physical processor
        //- number) (= reverse of baseProcNo)
        static label procNo(const label comm, const int baseProcID);

        //- Return processor number in communicator (given processor number
        //- and communicator)
        static label procNo
        (
            const label comm,
            const label currentComm,
            const int currentProcID
        );

        //- Add the valid option this type of communications library
        //- adds/requires on the command line
        static void addValidParOptions(HashTable<string>& validParOptions);

        //- Initialisation function called from main
        //  Spawns sub-processes and initialises inter-communication
        static bool init(int& argc, char**& argv, const bool needsThread);

        //- Special purpose initialisation function.
        //  Performs a basic MPI_Init without any other setup.
        //  Only used for applications that need MPI communication when
        //  OpenFOAM is running in a non-parallel mode.
        //  \note Behaves as a no-op if MPI has already been initialized.
        //      Fatal if MPI has already been finalized.
        static bool initNull();

        //- Impose a synchronisation barrier (optionally non-blocking)
        static void barrier
        (
            const int communicator,
            UPstream::Request* req = nullptr
        );

        //- Impose a point-to-point synchronisation barrier
        //- by sending a zero-byte \em "done" message to given rank.
        //  A no-op for non-parallel
        static void send_done
        (
            //! The destination rank
            const int toProcNo,
            //! The communicator index (eg, UPstream::worldComm)
            const int communicator,
            //! Message tag (must match on receiving side)
            const int tag = UPstream::msgType()+1970
        );

        //- Impose a point-to-point synchronisation barrier
        //- by receiving a zero-byte \em "done" message from given rank
        //  A no-op for non-parallel
        //  \returns the source rank (useful for ANY_SOURCE messages)
        //  or -1 for non-parallel
        static int wait_done
        (
            //! The source rank (negative == ANY_SOURCE)
            const int fromProcNo,
            //! The communicator index (eg, UPstream::worldComm)
            const int communicator,
            //! Message tag (must match on sending side)
            const int tag = UPstream::msgType()+1970
        );

        //- A standard or non-blocking probe for an incoming message.
        //
        //  \returns source rank and message size (bytes)
        //      and (-1, 0) on failure
        [[nodiscard]]
        static std::pair<int,int64_t> probeMessage
        (
            //! Non-blocking or otherwise
            const UPstream::commsTypes commsType,
            //! The source rank (negative == ANY_SOURCE)
            const int fromProcNo,
            //! The source message tag
            const int tag = UPstream::msgType(),
            //! The communicator index
            const int communicator = UPstream::worldComm
        );

        //- A standard probe for an incoming message.
        //
        //  \returns source rank and message size (bytes)
        //      and (-1, 0) on failure
        [[nodiscard]]
        static std::pair<int,int64_t> probeMessage
        (
            //! The source rank (negative == ANY_SOURCE)
            const int fromProcNo,
            //! The source message tag
            const int tag = UPstream::msgType(),
            //! The communicator index
            const int communicator = UPstream::worldComm
        )
        {
            return probeMessage
            (
                UPstream::commsTypes::scheduled,  // A standard probe
                fromProcNo,
                tag,
                communicator
            );
        }

        //- Probe incoming message sizes from specified ranks.
        static void probeMessages
        (
            //! [out] Message sizes (bytes)
            DynamicList<int64_t>& messageSizes,
            //! The source ranks (must be non-negative)
            const UList<int>& fromProcs,
            //! The source message tag: assumed to be identical for all sources
            const int tag = UPstream::msgType(),
            //! The communicator index
            const int communicator = UPstream::worldComm
        );


    // Information

        //- Report the node-communication settings
        static void printNodeCommsControl(Ostream& os);

        //- Report the topology routines settings
        static void printTopoControl(Ostream& os);


    // Requests (non-blocking comms).
    // Pending requests are usually handled as an internal (global) list,
    // since this simplifies the overall tracking and provides a convenient
    // wrapping to avoid exposing MPI-specific types, but can also handle
    // with a wrapped UPstream::Request as well.

        //- Number of outstanding requests (on the internal list of requests)
        static label nRequests() noexcept;

        //- Truncate outstanding requests to given length, which is
        //- expected to be in the range [0 to nRequests()].
        //  A no-op for out-of-range values.
        static void resetRequests(const label n);

        //- Transfer the (wrapped) MPI request to the internal global list
        //- and invalidate the parameter (ignores null requests)
        //  A no-op for non-parallel
        static void addRequest(UPstream::Request& req);

        //- Non-blocking comms: cancel and free outstanding request.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op if parRun() == false
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        static void cancelRequest(const label i);

        //- Non-blocking comms: cancel and free outstanding request.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op for non-parallel
        static void cancelRequest(UPstream::Request& req);

        //- Non-blocking comms: cancel and free outstanding requests.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op if parRun() == false or list is empty
        static void cancelRequests(UList<UPstream::Request>& requests);

        //- Non-blocking comms: cancel/free outstanding requests
        //- (from position onwards) and remove from internal list of requests.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op if parRun() == false,
        //  if the position is out-of-range [0 to nRequests()],
        //  or the internal list of requests is empty.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to remove (negative = until the end)
        static void removeRequests(label pos, label len = -1);

        //- Non-blocking comms: free outstanding request.
        //- Corresponds to MPI_Request_free()
        //  A no-op if parRun() == false
        static void freeRequest(UPstream::Request& req);

        //- Non-blocking comms: free outstanding requests.
        //- Corresponds to MPI_Request_free()
        //  A no-op if parRun() == false or list is empty
        static void freeRequests(UList<UPstream::Request>& requests);

        //- Wait until all requests (from position onwards) have finished.
        //- Corresponds to MPI_Waitall()
        //  A no-op if parRun() == false,
        //  if the position is out-of-range [0 to nRequests()],
        //  or the internal list of requests is empty.
        //
        //  If checking a trailing portion of the list, it will also trim
        //  the list of outstanding requests as a side-effect.
        //  This is a feature (not a bug) to conveniently manange the list.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to check (negative = until the end)
        static void waitRequests(label pos, label len = -1);

        //- Wait until all requests have finished.
        //- Corresponds to MPI_Waitall()
        //  A no-op if parRun() == false, or the list is empty.
        static void waitRequests(UList<UPstream::Request>& requests);

        //- Wait until any request (from position onwards) has finished.
        //- Corresponds to MPI_Waitany()
        //  A no-op and returns false if parRun() == false,
        //  if the position is out-of-range [0 to nRequests()],
        //  or the internal list of requests is empty.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to check (negative = until the end)
        //
        //  \returns true if any pending request completed.
        //  \returns false if all requests have already been handled.
        static bool waitAnyRequest(label pos, label len = -1);

        //- Wait until some requests (from position onwards) have finished.
        //- Corresponds to MPI_Waitsome()
        //  A no-op and returns false if parRun() == false,
        //  if the position is out-of-range [0 to nRequests],
        //  or the internal list of requests is empty.
        //
        //  \param[out] indices the completed request indices relative to the
        //      starting position. This is an optional parameter that can be
        //      used to recover the indices or simply to avoid reallocations
        //      when calling within a loop.
        //
        //  \returns true if some pending requests completed.
        //  \returns false if all requests have already been handled
        static bool waitSomeRequests
        (
            //! starting position within the internal list of requests
            label pos,
            //! length of slice to check (negative = until the end)
            label len = -1,
            DynamicList<int>* indices = nullptr
        );

        //- Wait until some requests (from position onwards) have finished.
        //- Corresponds to MPI_Waitsome()
        //  A no-op and returns false if parRun() == false,
        //  if the position is out-of-range [0 to nRequests],
        //  or the internal list of requests is empty.
        //
        //  \returns true if some pending requests completed.
        //  \returns false if all requests have already been handled
        static bool waitSomeRequests
        (
            //! starting position within the internal list of requests
            label pos,
            //! [out] completed request indices relative to \p pos
            DynamicList<int>& indices
        )
        {
            return waitSomeRequests(pos, -1, &indices);
        }

        //- Wait until some requests have finished.
        //- Corresponds to MPI_Waitsome()
        //  A no-op and returns false if parRun() == false,
        //  the list is empty,
        //  or if all the requests have already been handled.
        //
        //  \param[out] indices the completed request indices relative to the
        //      starting position. This is an optional parameter that can be
        //      used to recover the indices or simply to avoid reallocations
        //      when calling within a loop.
        static bool waitSomeRequests
        (
            //! The requests
            UList<UPstream::Request>& requests,
            DynamicList<int>* indices = nullptr
        );

        //- Wait until any request has finished and return its index.
        //- Corresponds to MPI_Waitany()
        //  Returns -1 if parRun() == false, or the list is empty,
        //  or if all the requests have already been handled
        static int waitAnyRequest(UList<UPstream::Request>& requests);

        //- Wait until request i has finished.
        //- Corresponds to MPI_Wait()
        //  A no-op if parRun() == false,
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        static void waitRequest(const label i);

        //- Wait until specified request has finished.
        //- Corresponds to MPI_Wait()
        //  A no-op if parRun() == false or for a null-request
        static void waitRequest(UPstream::Request& req);

        //- Is request \p i active (!= MPI_REQUEST_NULL)?
        //  False if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        static bool activeRequest(const label i);

        //- Is request active (!= MPI_REQUEST_NULL)?
        static bool activeRequest(const UPstream::Request& req);

        //- Non-blocking comms: has request i finished?
        //- Corresponds to MPI_Test()
        //  A no-op and returns true if parRun() == false,
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        static bool finishedRequest(const label i);

        //- Non-blocking comms: has request finished?
        //- Corresponds to MPI_Test()
        //  A no-op and returns true if parRun() == false
        //  or for a null-request
        static bool finishedRequest(UPstream::Request& req);

        //- Non-blocking comms: have all requests (from position onwards)
        //- finished?
        //- Corresponds to MPI_Testall()
        //  A no-op and returns true if parRun() == false,
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        //  or the addressed range is empty etc.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to check (negative = until the end)
        static bool finishedRequests(label pos, label len = -1);

        //- Non-blocking comms: have all requests finished?
        //- Corresponds to MPI_Testall()
        //  A no-op and returns true if parRun() == false or list is empty
        static bool finishedRequests(UList<UPstream::Request>& requests);

        //- Non-blocking comms: have both requests finished?
        //- Corresponds to pair of MPI_Test()
        //  A no-op and returns true if parRun() == false,
        //  if there are no pending requests,
        //  or if the indices are out-of-range (0 to nRequests)
        //  Each finished request parameter is set to -1 (ie, done).
        static bool finishedRequestPair(label& req0, label& req1);

        //- Non-blocking comms: wait for both requests to finish.
        //- Corresponds to pair of MPI_Wait()
        //  A no-op if parRun() == false,
        //  if there are no pending requests,
        //  or if the indices are out-of-range (0 to nRequests)
        //  Each finished request parameter is set to -1 (ie, done).
        static void waitRequestPair(label& req0, label& req1);


    // General

        //- Set as parallel run on/off.
        //  \return the previous value
        static bool parRun(bool on) noexcept
        {
            bool old(parRun_);
            parRun_ = on;
            return old;
        }

        //- Test if this a parallel run
        //  Modify access is deprecated
        static bool& parRun() noexcept { return parRun_; }

        //- Have support for threads
        static bool haveThreads() noexcept { return haveThreads_; }

        //- Relative rank for the master process - is always 0.
        static constexpr int masterNo() noexcept { return 0; }

        //- Number of ranks in parallel run (for given communicator).
        //- It is 1 for serial run
        static label nProcs(label communicator = worldComm)
        {
            return procIDs_[communicator].size();
        }

        //- Rank of this process in the communicator (starting from masterNo()).
        //- Negative if the process is not a rank in the communicator.
        static int myProcNo(label communicator = worldComm)
        {
            return myProcNo_[communicator];
        }

        //- True if process corresponds to the master rank in the communicator
        static bool master(label communicator = worldComm)
        {
            return myProcNo_[communicator] == masterNo();
        }

        //- True if process corresponds to the master rank in the communicator
        static bool is_master(label communicator = worldComm)
        {
            return UPstream::master(communicator);
        }

        //- True if process corresponds to \b any rank (master or sub-rank)
        //- in the given communicator
        static bool is_rank(label communicator = worldComm)
        {
            return myProcNo_[communicator] >= masterNo();
        }

        //- True if process corresponds to a sub-rank in the given communicator
        static bool is_subrank(label communicator = worldComm)
        {
            return myProcNo_[communicator] > masterNo();
        }

        //- True if parallel algorithm or exchange is required.
        //  This is when parRun() == true, the process corresponds to a rank
        //  in the communicator and there is more than one rank in the
        //  communicator
        static bool is_parallel(label communicator = worldComm)
        {
            return
            (
                parRun_ && is_rank(communicator) && nProcs(communicator) > 1
            );
        }

        //- The number of shared/host nodes in the (const) world communicator.
        static int numNodes() noexcept { return numNodes_; }

        //- The parent communicator
        static label parent(int communicator)
        {
            return parentComm_(communicator);
        }

        //- The list of ranks within a given communicator
        static List<int>& procID(int communicator)
        {
            return procIDs_[communicator];
        }


        //- Test for communicator equality.
        //  True if they have the same index or address the same ranks
        static bool sameProcs(int communicator1, int communicator2)
        {
            return sameProcs_impl(communicator1, communicator2);
        }

        //- Test equality of communicator procs with the given list of ranks.
        //- Includes a guard for the communicator index.
        template
        <
            typename T1,
            typename = std::void_t
            <std::enable_if_t<std::is_integral_v<T1>>>
        >
        static bool sameProcs(int communicator, const UList<T1>& procs)
        {
            return sameProcs_impl(communicator, procs);
        }

        //- Test the equality of two lists of ranks
        template
        <
            typename T1,
            typename T2,
            typename = std::void_t
            <
                std::enable_if_t<std::is_integral_v<T1>>,
                std::enable_if_t<std::is_integral_v<T2>>
            >
        >
        static bool sameProcs(const UList<T1>& procs1, const UList<T2>& procs2)
        {
            return sameProcs_impl(procs1, procs2);
        }


    // Worlds

        //- All worlds
        static const wordList& allWorlds() noexcept
        {
            return allWorlds_;
        }

        //- The indices into allWorlds for all processes
        static const labelList& worldIDs() noexcept
        {
            return worldIDs_;
        }

        //- My worldID
        static label myWorldID()
        {
            return worldIDs_[myProcNo(UPstream::commGlobal())];
        }

        //- My world
        static const word& myWorld()
        {
            return allWorlds_[worldIDs_[myProcNo(UPstream::commGlobal())]];
        }


    // Rank addressing

        //- Range of process indices for all processes
        static rangeType allProcs(label communicator = worldComm)
        {
            // Proc 0 -> nProcs (int value)
            return rangeType(static_cast<int>(nProcs(communicator)));
        }

        //- Range of process indices for sub-processes
        static rangeType subProcs(label communicator = worldComm)
        {
            // Proc 1 -> nProcs (int value)
            return rangeType(1, static_cast<int>(nProcs(communicator)-1));
        }

        //- Processor offsets corresponding to the inter-node communicator
        static const List<int>& interNode_offsets();

        //- The range (start/size) of the commLocalNode ranks in terms of the
        //- (const) world communicator processors.
        static const rangeType& localNode_parentProcs();

        //- Linear communication schedule (special case) for given communicator
        static const commsStructList& linearCommunication(int communicator);

        //- Tree communication schedule (standard case) for given communicator
        static const commsStructList& treeCommunication(int communicator);

        //- Communication schedule for all-to-master (proc 0) as
        //- linear/tree/none with switching based on UPstream::nProcsSimpleSum,
        //- the is_parallel() state and the optional \c linear parameter.
        static const commsStructList& whichCommunication
        (
            const int communicator,
            //! optionally select linear schedule only
            bool linear = false
        )
        {
            const label np
            (
                parRun_ && is_rank(communicator)  // cf. is_parallel()
              ? nProcs(communicator)
              : label(0)
            );

            return
            (
                np <= 1
              ? commsStructList::null()
              : (linear || np < UPstream::nProcsSimpleSum)
              ? linearCommunication(communicator)
              : treeCommunication(communicator)
            );
        }


        //- Message tag of standard messages
        static int& msgType() noexcept
        {
            return msgType_;
        }

        //- Set the message tag for standard messages
        //  \return the previous value
        static int msgType(int val) noexcept
        {
            int old(msgType_);
            msgType_ = val;
            return old;
        }

        //- Increment the message tag for standard messages
        //  \return the previous value
        static int incrMsgType(int val = 1) noexcept
        {
            int old(msgType_);
            msgType_ += val;
            return old;
        }

        //- Get the communications type of the stream
        commsTypes commsType() const noexcept
        {
            return commsType_;
        }

        //- Set the communications type of the stream
        //  \return the previous value
        commsTypes commsType(const commsTypes ct) noexcept
        {
            commsTypes old(commsType_);
            commsType_ = ct;
            return old;
        }


        //- Shutdown (finalize) MPI as required.
        //  Uses MPI_Abort instead of MPI_Finalize if errNo is non-zero
        static void shutdown(int errNo = 0);

        //- Call MPI_Abort with no other checks or cleanup
        static void abort(int errNo = 1);

        //- Shutdown (finalize) MPI as required and exit program with errNo.
        static void exit(int errNo = 1);


        #undef  Pstream_declareMethods
        #define Pstream_declareMethods(Type)                                  \
                                                                              \
        /*!\brief Exchange \c Type data with all ranks in communicator     */ \
        /*! \em non-parallel : simple copy of \p sendData to \p recvData   */ \
        static void allToAll                                                  \
        (                                                                     \
            /*! [in] The value at [proci] is sent to proci */                 \
            const UList<Type>& sendData,                                      \
            /*! [out] The data received from the other ranks */               \
            UList<Type>& recvData,                                            \
            const int communicator = UPstream::worldComm                      \
        );                                                                    \
                                                                              \
        /*!\brief Exchange \em non-zero \c Type data between ranks [NBX]   */ \
        /*! \p recvData is always initially assigned zero and no non-zero  */ \
        /*! values are sent/received from other ranks.                     */ \
        /*! \em non-parallel : simple copy of \p sendData to \p recvData   */ \
        /*! \note The message tag should be chosen to be a unique value    */ \
        /*! since the implementation uses probing with ANY_SOURCE !!       */ \
        /*! An initial barrier may help to avoid synchronisation problems  */ \
        /*! caused elsewhere (See "nbx.tuning" opt switch)                 */ \
        static void allToAllConsensus                                         \
        (                                                                     \
            /*! [in] The \em non-zero value at [proci] is sent to proci */    \
            const UList<Type>& sendData,                                      \
            /*! [out] The non-zero value received from each rank */           \
            UList<Type>& recvData,                                            \
            /*! Message tag for the communication */                          \
            const int tag,                                                    \
            const int communicator = UPstream::worldComm                      \
        );                                                                    \
                                                                              \
        /*!\brief Exchange \c Type data between ranks [NBX]                */ \
        /*! \p recvData map is always cleared initially so a simple check  */ \
        /*! of its keys is sufficient to determine connectivity.           */ \
        /*! \em non-parallel : copy own rank (if it exists)                */ \
        /*! See notes about message tags and "nbx.tuning" opt switch       */ \
        static void allToAllConsensus                                         \
        (                                                                     \
            /*! [in] The value at [proci] is sent to proci. */                \
            const Map<Type>& sendData,                                        \
            /*! [out] The values received from given ranks. */                \
            Map<Type>& recvData,                                              \
            /*! Message tag for the communication */                          \
            const int tag,                                                    \
            const int communicator = UPstream::worldComm                      \
        );                                                                    \
                                                                              \
        /*!\brief Exchange \c Type data between ranks [NBX]                */ \
        /*! \returns any received data as a Map                            */ \
        [[nodiscard]]                                                         \
        static Map<Type> allToAllConsensus                                    \
        (                                                                     \
            /*! [in] The value at [proci] is sent to proci. */                \
            const Map<Type>& sendData,                                        \
            /*! Message tag for the communication */                          \
            const int tag,                                                    \
            const int communicator = UPstream::worldComm                      \
        )                                                                     \
        {                                                                     \
            Map<Type> recvData;                                               \
            allToAllConsensus(sendData, recvData, tag, communicator);         \
            return recvData;                                                  \
        }

        Pstream_declareMethods(int32_t);
        Pstream_declareMethods(int64_t);

        #undef Pstream_declareMethods


    // Low-level gather/scatter routines

        //- Receive identically-sized (contiguous) data from all ranks
        template<class Type>
        static void mpiGather
        (
            //! All ranks: location of individual value to send
            const Type* sendData,
            //! Master: receive buffer with all values.
            //! Other ranks: ignored
            Type* recvData,
            //! Number of send/recv data per rank. Globally consistent!
            int count,
            const int communicator = UPstream::worldComm
        );

        //- Send identically-sized (contiguous) data to all ranks
        template<class Type>
        static void mpiScatter
        (
            //! Master: send buffer with all values.
            //! Other ranks: ignored
            const Type* sendData,
            //! All ranks: location to receive individual value
            Type* recvData,
            //! Number of send/recv data per rank. Globally consistent!
            int count,
            const int communicator = UPstream::worldComm
        );

        //- Gather/scatter identically-sized data
        //  Send data from proc slot, receive into all slots
        template<class Type>
        static void mpiAllGather
        (
            //! All ranks: the base of the data locations
            Type* allData,
            //! Number of send/recv data per rank. Globally consistent!
            int count,
            const int communicator = UPstream::worldComm
        );

        //- Receive variable length data from all ranks
        template<class Type>
        static void mpiGatherv
        (
            const Type* sendData,
            int sendCount,  //!< Ignored on master if recvCount[0] == 0
            Type* recvData,                 //!< Ignored on non-root rank
            const UList<int>& recvCounts,   //!< Ignored on non-root rank
            const UList<int>& recvOffsets,  //!< Ignored on non-root rank
            const int communicator = UPstream::worldComm
        );

        //- Send variable length data to all ranks
        template<class Type>
        static void mpiScatterv
        (
            const Type* sendData,           //!< Ignored on non-root rank
            const UList<int>& sendCounts,   //!< Ignored on non-root rank
            const UList<int>& sendOffsets,  //!< Ignored on non-root rank
            Type* recvData,
            int recvCount,
            const int communicator = UPstream::worldComm
        );


    // Gather single, contiguous value(s)

        //- Allgather individual values into list locations.
        //  The returned list has size nProcs, identical on all ranks.
        template<class T>
        [[nodiscard]]
        static List<T> allGatherValues
        (
            const T& localValue,
            const int communicator = UPstream::worldComm
        );

        //- Gather individual values into list locations.
        //  On master list length == nProcs, otherwise zero length.
        //  \n
        //  For \b non-parallel :
        //  the returned list length is 1 with localValue.
        template<class T>
        [[nodiscard]]
        static List<T> listGatherValues
        (
            const T& localValue,
            const int communicator = UPstream::worldComm
        );

        //- Scatter individual values from list locations.
        //  On master input list length == nProcs, ignored on other procs.
        //  \n
        //  For \b non-parallel :
        //  returns the first list element (or default initialized).
        template<class T>
        [[nodiscard]]
        static T listScatterValues
        (
            const UList<T>& allValues,
            const int communicator = UPstream::worldComm
        );


    // Broadcast Functions

        //- Broadcast buffer contents (contiguous types) to all ranks
        //- (default: from rank=0). The sizes must match on all processes.
        //  For \b non-parallel : do nothing.
        //  \return True on success
        template<class Type>
        inline static bool broadcast
        (
            Type* buffer,
            std::streamsize count,
            const int communicator,
            //! The broadcast root (0 == master, -1 = last rank)
            const int root = UPstream::masterNo()
        );

        //- Broadcast fixed-list content (contiguous types) to all ranks
        //- (default: from rank=0). The sizes must match on all processes.
        //  For \b non-parallel : do nothing.
        //  \return True on success
        template<class Type, unsigned N>
        inline static bool broadcast
        (
            FixedList<Type, N>& list,
            const int communicator,
            //! The broadcast root (0 == master, -1 = last rank)
            const int root = UPstream::masterNo()
        );


    // Reductions

        //- MPI_Reduce (blocking) for known operators
        //  For \b non-parallel : do nothing.
        template<class T>
        static void mpiReduce
        (
            T values[],
            int count,
            const UPstream::opCodes opCodeId,
            const int communicator = UPstream::worldComm
        );

        //- MPI_Reduce (blocking) for known operators
        //  For \b non-parallel : do nothing.
        template<UPstream::opCodes opCode, class T>
        static void mpiReduce
        (
            T values[],
            int count,
            const int communicator = UPstream::worldComm
        );

        //- MPI_Allreduce (blocking) for known operators
        //  For \b non-parallel : do nothing.
        template<class T>
        static void mpiAllReduce
        (
            T values[],
            int count,
            const UPstream::opCodes opCodeId,
            const int communicator = UPstream::worldComm
        );

        //- MPI_Allreduce (blocking) for known operators
        //  For \b non-parallel : do nothing.
        template<UPstream::opCodes opCode, class T>
        static void mpiAllReduce
        (
            T values[],
            int count,
            const int communicator = UPstream::worldComm
        );

        //- MPI_Iallreduce (non-blocking) for known operators
        //  For \b non-parallel : do nothing.
        template<class T>
        static void mpiAllReduce
        (
            T values[],
            int count,
            const UPstream::opCodes opCodeId,
            const int communicator,
            UPstream::Request& req
        );

        //- MPI_Iallreduce (non-blocking) for known operators
        //  For \b non-parallel : do nothing.
        template<UPstream::opCodes opCode, class T>
        static void mpiAllReduce
        (
            T values[],
            int count,
            const int communicator,
            UPstream::Request& req
        );

        #undef  Pstream_declareMethods
        #define Pstream_declareMethods(OpName, OpCode)                        \
                                                                              \
        /*!\brief MPI Reduce \b \c OpName */                                  \
        template<class T>                                                     \
        static void mpiReduce_##OpName                                        \
        (                                                                     \
            T values[],                                                       \
            int count,                                                        \
            const int communicator = UPstream::worldComm                      \
        )                                                                     \
        {                                                                     \
            UPstream::mpiReduce(values, count, OpCode, communicator);         \
        }                                                                     \
                                                                              \
        /*!\brief MPI Allreduce \b \c OpName */                               \
        template<class T>                                                     \
        static void mpiAllReduce_##OpName                                     \
        (                                                                     \
            T values[],                                                       \
            int count,                                                        \
            const int communicator = UPstream::worldComm                      \
        )                                                                     \
        {                                                                     \
            UPstream::mpiAllReduce(values, count, OpCode, communicator);      \
        }

        Pstream_declareMethods(min, UPstream::opCodes::op_min);
        Pstream_declareMethods(max, UPstream::opCodes::op_max);
        Pstream_declareMethods(sum, UPstream::opCodes::op_sum);
        #undef Pstream_declareMethods


    // Partial reductions

        //- Inclusive/exclusive scan (in-place).
        //  In exclusive mode, the degenerate value on rank=0 has no meaning
        //  and will normally be treated like non-exclusive mode
        //  (ie, original values).
        //  However, for the \b \c opCodes::op_sum type we can provide a
        //  sensible value and overwrite it with a zero-initialized value.
        //
        //  \note For \b non-parallel : do nothing.
        template<Foam::UPstream::opCodes opCode, class T>
        static void mpiScan
        (
            T values[],
            int count,
            const int communicator,
            //! Use exclusive scan (MPI_Exscan vs MPI_Scan)
            const bool exclusive = false
        );

        //- Inclusive/exclusive scan returning the result.
        //- In exclusive mode, the degenerate value on rank=0 has no meaning
        //- but is treated like non-exclusive mode (ie, original values)
        //
        //  \note For \b non-parallel : return own value
        template<Foam::UPstream::opCodes opCode, class T>
        [[nodiscard]]
        static T mpiScan
        (
            const T& localValue,
            const int communicator,
            //! Use exclusive scan (MPI_Exscan vs MPI_Scan)
            const bool exclusive = false
        );

        #undef  Pstream_declareMethods
        #define Pstream_declareMethods(OpName, OpCode)                        \
                                                                              \
        /*!\brief Inclusive/exclusive \b \c OpName scan (in-place) */         \
        /*! See UPstream::mpiScan for notes about the behaviour */            \
        template<class T>                                                     \
        static void mpiScan_##OpName                                          \
        (                                                                     \
            T values[],                                                       \
            int count,                                                        \
            const int communicator,                                           \
            /*! Use exclusive scan (MPI_Exscan vs MPI_Scan) */                \
            const bool exclusive = false                                      \
        )                                                                     \
        {                                                                     \
            UPstream::mpiScan<OpCode>(values, count, communicator, exclusive);\
        }                                                                     \
                                                                              \
        /*!\brief Exclusive \b \c OpName scan (in-place) */                   \
        /*! See UPstream::mpiScan for notes about the behaviour */            \
        template<class T>                                                     \
        static void mpiExscan_##OpName                                        \
        (                                                                     \
            T values[],                                                       \
            int count,                                                        \
            const int communicator                                            \
        )                                                                     \
        {                                                                     \
            UPstream::mpiScan<OpCode>(values, count, communicator, true);     \
        }                                                                     \
                                                                              \
        /*!\brief Inclusive/exclusive \b \c OpName scan returning result */   \
        /*! See UPstream::mpiScan for notes about the behaviour */            \
        template<class T>                                                     \
        [[nodiscard]]                                                         \
        static T mpiScan_##OpName                                             \
        (                                                                     \
            const T& value,                                                   \
            const int communicator,                                           \
            /*! Use exclusive scan (MPI_Exscan vs MPI_Scan) */                \
            const bool exclusive = false                                      \
        )                                                                     \
        {                                                                     \
            return UPstream::mpiScan<OpCode>(value, communicator, exclusive); \
        }                                                                     \
                                                                              \
        /*!\brief Exclusive \b \c OpName scan returning result */             \
        /*! See UPstream::mpiScan for notes about the behaviour */            \
        template<class T>                                                     \
        [[nodiscard]]                                                         \
        static T mpiExscan_##OpName                                           \
        (                                                                     \
            const T& value,                                                   \
            const int communicator                                            \
        )                                                                     \
        {                                                                     \
            return UPstream::mpiScan<OpCode>(value, communicator, true);      \
        }

        Pstream_declareMethods(min, UPstream::opCodes::op_min);
        Pstream_declareMethods(max, UPstream::opCodes::op_max);
        Pstream_declareMethods(sum, UPstream::opCodes::op_sum);
        #undef Pstream_declareMethods


    // Logical reductions

        //- Logical (and) reduction (MPI_AllReduce)
        //  For \b non-parallel : do nothing
        static void reduceAnd
        (
            bool& value,
            const int communicator = UPstream::worldComm
        );

        //- Logical (or) reduction (MPI_AllReduce)
        //  For \b non-parallel : do nothing
        static void reduceOr
        (
            bool& value,
            const int communicator = UPstream::worldComm
        );

        //- Locate the first rank for which the \p condition is true,
        //- or -1 if no ranks satisfy the condition.
        static int find_first(bool condition, int communicator);

        //- Locate the last rank for which the \p condition is true,
        //- or -1 if no ranks satisfy the condition.
        static int find_last(bool condition, int communicator);


    // Housekeeping

        //- \deprecated(2025-02) prefer newCommunicator
        FOAM_DEPRECATED_FOR(2025-02, "newCommunicator()")
        static label allocateCommunicator
        (
            const label parent,
            const labelRange& subRanks,
            const bool withComponents = true
        )
        {
            return newCommunicator(parent, subRanks, withComponents);
        }

        //- \deprecated(2025-02) prefer newCommunicator
        FOAM_DEPRECATED_FOR(2025-02, "newCommunicator()")
        static label allocateCommunicator
        (
            const label parent,
            const labelUList& subRanks,
            const bool withComponents = true
        )
        {
            return newCommunicator(parent, subRanks, withComponents);
        }

        //- Wait for all requests to finish.
        //  \deprecated(2023-01) Probably not what you want.
        //     Should normally be restricted to a particular starting request.
        FOAM_DEPRECATED_FOR(2023-01, "waitRequests(int) method")
        static void waitRequests() { waitRequests(0); }

        //- \deprecated(2025-02) prefer mpiGatherv
        template<class Type>
        FOAM_DEPRECATED_FOR(2025-02, "mpiGatherv()")
        static void gather
        (
            const Type* send,
            int count,
            Type* recv,
            const UList<int>& counts,
            const UList<int>& offsets,
            const int comm = UPstream::worldComm
        )
        {
            UPstream::mpiGatherv(send, count, recv, counts, offsets, comm);
        }

        //- \deprecated(2025-02) prefer mpiScatterv
        template<class Type>
        FOAM_DEPRECATED_FOR(2025-02, "mpiScatterv()")
        static void scatter
        (
            const Type* send,
            const UList<int>& counts,
            const UList<int>& offsets,
            Type* recv,
            int count,
            const int comm = UPstream::worldComm
        )
        {
            UPstream::mpiScatterv(send, counts, offsets, recv, count, comm);
        }
};

/*---------------------------------------------------------------------------*\
                   Class UPstream::communicator Declaration
\*---------------------------------------------------------------------------*/

//- Wrapper for internally indexed communicator label.
//- Always invokes UPstream::allocateCommunicatorComponents()
//- and UPstream::freeCommunicatorComponents()
class UPstream::communicator
{
    // Private Data

        //- The communicator label
        label comm_;

public:

    // Generated Methods

        //- No copy construct
        communicator(const communicator&) = delete;

        //- No copy assignment
        void operator=(const communicator&) = delete;


    // Constructors

        //- Default construct (a placeholder communicator)
        communicator() noexcept : comm_(-1) {}

        //- Move construct, takes ownership
        communicator(communicator&& c) : comm_(c.comm_) { c.comm_ = -1; }

        //- Allocate communicator for contiguous sub-ranks on given parent
        communicator
        (
            //! The parent communicator
            const label parentComm,
            //! The contiguous sub-ranks of parent to use
            const labelRange& subRanks
        )
        :
            comm_(UPstream::newCommunicator(parentComm, subRanks))
        {}

        //- Allocate communicator for sub-ranks on given parent
        communicator
        (
            //! The parent communicator
            const label parentComm,
            //! The sub-ranks of parent to use (negative values ignored)
            const labelUList& subRanks
        )
        :
            comm_(UPstream::newCommunicator(parentComm, subRanks))
        {}


    // Destructor

        //- Free allocated communicator
        ~communicator() { UPstream::freeCommunicator(comm_); }


    // Factory Methods

        //- Duplicate the given communicator
        [[nodiscard]]
        static communicator duplicate(const label parentComm)
        {
            communicator c;
            c.comm_ = UPstream::dupCommunicator(parentComm);
            return c;
        }

        //- Factory Method :
        //- Split the communicator on the given \em colour.
        [[nodiscard]]
        static communicator split
        (
            //! The parent communicator
            const label parentComm,
            //! The colouring to select which ranks to include.
            //! Negative values correspond to 'ignore'
            const int colour,
            //! Use MPI_Allgather+MPI_Comm_create_group vs MPI_Comm_split
            const bool two_step = true
        )
        {
            communicator c;
            c.comm_ = UPstream::splitCommunicator(parentComm, colour, two_step);
            return c;
        }


    // Member Functions

        //- True if communicator is non-negative (ie, was allocated)
        bool good() const noexcept { return (comm_ >= 0); }

        //- The communicator label
        label comm() const noexcept { return comm_; }

        //- Return non-const reference to this
        communicator& constCast() const noexcept
        {
            return const_cast<communicator&>(*this);
        }

        //- Release ownership of the communicator, return old value.
        //  Leave further management to the caller
        label release() noexcept { label c(comm_); comm_ = -1; return c; }

        //- Free allocated communicator
        void reset() { UPstream::freeCommunicator(comm_); comm_ = -1; }

        //- Allocate with contiguous sub-ranks of parent communicator
        void reset(label parent, const labelRange& subRanks)
        {
            UPstream::freeCommunicator(comm_);
            comm_ = UPstream::newCommunicator(parent, subRanks);
        }

        //- Allocate with sub-ranks of parent communicator
        void reset(label parent, const labelUList& subRanks)
        {
            UPstream::freeCommunicator(comm_);
            comm_ = UPstream::newCommunicator(parent, subRanks);
        }

        //- Take ownership, free managed communicator
        void reset(communicator&& c)
        {
            if (this == &c) return;  // No self-assignment
            if (comm_ != c.comm_) UPstream::freeCommunicator(comm_);
            comm_ = c.comm_;
            c.comm_ = -1;
        }

        //- Swap communicator labels
        void swap(communicator& c) { std::swap(comm_, c.comm_); }


    // Member Operators

        //- Implicit cast to label - the same as comm()
        operator label() const noexcept { return comm_; }

        //- Move assignment, takes ownership
        void operator=(communicator&& c) { reset(std::move(c)); }

        //- Test for equality
        bool operator==(const communicator& c) const noexcept
        {
            return (comm_ == c.comm_);
        }

        //- Test for inequality
        bool operator!=(const communicator& c) const noexcept
        {
            return (comm_ != c.comm_);
        }
};


/*---------------------------------------------------------------------------*\
                   Class UPstream::Communicator Declaration
\*---------------------------------------------------------------------------*/

//- An opaque wrapper for MPI_Comm with a vendor-independent
//- representation without any \c <mpi.h> header.
//  The MPI standard states that MPI_Comm is always an opaque object.
//  Generally it is either an integer (eg, mpich) or a pointer (eg, openmpi).
class UPstream::Communicator
{
public:

    // Public Types

        //- Storage for MPI_Comm (as integer or pointer)
        typedef std::intptr_t value_type;


private:

    // Private Data

        //- The MPI_Comm (as wrapped value)
        value_type value_;

public:

    // Generated Methods

        //- Copy construct
        Communicator(const Communicator&) noexcept = default;

        //- Move construct
        Communicator(Communicator&&) noexcept = default;

        //- Copy assignment
        Communicator& operator=(const Communicator&) noexcept = default;

        //- Move assignment
        Communicator& operator=(Communicator&&) noexcept = default;


    // Member Operators

        //- Test for equality
        bool operator==(const Communicator& rhs) const noexcept
        {
            return (value_ == rhs.value_);
        }

        //- Test for inequality
        bool operator!=(const Communicator& rhs) const noexcept
        {
            return (value_ != rhs.value_);
        }


    // Constructors

        //- Default construct as MPI_COMM_NULL
        Communicator() noexcept;

        //- Construct from MPI_Comm (as pointer type)
        explicit Communicator(const void* p) noexcept
        :
            value_(reinterpret_cast<value_type>(p))
        {}

        //- Construct from MPI_Comm (as integer type)
        explicit Communicator(value_type val) noexcept
        :
            value_(val)
        {}


    // Factory Methods

        //- Transcribe internally indexed communicator to wrapped value.
        //
        //  Example,
        //  \code
        //  PstreamUtils::Cast::to_mpi
        //  (
        //      UPstream::Communicator::lookup(UPstream::worldComm)
        //  )
        //  \endcode
        static Communicator lookup
        (
            //! The internally indexed communicator.
            //! The default (-1) is the same as UPstream::worldComm.
            int communicator = -1
        );


    // Member Functions

    // Access

        //- Return raw value
        value_type value() const noexcept { return value_; }

        //- Return as pointer value
        const void* pointer() const noexcept
        {
            return reinterpret_cast<const void*>(value_);
        }


    // Query

        //- True if not equal to MPI_COMM_NULL
        bool good() const noexcept;

        //- Reset to default constructed value (MPI_COMM_NULL)
        void reset() noexcept;

        //- The number of ranks associated with the communicator
        int size() const;
};


/*---------------------------------------------------------------------------*\
                      Class UPstream::Request Declaration
\*---------------------------------------------------------------------------*/

//- An opaque wrapper for MPI_Request with a vendor-independent
//- representation without any \c <mpi.h> header.
//  The MPI standard states that MPI_Request is always an opaque object.
//  Generally it is either an integer (eg, mpich) or a pointer (eg, openmpi).
class UPstream::Request
{
public:

    // Public Types

        //- Storage for MPI_Request (as integer or pointer)
        typedef std::intptr_t value_type;


private:

    // Private Data

        //- The MPI_Request (as wrapped value)
        value_type value_;

public:

    // Generated Methods

        //- Copy construct
        Request(const Request&) noexcept = default;

        //- Move construct
        Request(Request&&) noexcept = default;

        //- Copy assignment
        Request& operator=(const Request&) noexcept = default;

        //- Move assignment
        Request& operator=(Request&&) noexcept = default;


    // Member Operators

        //- Test for equality
        bool operator==(const Request& rhs) const noexcept
        {
            return (value_ == rhs.value_);
        }

        //- Test for inequality
        bool operator!=(const Request& rhs) const noexcept
        {
            return (value_ != rhs.value_);
        }


    // Constructors

        //- Default construct as MPI_REQUEST_NULL
        Request() noexcept;

        //- Construct from MPI_Request (as pointer type)
        explicit Request(const void* p) noexcept
        :
            value_(reinterpret_cast<value_type>(p))
        {}

        //- Construct from MPI_Request (as integer type)
        explicit Request(value_type val) noexcept
        :
            value_(val)
        {}


    // Member Functions

    // Access

        //- Return raw value
        value_type value() const noexcept { return value_; }

        //- Return as pointer value
        const void* pointer() const noexcept
        {
            return reinterpret_cast<const void*>(value_);
        }


    // Basics

        //- True if not equal to MPI_REQUEST_NULL
        bool good() const noexcept;

        //- Reset to default constructed value (MPI_REQUEST_NULL)
        void reset() noexcept;

        //- True if request is active (!= MPI_REQUEST_NULL)
        //- Same as good(). Same as calling UPstream::activeRequest()
        bool active() const noexcept { return good(); }

        //- Same as calling UPstream::finishedRequest()
        //  Uses MPI_Test()
        bool finished() { return UPstream::finishedRequest(*this); }

        //- Same as calling UPstream::waitRequest().
        //  Uses MPI_Wait()
        void wait() { UPstream::waitRequest(*this); }


    // Other

        //- Same as calling UPstream::cancelRequest().
        //  Uses MPI_Cancel(), MPI_Request_free()
        void cancel() { UPstream::cancelRequest(*this); }

        //- Same as calling UPstream::freeRequest().
        //  Uses MPI_Request_free()
        void free() { UPstream::freeRequest(*this); }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream&, const UPstream::commsStruct&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Traits, nested classes etc
#include "UPstreamTraits.H"
#include "UPstreamWindow.H"

// Locally include the following where required:
//   - UPstreamFile.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPstream.txx"
    #include "UPstreamReduceOffsets.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
