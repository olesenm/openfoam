/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
    Foam::Pstream

Description
    Inter-processor communications stream.

SourceFiles
    Pstream.C
    PstreamBroadcast.txx
    PstreamGather.txx
    PstreamGatherList.txx
    PstreamGatherValues.txx
    PstreamExchange.txx
    PstreamExchangeConsensus.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_Pstream_H
#define Foam_Pstream_H

#include "UPstream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;

/*---------------------------------------------------------------------------*\
                           Class Pstream Declaration
\*---------------------------------------------------------------------------*/

class Pstream
:
    public UPstream
{
protected:

    // Protected Data

        //- Allocated transfer buffer (can be used for send or receive)
        DynamicList<char> transferBuf_;


public:

    //- Declare name of the class and its debug switch
    ClassName("Pstream");


    // Constructors

        //- Construct for communication type with empty buffer
        explicit Pstream(UPstream::commsTypes commsType) noexcept
        :
            UPstream(commsType)
        {}

        //- Construct for communication type with given buffer size
        Pstream(UPstream::commsTypes commsType, int bufferSize)
        :
            UPstream(commsType)
        {
            if (bufferSize > 0)
            {
                transferBuf_.setCapacity(bufferSize + 2*sizeof(double) + 1);
            }
        }


    // Static Functions

    // Broadcast

        //- Broadcast buffer content to all processes in communicator.
        using UPstream::broadcast;

        //- Broadcast content (contiguous or non-contiguous) to all
        //- communicator ranks. Does nothing in \b non-parallel.
        template<class Type>
        static void broadcast
        (
            Type& value,
            const int communicator = UPstream::worldComm
        );

        //- Broadcast fixed-list content (contiguous or non-contiguous) to all
        //- communicator ranks. Does nothing in \b non-parallel.
        template<class Type, unsigned N>
        static void broadcast
        (
            FixedList<Type, N>& list,
            const int communicator = UPstream::worldComm
        );

        //- Broadcast multiple items to all communicator ranks.
        //- Does nothing in \b non-parallel.
        template<class Type, class... Args>
        static void broadcasts
        (
            const int communicator,
            Type& value,
            Args&&... values
        );

        //- Broadcast list content (contiguous or non-contiguous) to all
        //- communicator ranks. Does nothing in \b non-parallel.
        //  For contiguous list data, this avoids serialization overhead,
        //  but at the expense of an additional broadcast call.
        template<class ListType>
        static void broadcastList
        (
            ListType& list,
            const int communicator = UPstream::worldComm
        );


    // Gather/scatter : single value

        //- Gather (reduce) data, applying \c bop to combine \c value
        //- from different processors. The basis for Foam::reduce().
        //  A no-op for non-parallel.
        //
        //  \tparam InplaceMode indicates that the binary operator
        //      modifies values in-place, not using assignment
        template<class T, class BinaryOp, bool InplaceMode=false>
        static void gather
        (
            //! [in,out] the result is only reliable on rank=0
            T& value,
            BinaryOp bop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Gather individual values into list locations.
        //  On master list length == nProcs, otherwise zero length.
        //  \n
        //  For \b non-parallel :
        //  the returned list length is 1 with localValue.
        template<class T>
        static List<T> listGatherValues
        (
            const T& localValue,
            const int communicator = UPstream::worldComm,
            //! Only used for non-contiguous types
            const int tag = UPstream::msgType()
        );

        //- Scatter individual values from list locations.
        //  On master input list length == nProcs, ignored on other procs.
        //  \n
        //  For \b non-parallel :
        //  returns the first list element (or default initialized).
        template<class T>
        static T listScatterValues
        (
            const UList<T>& allValues,
            const int communicator = UPstream::worldComm,
            //! Only used for non-contiguous types
            const int tag = UPstream::msgType()
        );


    // Inplace combine (gather) : single value

        //- Forwards to Pstream::gather with an \em in-place \c cop
        template<class T, class CombineOp>
        static void combineGather
        (
            //! [in,out] the result is only reliable on rank=0
            T& value,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Reduce inplace (cf. MPI Allreduce)
        //- applying \c cop to inplace combine \c value
        //- from different processors.
        template<class T, class CombineOp>
        static void combineReduce
        (
            //! [in,out] the result is consistent on all ranks
            T& value,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Same as Pstream::combineReduce
        template<class T, class CombineOp>
        static void combineAllGather
        (
            //! [in,out] the result is consistent on all ranks
            T& value,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        )
        {
            Pstream::combineReduce(value, cop, tag, communicator);
        }


    // Gather/combine variants working on entire List

        //- Gather (reduce) list elements,
        //- applying \c bop to each list element
        //  Respects UPstream::listChunkLength()
        //
        //  \tparam InplaceMode indicates that the binary operator
        //      modifies values in-place, not using assignment
        template<class T, class BinaryOp, bool InplaceMode=false>
        static void listGather
        (
            //! [in,out] the result is only reliable on rank=0
            UList<T>& values,
            BinaryOp bop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Forwards to Pstream::listGather with an \em in-place \c cop
        //  Respects UPstream::listChunkLength()
        template<class T, class CombineOp>
        static void listCombineGather
        (
            //! [in,out] the result is only reliable on rank=0
            UList<T>& values,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Reduce list elements (list must be equal size on all ranks),
        //- applying \c bop to each list element.
        //  Respects UPstream::listChunkLength()
        //
        //  \tparam InplaceMode indicates that the binary operator
        //      modifies values in-place, not using assignment
        template<class T, class BinaryOp, bool InplaceMode=false>
        static void listReduce
        (
            //! [in,out] the result is consistent on all ranks
            UList<T>& values,
            BinaryOp bop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Forwards to Pstream::listReduce with an \em in-place \c cop
        //  Respects UPstream::listChunkLength()
        template<class T, class CombineOp>
        static void listCombineReduce
        (
            //! [in,out] the result is consistent on all ranks
            UList<T>& values,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );


    // Gather/combine variants working on Map/HashTable containers

        //- Gather (reduce) Map/HashTable containers,
        //- applying \c bop to combine entries from different processors.
        //
        //  \tparam InplaceMode indicates that the binary operator
        //      modifies values in-place, not using assignment
        template<class Container, class BinaryOp, bool InplaceMode=false>
        static void mapGather
        (
            //! [in,out] the result is only reliable on rank=0
            Container& values,
            BinaryOp bop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Forwards to Pstream::mapGather with an \em in-place \c cop
        template<class Container, class CombineOp>
        static void mapCombineGather
        (
            //! [in,out] the result is only reliable on rank=0
            Container& values,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Reduce inplace (cf. MPI Allreduce)
        //- applying \c bop to combine map \c values
        //- from different processors.
        //- After completion all processors have the same data.
        //
        //  Wraps mapCombineGather/broadcast (may change in the future).
        template<class Container, class BinaryOp, bool InplaceMode=false>
        static void mapReduce
        (
            //! [in,out] the result is consistent on all ranks
            Container& values,
            BinaryOp bop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Forwards to Pstream::mapReduce with an \em in-place \c cop
        template<class Container, class CombineOp>
        static void mapCombineReduce
        (
            //! [in,out] the result is consistent on all ranks
            Container& values,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Same as Pstream::mapCombineReduce
        template<class Container, class CombineOp>
        static void mapCombineAllGather
        (
            //! [in,out] the result is consistent on all ranks
            Container& values,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        )
        {
            Pstream::mapCombineReduce(values, cop, tag, communicator);
        }


    // Gather/scatter keeping the individual processor data separate.
    // The values is a List of size UPstream::nProcs() where
    // values[UPstream::myProcNo()] is the data for the current processor.

        //- Gather data, but keep individual values separate.
        template<class T>
        static void gatherList
        (
            //! [in,out]
            UList<T>& values,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );

        //- Gather data, but keep individual values separate.
        //- Uses MPI_Allgather or manual communication.
        //  After completion all processors have the same data.
        //  Wraps gatherList/scatterList (may change in the future).
        template<class T>
        static void allGatherList
        (
            //! [in,out]
            UList<T>& values,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );


    // Exchange

        //- Helper: exchange sizes of sendBufs for specified send/recv ranks
        template<class Container>
        static void exchangeSizes
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            const Container& sendBufs,
            labelList& sizes,
            const int tag = UPstream::msgType(),
            const int comm = UPstream::worldComm
        );

        //- Helper: exchange sizes of sendBufs for specified neighbour ranks
        template<class Container>
        static void exchangeSizes
        (
            const labelUList& neighProcs,
            const Container& sendBufs,
            labelList& sizes,
            const int tag = UPstream::msgType(),
            const int comm = UPstream::worldComm
        );

        //- Helper: exchange sizes of sendBufs.
        //- The sendBufs is the data per processor (in the communicator).
        //  Returns sizes of sendBufs on the sending processor.
        //  \n
        //  For \b non-parallel : copy sizes from sendBufs directly.
        template<class Container>
        static void exchangeSizes
        (
            const Container& sendBufs,
            labelList& recvSizes,
            const int comm = UPstream::worldComm
        );

        //- Exchange the \b non-zero sizes of sendBufs entries (sparse map)
        //- with other ranks in the communicator
        //- using non-blocking consensus exchange.
        //
        //  Since the recvData map always cleared before receipt and sizes
        //  of zero are never transmitted, a simple check
        //  of its keys is sufficient to determine connectivity.
        //
        //  For \b non-parallel : copy size of rank (if it exists and non-empty)
        //  from sendBufs to recvSizes.
        //
        //  \note The message tag is adjusted internally to improve uniqueness
        template<class Container>
        static void exchangeSizes
        (
            const Map<Container>& sendBufs,
            Map<label>& recvSizes,
            const int tag = UPstream::msgType(),
            const int comm = UPstream::worldComm
        );

        //- Helper: exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs using predetermined receive
        //- sizing.
        //  If wait=true will wait for all transfers to finish.
        template<class Container, class Type>
        static void exchange
        (
            const UList<Container>& sendBufs,
            const labelUList& recvSizes,
            List<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const int comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );

        //- Exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs.
        //  Data provided and received as container.
        //
        // No internal guards or resizing.
        template<class Container, class Type>
        static void exchange
        (
            const Map<Container>& sendBufs,
            const Map<label>& recvSizes,  //!< Num of recv elements (not bytes)
            Map<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const int comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );

        //- Exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs.
        //- Determines sizes to receive.
        //  If wait=true will wait for all transfers to finish.
        template<class Container, class Type>
        static void exchange
        (
            const UList<Container>& sendBufs,
            List<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const int comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );

        //- Exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs.
        //- Determines sizes to receive.
        //  If wait=true will wait for all transfers to finish.
        template<class Container, class Type>
        static void exchange
        (
            const Map<Container>& sendBufs,
            Map<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const int comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );


    // Non-blocking exchange

        //- Exchange \em contiguous data using non-blocking consensus (NBX)
        //- Sends sendData, receives into recvData.
        //
        //  Each entry of the recvBufs list is cleared before receipt.
        //  For \b non-parallel : copy own rank from sendBufs to recvBufs.
        //
        //  \note The message tag should be chosen to be a unique value
        //  since the implementation uses probing with ANY_SOURCE !!
        template<class Container, class Type>
        static void exchangeConsensus
        (
            const UList<Container>& sendBufs,
            List<Container>& recvBufs,
            const int tag,
            const int comm,
            const bool wait = true  //!< (ignored)
        );

        //- Exchange \em contiguous data using non-blocking consensus (NBX)
        //- Sends sendData, receives into recvData.
        //
        //  Each \em entry of the recvBufs map is cleared before receipt,
        //  but the map itself if not cleared. This allows the map to preserve
        //  allocated space (eg DynamicList entries) between calls.
        //
        //  For \b non-parallel : copy own rank (if it exists and non-empty)
        //  from sendBufs to recvBufs.
        //
        //  \note The message tag should be chosen to be a unique value
        //  since the implementation uses probing with ANY_SOURCE !!
        template<class Container, class Type>
        static void exchangeConsensus
        (
            const Map<Container>& sendBufs,
            Map<Container>& recvBufs,
            const int tag,
            const int comm,
            const bool wait = true  //!< (ignored)
        );

        //- Exchange \em contiguous data using non-blocking consensus (NBX)
        //- Sends sendData returns receive information.
        //
        //  For \b non-parallel : copy own rank (if it exists and non-empty)
        //
        //  \note The message tag should be chosen to be a unique value
        //  since the implementation uses probing with ANY_SOURCE !!
        template<class Container, class Type>
        static Map<Container> exchangeConsensus
        (
            const Map<Container>& sendBufs,
            const int tag,
            const int comm,
            const bool wait = true  //!< (ignored)
        );


    // Housekeeping

        //- \deprecated(2024-01) Broadcast data
        template<class T>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void scatter
        (
            T& value,
            const int tag = UPstream::msgType(),  //!< ignored
            const int comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(value, comm);
        }

        //- \deprecated(2024-01) Broadcast data
        template<class T>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void combineScatter
        (
            T& value,
            const int tag = UPstream::msgType(),  //!< ignored
            const int comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(value, comm);
        }

        //- \deprecated(2024-01) Broadcast data
        template<class T>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void listCombineScatter
        (
            List<T>& value,
            const int tag = UPstream::msgType(),  //!< ignored
            const int comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(value, comm);
        }

        //- \deprecated(2022-08) == listCombineReduce
        template<class T, class CombineOp>
        FOAM_DEPRECATED_FOR(2022-08, "Pstream::listCombineReduce()")
        static void listCombineAllGather
        (
            //! [in,out] the result is consistent on all ranks
            UList<T>& values,
            CombineOp cop,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        )
        {
            Pstream::listCombineReduce(values, cop, tag, communicator);
        }

        //- \deprecated(2024-01) Broadcast data
        template<class Container>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void mapCombineScatter
        (
            Container& values,
            const int tag = UPstream::msgType(),  //!< ignored
            const int comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(values, comm);
        }

        //- The inverse of gatherList, but when combined with gatherList
        //- it effectively acts like a partial broadcast...
        //  \deprecated(2025-03) Normally prefer broadcast
        template<class T>
        FOAM_DEPRECATED_FOR(2025-03, "broadcast() or broadcastList()")
        static void scatterList
        (
            //! [in,out]
            UList<T>& values,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PstreamBroadcast.txx"
    #include "PstreamGather.txx"
    #include "PstreamGatherList.txx"
    #include "PstreamGatherValues.txx"
    #include "PstreamExchange.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
