/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::UPstream::broadcast
(
    Type* buffer,
    std::streamsize count,
    const int communicator,
    const int root
)
{
    // Likely no reason to check for nullptr
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return true;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        // Also report parameters to silence compiler warnings about unused
        FatalErrorInFunction
            << "Invalid for non-contiguous data types."
            << " buffer:" << (buffer != nullptr)
            << " count:" << count
            << " root:" << root
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        return UPstream::mpi_broadcast
        (
            buffer,     // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            communicator,
            root
        );
    }
}


template<class Type, unsigned N>
bool Foam::UPstream::broadcast
(
    FixedList<Type, N>& list,
    const int communicator,
    const int root
)
{
    return UPstream::broadcast(list.data(), list.size(), communicator, root);
}


template<class Type>
void Foam::UPstream::mpiGather
(
    const Type* sendData,
    Type* recvData,
    int count,
    const int communicator
)
{
    if (!count || !UPstream::is_rank(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }
    else if (!UPstream::is_parallel(communicator))
    {
        // Perform any fallback copying here, while we still know the Type
        if (sendData && recvData && (sendData != recvData))
        {
            std::copy_n(sendData, count, recvData);
        }
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        UPstream::mpi_gather
        (
            sendData,   // The data or cmpt pointer
            recvData,   // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
}


template<class Type>
void Foam::UPstream::mpiScatter
(
    const Type* sendData,
    Type* recvData,
    int count,
    const int communicator
)
{
    if (!count || !UPstream::is_rank(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }
    else if (!UPstream::is_parallel(communicator))
    {
        // Perform any fallback copying here, while we still know the Type
        if (sendData && recvData && (sendData != recvData))
        {
            std::copy_n(sendData, count, recvData);
        }
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        UPstream::mpi_scatter
        (
            sendData,   // The data or cmpt pointer
            recvData,   // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
}


template<class Type>
void Foam::UPstream::mpiAllGather
(
    Type* allData,
    int count,
    const int communicator
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing sensible to do
        return;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        UPstream::mpi_allgather
        (
            allData,    // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
}


template<class T>
Foam::List<T> Foam::UPstream::allGatherValues
(
    const T& localValue,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(comm) as well?
        List<T> allValues(1);
        allValues[0] = localValue;
        return allValues;
    }
    else if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot all-gather values for non-contiguous types"
            << " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);
        return List<T>();
    }
    else
    {
        // Standard gather with direct MPI communication
        List<T> allValues;

        allValues.resize(UPstream::nProcs(communicator));
        allValues[UPstream::myProcNo(communicator)] = localValue;

        UPstream::mpiAllGather(allValues.data(), 1, communicator);
        return allValues;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T> Foam::UPstream::listGatherValues
(
    const T& localValue,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(communicator) as well?
        List<T> allValues(1);
        allValues[0] = localValue;
        return allValues;
    }
    else if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot gather values for non-contiguous types"
               " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);
        return List<T>();
    }
    else
    {
        // Local sizes are identical, can use MPI_Gather
        List<T> allValues;

        if (UPstream::master(communicator))
        {
            allValues.resize(UPstream::nProcs(communicator));
        }

        UPstream::mpiGather(&localValue, allValues.data(), 1, communicator);
        return allValues;
    }
}


template<class T>
T Foam::UPstream::listScatterValues
(
    const UList<T>& allValues,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(communicator) as well?

        if (!allValues.empty())
        {
            return allValues[0];
        }

        return T{};  // Fallback value
    }
    else if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot scatter non-contiguous values"
               " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);

        return T{};  // Fallback value
    }
    else
    {
        // Local sizes are identical, can use MPI_Scatter

        const label nProcs = UPstream::nProcs(communicator);

        if
        (
            FOAM_UNLIKELY
            (
                UPstream::master(communicator)
             && allValues.size() < nProcs
            )
        )
        {
            FatalErrorInFunction
                << "Attempting to send " << allValues.size()
                << " values to " << nProcs << " processors" << endl
                << Foam::abort(FatalError);
        }

        T localValue{};
        UPstream::mpiScatter(allValues.cdata(), &localValue, 1, communicator);
        return localValue;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::UPstream::mpiGatherv
(
    const Type* sendData,
    int sendCount,
    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        if constexpr (is_contiguous_v<Type>)
        {
            if (sendData && recvData && (sendData != recvData))
            {
                // recvCounts[0] may be invalid - use sendCount instead
                std::copy_n(sendData, sendCount, recvData);
            }
        }
        // Nothing further to do
    }
    else if constexpr (UPstream_dataType<Type>::value)
    {
        // Restrict to basic (or aliased) MPI types to avoid recalculating
        // the list of counts/offsets.

        // The sizing factor (constexpr) must be 1 otherwise
        // [recvCounts,recvOffsets] are likely incorrect

        constexpr std::streamsize count = UPstream_dataType<Type>::size(1);

        static_assert
        (
            (count == 1),
            "Code does not (yet) work with aggregate types"
        );

        UPstream::mpi_gatherv
        (
            sendData,
            sendCount,
            recvData,
            recvCounts,
            recvOffsets,

            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
    else
    {
        static_assert
        (
            stdFoam::dependent_false_v<Type>,
            "Only basic and user data types"
        );
    }
}


template<class Type>
void Foam::UPstream::mpiScatterv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,
    Type* recvData,
    int recvCount,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        if constexpr (is_contiguous_v<Type>)
        {
            if (sendData && recvData && (sendData != recvData))
            {
                std::copy_n(sendData, recvCount, recvData);
            }
        }
        // Nothing further to do
    }
    else if constexpr (UPstream_dataType<Type>::value)
    {
        // Restrict to basic (or aliased) MPI types to avoid recalculating
        // the list of counts/offsets.

        // The sizing factor (constexpr) must be 1 otherwise
        // [sendCounts,sendOffsets] are likely incorrect

        constexpr std::streamsize count = UPstream_dataType<Type>::size(1);

        static_assert
        (
            (count == 1),
            "Code does not (yet) work with aggregate types"
        );

        UPstream::mpi_scatterv
        (
            sendData,
            sendCounts,
            sendOffsets,
            recvData,
            recvCount,

            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
    else
    {
        static_assert
        (
            stdFoam::dependent_false_v<Type>,
            "Only basic and user data types"
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void Foam::UPstream::mpiReduce
(
    T values[],
    int count,
    const UPstream::opCodes opCodeId,
    const int communicator
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        // Restricted to basic data types
        UPstream::mpi_reduce
        (
            values,     // The data or cmpt pointer
            UPstream_basic_dataType<T>::size(count),
            UPstream_basic_dataType<T>::datatype_id,
            opCodeId,
            communicator
        );
    }
}


template<Foam::UPstream::opCodes opCode, class T>
void Foam::UPstream::mpiReduce
(
    T values[],
    int count,
    const int communicator
)
{
    UPstream::mpiReduce(values, count, opCode, communicator);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void Foam::UPstream::mpiAllReduce
(
    T values[],
    int count,
    const UPstream::opCodes opCodeId,
    const int communicator
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        // Restricted to basic data types
        UPstream::mpi_allreduce
        (
            values,     // The data or cmpt pointer
            UPstream_basic_dataType<T>::size(count),
            UPstream_basic_dataType<T>::datatype_id,
            opCodeId,
            communicator
        );
    }
}


template<Foam::UPstream::opCodes opCode, class T>
void Foam::UPstream::mpiAllReduce
(
    T values[],
    int count,
    const int communicator
)
{
    UPstream::mpiAllReduce(values, count, opCode, communicator);
}


template<class T>
void Foam::UPstream::mpiAllReduce
(
    T values[],
    int count,
    const UPstream::opCodes opCodeId,
    const int communicator,
    UPstream::Request& req
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        // Restricted to basic data types
        UPstream::mpi_allreduce
        (
            values,     // The data or cmpt pointer
            UPstream_basic_dataType<T>::size(count),
            UPstream_basic_dataType<T>::datatype_id,
            opCodeId,
            communicator,
           &req
        );
    }
}


template<Foam::UPstream::opCodes opCode, class T>
void Foam::UPstream::mpiAllReduce
(
    T values[],
    int count,
    const int communicator,
    UPstream::Request& req
)
{
    UPstream::mpiAllReduce(values, count, opCode, communicator);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<Foam::UPstream::opCodes opCode, bool exclusive, class T>
void Foam::UPstream::mpiScan
(
    T values[],
    int count,
    const int communicator
)
{
    if (count <= 0)
    {
        // Nothing to do
        return;
    }

    if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot " << (exclusive ? "exscan" : "scan")
            << " values for non-contiguous types."
            << " op-code: " << int(opCode)
            << " buffer: " << Foam::name(values) << endl
            << Foam::abort(FatalError);
    }
    else if (UPstream::is_parallel(communicator))  // && FOAM_LIKELY(values)
    {
        // Use element or component type (or byte-wise) for data type
        // Restricted to basic data types
        UPstream::mpi_scan_reduce
        (
            values,     // The data or cmpt pointer
            UPstream_basic_dataType<T>::size(count),
            UPstream_basic_dataType<T>::datatype_id,
            opCode,
            communicator,
            exclusive
        );
    }

    // Any special handling of exscan degenerate value at rank=0

    if constexpr (exclusive && UPstream::opCodes::op_sum == opCode)
    {
        // Exscan and sum:
        // - replace degenerate value (rank=0) with zero value
        //   even if non-parallel

        if (!UPstream::parRun() || UPstream::master(communicator))
        {
            T zeroValue{};
            if constexpr (is_vectorspace_v<T>)
            {
                zeroValue.fill(0);
            }
            std::fill_n(values, count, zeroValue);
        }
    }
}


template<Foam::UPstream::opCodes opCode, bool exclusive, class T>
T Foam::UPstream::mpiScan
(
    const T& localValue,
    int communicator
)
{
    T work(localValue);
    UPstream::mpiScan<opCode, exclusive>(&work, 1, communicator);
    return work;
}


template<class IntType>
IntType Foam::UPstream::find_offset
(
    IntType localSize,
    int communicator
)
{
    static_assert(std::is_integral_v<IntType>, "Integral required");

    // Exclusive scan sum
    IntType work(localSize);
    UPstream::mpiScan<UPstream::opCodes::op_sum, true>(&work, 1, communicator);

    // Always return offset=0 (no offset) for master rank or non-parallel.
    // Our mpiExscan_sum() should already do that [implemented above],
    // but force it again for good measure.

    return
    (
        (UPstream::parRun() && UPstream::is_subrank(communicator)) ? work : 0
    );
}


// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// Dispatch for allToAll_n or allToAllConsensus_n etc.
// The Container is FixedList<T,Count> or Pair<T> etc.
//
// The AllToAllOp is invoked as op(sends, recvs, Count) with the flattened
// lists, and performs the actual communication.
template<class Container, class T, int Count, class AllToAllOp>
void allToAll_fixedWidth
(
    const UList<Container>& sendData,
    UList<Container>& recvData,
    const int communicator,
    AllToAllOp&& op
)
{
    // Restrict to known and supported integral types.
    // This could be improved/expanded in the future if needed

    if constexpr
    (
        std::is_same_v<int32_t, T>
     || std::is_same_v<int64_t, T>
    )
    {
        // Not involved
        if (!UPstream::is_rank(communicator))
        {
            recvData = Foam::zero{};
            return;
        }
        else if (!UPstream::is_parallel(communicator))
        {
            recvData.deepCopy(sendData);
            return;
        }

        // Wrap content into the flattened UList representation
        // used by the base-level operation.
        // Note that pointer to Container content is also the pointer
        // to its 0th element for these fixed-width types.
        // The const_cast is only due to the UList API requirements

        const UList<T> sends
        (
            const_cast<T*>
            (reinterpret_cast<const T*>(sendData.cdata())),
            (Count*sendData.size())
        );
        UList<T> recvs
        (
            reinterpret_cast<T*>(recvData.data()),
            (Count*recvData.size())
        );

        op(sends, recvs, Count);
    }
    else
    {
        static_assert
        (
            stdFoam::dependent_false_v<T>,
            "Currently restricted to integral types (int32|int64)"
        );
    }
}

} // End namespace PstreamDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The allToAll variants

template<class T, unsigned N>
void Foam::UPstream::allToAll
(
    const UList<FixedList<T, N>>& sendData,
    UList<FixedList<T, N>>& recvData,
    const int communicator
)
{
    PstreamDetail::allToAll_fixedWidth<FixedList<T, N>, T, int(N)>
    (
        sendData, recvData, communicator,
        [=](const UList<T>& sends, UList<T>& recvs, int count)
        {
            UPstream::allToAll_n(sends, recvs, count, communicator);
        }
    );
}


template<class T>
void Foam::UPstream::allToAll
(
    const UList<Pair<T>>& sendData,
    UList<Pair<T>>& recvData,
    const int communicator
)
{
    PstreamDetail::allToAll_fixedWidth<Pair<T>, T, 2>
    (
        sendData, recvData, communicator,
        [=](const UList<T>& sends, UList<T>& recvs, int count)
        {
            UPstream::allToAll_n(sends, recvs, count, communicator);
        }
    );
}


template<class T, unsigned N>
void Foam::UPstream::allToAllConsensus
(
    const UList<FixedList<T, N>>& sendData,
    UList<FixedList<T, N>>& recvData,
    const int tag,
    const int communicator
)
{
    PstreamDetail::allToAll_fixedWidth<FixedList<T, N>, T, int(N)>
    (
        sendData, recvData, communicator,
        [=](const UList<T>& sends, UList<T>& recvs, int count)
        {
            UPstream::allToAllConsensus_n
            (
                sends, recvs, count, tag, communicator
            );
        }
    );
}


template<class T>
void Foam::UPstream::allToAllConsensus
(
    const UList<Pair<T>>& sendData,
    UList<Pair<T>>& recvData,
    const int tag,
    const int communicator
)
{
    PstreamDetail::allToAll_fixedWidth<Pair<T>, T, 2>
    (
        sendData, recvData, communicator,
        [=](const UList<T>& sends, UList<T>& recvs, int count)
        {
            UPstream::allToAllConsensus_n
            (
                sends, recvs, count, tag, communicator
            );
        }
    );
}


// ************************************************************************* //
