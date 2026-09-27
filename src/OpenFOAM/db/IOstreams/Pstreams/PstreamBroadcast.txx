/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// The max number of elements that can be sent (contiguous or serialized)
// without exceeding MPI <int> buffer limits.
// Never larger than the list itself, so the result can be used directly
// as a chunk length (no chunking needed if result == values.size()).
template<class T>
label listMaxChunkSize(const UList<T>& values)
{
    // The MPI (int) message-size limit with a small safety margin
    constexpr auto maxBufferSize = label(INT_MAX - 63);

    label maxElements = maxBufferSize;

    if constexpr (is_contiguous_v<T>)
    {
        // Using element or component type (or byte-wise) size
        maxElements = (maxBufferSize / label(UPstream_dataType<T>::size(1)));
    }
    else if (values.size() >= 100000)
    {
        // Larger list, may need chunking for serialization.
        // - estimate the size by serializing a few items

        DynamicList<char> scratch(1024);
        UOPstream os(scratch);

        constexpr int averaging = 8;

        for (int i = 0; i < averaging; ++i)
        {
            os << values[i];
        }
        os.alignp(8);

        maxElements = (maxBufferSize / label(os.tellp()/averaging));
    }
    // else: smaller list, very unlikely to need chunking

    return Foam::min(maxElements, values.size());
}


// Broadcast (sync) UList content without exceeding MPI <int> buffer limits.
// Lists are same length on all ranks!
//
// Respects UPstream::listChunkLength()
template<class T>
void UList_chunked_broadcast(UList<T>& values, int communicator)
{
    if (values.empty())
    {
        return;
    }
    else if (values.size() == 1)
    {
        // Single value - optimized version
        Pstream::broadcast(values[0], communicator);
        return;
    }

    // Multiple values - same length on all ranks
    const label listLen = values.size();

    // Chunk length. User defined, or MPI max limits.
    // (will simply be the listLen if no chunking is required)
    const label chunkLen =
    (
        (UPstream::listChunkLength() >= UPstream::listChunkLength_min())
      ? UPstream::listChunkLength()
      : PstreamDetail::listMaxChunkSize(values)
    );

    // Operate on chunks, single pass when no chunking is needed
    for (label start = 0; start < listLen; start += chunkLen)
    {
        const auto len = Foam::min(chunkLen, listLen - start);

        auto slice = values.slice(start, len);

        if constexpr (is_contiguous_v<T>)
        {
            // Broadcast contiguous content
            UPstream::broadcast(slice.data(), slice.size(), communicator);
        }
        else
        {
            // Non-contiguous content - serialize it
            // Same as Pstream::broadcast(slice, communicator);

            if (UPstream::master(communicator))
            {
                OPBstream::send(slice, communicator);
            }
            else
            {
                IPBstream::recv(slice, communicator);
            }
        }
    }
}

} // End namespace PstreamDetail
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::Pstream::broadcast
(
    Type& value,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (std::is_same_v<bitSet, Type>)
    {
        // Handled properly by bitSet itself (is internal contiguous)
        value.broadcast(communicator);
    }
    else if constexpr (is_contiguous_v<Type>)
    {
        // A single-value and contiguous
        UPstream::broadcast(&value, 1, communicator);
    }
    else
    {
        // Non-contiguous content - serialize it
        if (UPstream::master(communicator))
        {
            OPBstream::send(value, communicator);
        }
        else
        {
            IPBstream::recv(value, communicator);
        }
    }
}


template<class Type, unsigned N>
void Foam::Pstream::broadcast
(
    FixedList<Type, N>& list,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (is_contiguous_v<Type>)
    {
        UPstream::broadcast(list.data(), list.size(), communicator);
    }
    else
    {
        // Non-contiguous content - serialize it
        if (UPstream::master(communicator))
        {
            OPBstream::send(list, communicator);
        }
        else
        {
            IPBstream::recv(list, communicator);
        }
    }
}


template<class Type, class... Args>
void Foam::Pstream::broadcasts
(
    const int communicator,
    Type& value,
    Args&&... values
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (!sizeof...(values) && is_contiguous_v<Type>)
    {
        // A single-value and contiguous
        UPstream::broadcast(&value, 1, communicator);
    }
    else
    {
        // Non-contiguous data, or multiple data - needs serialization

        if (UPstream::master(communicator))
        {
            OPBstream::sends
            (
                communicator,
                value,
                std::forward<Args>(values)...
            );
        }
        else
        {
            IPBstream::recvs
            (
                communicator,
                value,
                std::forward<Args>(values)...
            );
        }
    }
}


template<class ListType>
void Foam::Pstream::broadcastList
(
    ListType& list,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (std::is_same_v<bitSet, ListType>)
    {
        // Handled properly by bitSet itself (is internal contiguous)
        list.broadcast(communicator);
    }
    else if constexpr (is_contiguous_v<typename ListType::value_type>)
    {
        // List data are contiguous
        // 1. broadcast the size
        // 2. resize for receiver list
        // 3. broadcast contiguous contents

        int64_t len(list.size());

        UPstream::mpi_broadcast
        (
            &len,
            1,
            UPstream::dataTypes::type_int64,
            communicator
        );

        if (UPstream::is_subrank(communicator))
        {
            list.resize_nocopy(len);
        }

        if (len)
        {
            // Only broadcast non-empty content
            UPstream::broadcast(list.data(), list.size(), communicator);
        }
    }
    else
    {
        // List data are non-contiguous - serialize/de-serialize

        if (UPstream::master(communicator))
        {
            if (list.empty())
            {
                // Do not serialize if empty.
                // Just broadcast zero-size in a form that IPBstream can expect
                OPBstream::send(Foam::zero{}, communicator);
            }
            else
            {
                OPBstream::send(list, communicator);
            }
        }
        else
        {
            IPBstream is(communicator);
            if (is.remaining() > 0)  // Received a non-empty buffer
            {
                is >> list;
            }
            else
            {
                list.clear();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convenience wrappers - defined after all specialisations are known

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a broadcasted value (uses a copy internally)
template<class Type>
Type returnBroadcast
(
    const Type& value,
    const int communicator = UPstream::worldComm
)
{
    Type work(value);
    Pstream::broadcast(work, communicator);
    return work;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
