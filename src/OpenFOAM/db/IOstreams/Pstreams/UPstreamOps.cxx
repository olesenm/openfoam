/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Reduction * * * * * * * * * * * * * * * * //

int Foam::UPstream::find_first(bool condition, int communicator)
{
    if (UPstream::is_parallel(communicator))
    {
        // Ensure we use exact and known data types
        typedef int32_t IntType;
        constexpr auto dataType = UPstream::dataTypes::type_int32;
        constexpr auto opCodeId = UPstream::opCodes::op_min;

        auto proci = static_cast<IntType>(UPstream::myProcNo(communicator));
        auto nproc = static_cast<IntType>(UPstream::nProcs(communicator));

        if (!condition)
        {
            // The out-of-range value (min reduction)
            proci = nproc;
        }

        UPstream::mpi_allreduce(&proci, 1, dataType, opCodeId, communicator);

        // Found: the first rank (min). Not found: -1
        return (proci < nproc ? proci : -1);
    }
    else
    {
        // Not parallel (serial), so we are rank=0
        // return 0 (found) or -1 (not found)
        return (condition ? 0 : -1);
    }
}


int Foam::UPstream::find_last(bool condition, int communicator)
{
    if (UPstream::is_parallel(communicator))
    {
        // Ensure we use exact and known data types
        typedef int32_t IntType;
        constexpr auto dataType = UPstream::dataTypes::type_int32;
        constexpr auto opCodeId = UPstream::opCodes::op_max;

        auto proci = static_cast<IntType>(UPstream::myProcNo(communicator));

        if (!condition)
        {
            // The out-of-range value (max reduction)
            proci = -1;
        }

        UPstream::mpi_allreduce(&proci, 1, dataType, opCodeId, communicator);

        // Found: the last rank (max). Not found: -1
        return (proci >= 0 ? proci : -1);
    }
    else
    {
        // Not parallel (serial), so we are rank=0
        // return 0 (found) or -1 (not found)
        return (condition ? 0 : -1);
    }
}


// Implementation detail for int64 minLocation/maxLocation on LLP64
// (eg, mingw/Windows, where 'long' is only 32-bit, so MPI_LONG_INT can't
// carry a full int64). Splits each value into (hi32, lo32) and does two
// portable (int,int) MINLOC/MAXLOC passes instead:
//
// - Pass 1 reduces (hi32, lo32-as-tiebreak).
//   NB: minloc/maxloc *always* pick the smaller of the two tie-break
//   fields (even for maxloc!), so lo32 is complemented for maxloc
//   to flip that sense. This makes the combined compare an exact
//   lexicographic 64-bit compare, giving the correct global value directly.
//
// - Pass 2 recovers the index.
//   Candidates whose local hi32/lo32 don't match the pass-1 winner are
//   forced to a losing key; the rest (which must include the true winner)
//   are compared on a plain match/no-match flag, so the usual
//   smaller-index tie-break then yields the smallest original index among
//   any exact ties.

template<Foam::UPstream::opCodes OpCode>
void Foam::UPstream::allReduceLocation_int64_split
(
    Tuple2<int64_t, int>* tuples,
    int64_t count,
    int communicator
)
{
    static_assert
    (
        (
            (OpCode == UPstream::opCodes::op_minloc)
         || (OpCode == UPstream::opCodes::op_maxloc)
        ),
        "Only supports (minloc, maxloc) operations"
    );

    if (!count || !tuples || !UPstream::is_parallel(communicator))
    {
        return;
    }

    constexpr bool isMaxOp = (OpCode == UPstream::opCodes::op_maxloc);

    // For the new comparison key.
    // match = 0, to maintain for lexicographical order the non-match is
    // (-1, +1) for max, min operations, respectively.

    constexpr int nonMatch = (isMaxOp ? -1 : 1);

    // The MPI struct for both (portable, 32-bit) reduction passes
    struct mpiTuple { int first; int second; };

    const auto reduceChunked = [communicator](mpiTuple* data, int64_t count)
    {
        // The MPI (int) message-size limit, but since we are chunking
        // anyhow can take a smaller size too. Using half of that range.
        constexpr auto maxChunk = int64_t(INT_MAX/2+1);

        for (int64_t start = 0; start < count; start += maxChunk)
        {
            UPstream::mpi_allreduce
            (
                (data + start),
                int(std::min(maxChunk, count - start)),
                UPstream::dataTypes::type_int_int,  //< (int int)
                OpCode,
                communicator
            );
        }
    };

    // Work array of mpiTuple
    List<mpiTuple> work(count);

    // Note: the lower 32 bits are bias-flipped (the xor 0x80000000u)
    // so that unsigned ordering matches signed int comparison
    // (needed since MPI's index/tie-break field compares as a signed int)

    for (int64_t i = 0; i < count; ++i)
    {
        // From the original value into (high, biased-low)
        const int64_t origValue = tuples[i].first();
        const int32_t origUpper = int32_t(uint64_t(origValue) >> 32);
        const int32_t origLowerB = int32_t(uint32_t(origValue) ^ 0x80000000u);

        work[i].first  = origUpper;

        // MPI always breaks ties on the smaller second field, even for
        // MAXLOC, so invert the bits for the max case to flip that sense
        work[i].second = isMaxOp ? ~origLowerB : origLowerB;
    }

    reduceChunked(work.data(), count);

    for (int64_t i = 0; i < count; ++i)
    {
        // From the original value into (high, biased-low)
        const int64_t origValue = tuples[i].first();
        const int32_t origUpper = int32_t(uint64_t(origValue) >> 32);
        const int32_t origLowerB = int32_t(uint32_t(origValue) ^ 0x80000000u);

        // From the encoded, reduced value:
        const int32_t encUpper  = work[i].first;
        const int32_t encLowerB = isMaxOp ? ~work[i].second : work[i].second;

        // Reconstruct the true global extreme value and store it now
        tuples[i].first() =
        (
            (int64_t(encUpper) << 32)
          | (uint32_t(encLowerB) ^ 0x80000000u)
        );

        // The new comparison key. match = 0, non-match is (-1, +1)
        // for max, min operations respectively

        work[i].first =
        (
            (origUpper == encUpper && origLowerB == encLowerB)
          ? 0 //< match
          : nonMatch
        );

        work[i].second = tuples[i].second();  // original payload index
    }

    // repeat
    reduceChunked(work.data(), count);

    for (int64_t i = 0; i < count; ++i)
    {
        tuples[i].second() = work[i].second;
    }
}


// Implementation detail for minLocation/maxLocation
template<class Type, Foam::UPstream::opCodes OpCode>
void Foam::UPstream::allReduceLocation_impl
(
    Tuple2<Type, int>* tuples,
    int64_t count,
    int communicator
)
{
    if (!count || !tuples || !UPstream::is_parallel(communicator))
    {
        return;
    }

    if constexpr
    (
        std::is_same_v<Type, int64_t>
     && (sizeof(long) < sizeof(int64_t))
    )
    {
        // LLP64 (eg, mingw/Windows): 'long' is only 32-bit, so MPI_LONG_INT
        // cannot carry a full int64 value without truncation. Use the
        // portable two-pass 32-bit split instead (see above).
        allReduceLocation_int64_split<OpCode>(tuples, count, communicator);
        return;
    }
    else
    {

    using basic_type = typename UPstream_location_dataType<Type>::type;
    constexpr auto datatype_id = UPstream_location_dataType<Type>::datatype_id;

    // The C++ container
    using  cppTuple = Tuple2<Type, int>;

    // The struct for MPI
    struct mpiTuple { basic_type first; int second; };

    static_assert
    (
        (
            (datatype_id == UPstream::dataTypes::type_int_int)
         || (datatype_id == UPstream::dataTypes::type_long_int)
         || (datatype_id == UPstream::dataTypes::type_float_int)
         || (datatype_id == UPstream::dataTypes::type_double_int)
        ),
        "Only supports (int_int, long_int, float_int, double_int) tuples"
    );

    static_assert
    (
        (
            (OpCode == UPstream::opCodes::op_minloc)
         || (OpCode == UPstream::opCodes::op_maxloc)
        ),
        "Only supports (minloc, maxloc) operations"
    );

    static_assert
    (
        (sizeof(cppTuple) >= sizeof(mpiTuple)),
        "Tuple sizeof is not >= sizeof MPI struct. Transcribing will fail!"
    );

    if (count == 1)
    {
        // Fast(er) path
        auto& data = *tuples;

        // Copy into the C-struct for exact alignment etc.
        mpiTuple blob { basic_type(data.first()), data.second() };

        UPstream::mpi_allreduce(&blob, 1, datatype_id, OpCode, communicator);

        // Copy back from C-struct to original locations
        data.first() = Type(blob.first);
        data.second() = blob.second;
        return;
    }

    // Regular path
    auto* blobs = reinterpret_cast<mpiTuple*>(tuples);

    // Transcribe to ensure data layout is exact, without padding
    for (int64_t i = 0; i < count; ++i)
    {
        mpiTuple t { basic_type(tuples[i].first()), tuples[i].second() };
        blobs[i] = t;
    }

    // The MPI (int) message-size limit, but since we are chunking
    // anyhow can take a smaller size too. Using half of that range.
    constexpr auto maxChunk = int64_t(INT_MAX/2+1);

    // Operate on chunks, single pass when no chunking is needed
    for (int64_t start = 0; start < count; start += maxChunk)
    {
        UPstream::mpi_allreduce
        (
            (blobs + start),
            int(std::min(maxChunk, count - start)),
            datatype_id,
            OpCode,
            communicator
        );
    }

    // Copy back from C-struct to original locations. Reverse order!
    for (int64_t i = count-1; i >= 0; --i)
    {
        cppTuple t( Type(blobs[i].first), blobs[i].second );
        tuples[i] = t;
    }

    } // end else (native fast path)
}


template<class Type, Foam::UPstream::opCodes OpCode>
int Foam::UPstream::findLocation_impl(Type value, int communicator)
{
    constexpr auto datatype_id = UPstream_location_dataType<Type>::datatype_id;

    static_assert
    (
        (
            (datatype_id == UPstream::dataTypes::type_int_int)
         || (datatype_id == UPstream::dataTypes::type_long_int)
         || (datatype_id == UPstream::dataTypes::type_float_int)
         || (datatype_id == UPstream::dataTypes::type_double_int)
        ),
        "Only supports (int_int, long_int, float_int, double_int) tuples"
    );

    static_assert
    (
        (
            (OpCode == UPstream::opCodes::op_minloc)
         || (OpCode == UPstream::opCodes::op_maxloc)
        ),
        "Only supports (minloc, maxloc) operations"
    );

    if (!UPstream::is_parallel(communicator))
    {
        // For non-parallel (serial) we are rank=0
        return 0;
    }

    int proci = UPstream::myProcNo(communicator);

    if constexpr
    (
        std::is_same_v<Type, int64_t>
     && (sizeof(long) < sizeof(int64_t))
    )
    {
        // See note about LLP64 (mingw/Windows) elsewhere in this file
        Tuple2<int64_t, int> tup2(value, proci);

        allReduceLocation_int64_split<OpCode>(&tup2, 1, communicator);

        proci = tup2.second();
    }
    else
    {
        using basic_type = typename UPstream_location_dataType<Type>::type;

        // The struct for MPI
        struct mpiTuple { basic_type first; int second; };

        mpiTuple blob { basic_type(value), proci };

        UPstream::mpi_allreduce(&blob, 1, datatype_id, OpCode, communicator);

        proci = blob.second;
    }

    return proci;
}


#undef  Pstream_defineMethods
#define Pstream_defineMethods(Type)                                           \
                                                                              \
int Foam::UPstream::find_minLocation(Type value, int communicator)            \
{                                                                             \
    return                                                                    \
    findLocation_impl<Type,UPstream::opCodes::op_minloc>(value, communicator);\
}                                                                             \
                                                                              \
int Foam::UPstream::find_maxLocation(Type value, int communicator)            \
{                                                                             \
    return                                                                    \
    findLocation_impl<Type,UPstream::opCodes::op_maxloc>(value, communicator);\
}

// NB: only four types are ever possible (int, long, float, double)
Pstream_defineMethods(double);

#undef Pstream_defineMethods


#undef  Pstream_defineMethods
#define Pstream_defineMethods(Type)                                           \
                                                                              \
void Foam::UPstream::minLocation                                              \
(                                                                             \
    Tuple2<Type, int>* data_tuples, int64_t count, int communicator           \
)                                                                             \
{                                                                             \
    allReduceLocation_impl<Type, UPstream::opCodes::op_minloc>                \
    ( data_tuples, count, communicator );                                     \
}                                                                             \
                                                                              \
void Foam::UPstream::maxLocation                                              \
(                                                                             \
    Tuple2<Type, int>* data_tuples, int64_t count, int communicator           \
)                                                                             \
{                                                                             \
    allReduceLocation_impl<Type, UPstream::opCodes::op_maxloc>                \
    ( data_tuples, count, communicator );                                     \
}                                                                             \
                                                                              \
void Foam::UPstream::minLocation(UList<Tuple2<Type, int>>& tuples, int comm)  \
{                                                                             \
    UPstream::minLocation(tuples.data(), tuples.size(), comm);                \
}                                                                             \
                                                                              \
void Foam::UPstream::maxLocation(UList<Tuple2<Type, int>>& tuples, int comm)  \
{                                                                             \
    UPstream::maxLocation(tuples.data(), tuples.size(), comm);                \
}                                                                             \
                                                                              \
int Foam::UPstream::minLocation(Tuple2<Type, int>& data_tuple, int comm)      \
{                                                                             \
    UPstream::minLocation(&data_tuple, 1, comm);                              \
    return data_tuple.second();                                               \
}                                                                             \
                                                                              \
int Foam::UPstream::maxLocation(Tuple2<Type, int>& data_tuple, int comm)      \
{                                                                             \
    UPstream::maxLocation(&data_tuple, 1, comm);                              \
    return data_tuple.second();                                               \
}

// NB: only four types are ever possible (int, long, float, double)
Pstream_defineMethods(int32_t);
Pstream_defineMethods(int64_t);
Pstream_defineMethods(float);
Pstream_defineMethods(double);

#undef Pstream_defineMethods


// ************************************************************************* //
