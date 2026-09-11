/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Inter-processor communication reduction functions.

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamReduceOps_H
#define Foam_PstreamReduceOps_H

#include "Pstream.H"
#include "FixedList.H"
#include "MinMax.H"
#include "ops.H"
#include "VectorSpaceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Single value: mpiAllReduce or manual

//- Reduce inplace (cf. MPI Allreduce)
template<class T, class BinaryOp>
void reduce
(
    T& value,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiAllReduce
        (
            &value,
            1,
            UPstream_opType<BinaryOp>::opcode_id,
            communicator
        );
    }
    else
    {
        if (UPstream::warnComm >= 0 && communicator != UPstream::warnComm)
        {
            Perr<< "** reducing:" << value << " comm:" << communicator << endl;
            error::printStack(Perr);
        }
        Pstream::gather(value, bop, tag, communicator);
        Pstream::broadcast(value, communicator);
    }
}


// Multiple values: mpiAllReduce only!

//- Reduce inplace (cf. MPI Allreduce)
//- multiple values (same size on all ranks!)
template<class T, class BinaryOp>
void reduce
(
    T values[],
    const int count,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int communicator = UPstream::worldComm
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiAllReduce
        (
            values,
            count,
            UPstream_opType<BinaryOp>::opcode_id,
            communicator
        );
    }
    else
    {
        // static_assert(false, "unsupported data type");
        static_assert
        (
            stdFoam::dependent_false_v<T>,
            "only for specific, contiguous, known data types"
        );
    }
}


//- Reduce multiple values
//  Multiple values: mpiAllReduce only!
template<class T, unsigned N, class BinaryOp>
inline void reduce
(
    FixedList<T, N>& values,
    BinaryOp bop,
    const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    reduce(values.data(), int(values.size()), bop, tag, communicator);
}


//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- multiple values. Sets request.
template<class T, class BinaryOp>
void reduce
(
    T values[],
    int count,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int communicator,
    UPstream::Request& req
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiAllReduce
        (
            values,
            count,
            UPstream_opType<BinaryOp>::opcode_id,
            communicator,
            req
        );
    }
    else
    {
        // static_assert(false, "unsupported data type");
        static_assert
        (
            stdFoam::dependent_false_v<T>,
            "only for specific, contiguous, known data types"
        );
    }
}


//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- single value. Sets request.
template<class T, class BinaryOp>
void reduce
(
    T& value,
    BinaryOp bop,
    const int tag,
    const int communicator,
    UPstream::Request& req
)
{
    // single value
    reduce(&value, 1, bop, tag, communicator, req);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Special reductions for bool

//- Logical (and) inplace reduction. Uses UPstream::reduceAnd
inline void reduce
(
    bool& value,
    Foam::andOp<bool>,
    [[maybe_unused]] const int tag = UPstream::msgType(),  /*!< (ignored) */
    const int communicator = UPstream::worldComm
)
{
    UPstream::reduceAnd(value, communicator);
}


//- Logical (or) inplace reduction. Uses UPstream::reduceOr
inline void reduce
(
    bool& value,
    Foam::orOp<bool>,
    [[maybe_unused]] const int tag = UPstream::msgType(),  /*!< (ignored) */
    const int communicator = UPstream::worldComm
)
{
    UPstream::reduceOr(value, communicator);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Special reductions for MinMax

//- Reduce inplace (Allreduce) for MinMax with sum/plus semantics, which
//- corresponds to op_min on the min() value and op_max on the max() value
//
//  Note: similar to boundBox::reduce()
template<class T>
void reduce
(
    MinMax<T>& value,
    Foam::sumOp<MinMax<T>>, //<- or std::plus<void>
    [[maybe_unused]] const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    constexpr auto nCmpts = int(pTraits_nComponents<T>::value);

    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (nCmpts > 0 && UPstream_basic_dataType<T>::value)
    {
        // MinMax (directly/indirectly) uses basic dataType.
        // sumOp/plusOp combining is just a min-reduce and max-reduce

        #ifdef FULLDEBUG
        if (UPstream::warnComm >= 0 && communicator != UPstream::warnComm)
        {
            Perr<< "** min/max mpi-reduce:" << value << " comm:"
                << communicator << endl;
            error::printStack(Perr, 4);
        }
        #endif

        UPstream::mpiAllReduce
        (
            vectorspace_data(value.min()),
            nCmpts,
            UPstream::opCodes::op_min,
            communicator
        );
        UPstream::mpiAllReduce
        (
            vectorspace_data(value.max()),
            nCmpts,
            UPstream::opCodes::op_max,
            communicator
        );
    }
    else
    {
        // Unusual case: MinMax with some non-standard type of content
        #ifdef FULLDEBUG
        if (UPstream::warnComm >= 0 && communicator != UPstream::warnComm)
        {
            Perr<< "** min/max reduce:" << value << " comm:"
                << communicator << endl;
            error::printStack(Perr, 4);
        }
        #endif
        Pstream::gather(value, sumOp<MinMax<T>>(), tag, communicator);
        Pstream::broadcast(value, communicator);
    }
}


//- Reduce inplace (Allreduce) for MinMax with sum/plus semantics
//  Note: plusOp and sumOp are synonymous
template<class T>
inline void reduce
(
    MinMax<T>& value,
    Foam::plusOp<MinMax<T>>,
    const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    Foam::reduce(value, Foam::sumOp<MinMax<T>>(), tag, communicator);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Reduce inplace (cf. MPI Allreduce)
//- the sum of value and counter (eg, for averaging)
//  \note the summed counter value may lose some precision if components
//  of the main value are narrower. This is usually not the case.
template<class T, class CountType = int>
void sumReduce
(
    //! The value to sum
    T& value,
    //! The integer count or floating-point weight etc
    CountType& count,
    [[maybe_unused]] const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    static_assert(std::is_arithmetic_v<CountType>, "Counter not arithmetic");

    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (std::is_floating_point_v<T>)
    {
        // Bundle value and count together
        T work[2];
        work[0] = value;
        work[1] = static_cast<T>(count);

        UPstream::mpiAllReduce
        (
            work,
            2,
            UPstream::opCodes::op_sum,
            communicator
        );

        // Unbundle
        value = work[0];
        count = static_cast<CountType>(work[1]);
    }
    else if constexpr
    (
        is_vectorspace_v<T>
     && std::is_floating_point_v<typename pTraits_cmptType<T>::type>
    )
    {
        constexpr auto nCmpts = pTraits_nComponents<T>::value;
        using cmpt = typename pTraits_cmptType<T>::type;

        // Bundle all components and count together
        cmpt work[nCmpts+1];
        VectorSpaceOps<nCmpts>::copy_n(value.begin(), work);
        work[nCmpts] = static_cast<cmpt>(count);

        UPstream::mpiAllReduce
        (
            work,
            nCmpts+1,
            UPstream::opCodes::op_sum,
            communicator
        );

        // Unbundle
        VectorSpaceOps<nCmpts>::copy_n(work, value.begin());
        count = static_cast<CountType>(work[nCmpts]);
    }
    else
    {
        Foam::reduce(value, sumOp<T>(), tag, communicator);
        Foam::reduce(count, sumOp<CountType>(), tag, communicator);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convenience wrappers - defined after all specialisations are known

//- Perform reduction on a copy, using specified binary operation
//  \return the resulting value
template<class T, class BinaryOp>
T returnReduce
(
    const T& value,
    BinaryOp bop,
    const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    T work(value);
    Foam::reduce(work, bop, tag, communicator);
    return work;
}


//- Perform logical (and) MPI Allreduce on a copy. Uses UPstream::reduceAnd
//  \return the resulting value
inline bool returnReduceAnd
(
    const bool value,
    const int communicator = UPstream::worldComm
)
{
    bool work(value);
    UPstream::reduceAnd(work, communicator);
    return work;
}


//- Perform logical (or) MPI Allreduce on a copy. Uses UPstream::reduceOr
//  \return the resulting value
inline bool returnReduceOr
(
    const bool value,
    const int communicator = UPstream::worldComm
)
{
    bool work(value);
    UPstream::reduceOr(work, communicator);
    return work;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
