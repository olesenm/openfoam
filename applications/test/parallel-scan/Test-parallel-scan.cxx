/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

Application
    Test-parallel-scan

Description
    Simple tests for MPI_Scan/MPI_Exscan

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "vector.H"
#include "Pstream.H"
#include "globalIndex.H"
#include "globalOffset.H"

using namespace Foam;


template<class Type>
void testGlobalIndex(const Type& localSize)
{
    const auto comm = UPstream::worldComm;

    // With find_offset
    {
        Type start = UPstream::find_offset(localSize, comm);

        Type total(start + localSize);

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, UPstream::nProcs(comm)-1);
        }

        Pout<< "size: " << localSize
            << " offset: " << start
            << " total: " << total
            << " (UPstream::find_offset)" << endl;
    }

    // With calcRange
    {
        IntRange<Type> range = globalIndex::calcRange(localSize, comm);

        Type total = range.end_value();

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, UPstream::nProcs(comm)-1);
        }

        Pout<< "range: " << range
            << " total: " << total
            << " (globalIndex calcRange)" << endl;
    }

    // With calculate OffsetRange
    {
        // Bad: possible narrowing mismatch
        // OffsetRange<Type> range =
        //     globalOffset::calculate(localSize, comm);

        // Possible implicit narrowing
        // auto range = globalOffset::calculate(localSize, comm);

        // OK: no possible mismatch
        auto range = GlobalOffset<Type>::calculate(localSize, comm);

        GlobalOffset<Type> other(10);

        Pout<< "range: " << range
            << " (globalOffset calculate)" << endl;

        if (!range.empty())
        {
            other = range;
            range.clear();
        }
    }
}


template<class Type, unsigned N>
void testOffsetTotals(const FixedList<Type, N>& localSizes)
{
    const auto comm = UPstream::worldComm;

    FixedList<Type, N> starts(localSizes);
    FixedList<Type, N> totals(localSizes);

    // The starting offsets

    // This actually compiles and does the right thing!
    // auto starts = UPstream::mpiExscan_sum(localSizes, comm);

    UPstream::mpiExscan_sum(starts.data(), starts.size(), comm);

    if (UPstream::is_parallel(comm))
    {
        // Broadcast totals from rank=N-1
        const int root = (UPstream::nProcs(comm)-1);

        // Update the totals.
        // - can actually do this on all ranks, but they
        // will be overwritten by the broadcast anyhow.
        if (root == UPstream::myProcNo(comm))
        {
            for (unsigned i = 0; i < N; ++i)
            {
                totals[i] += starts[i];
            }
        }

        UPstream::broadcast(totals, comm, root);
    }

    Pout<< "sizes: " << localSizes
        << " offsets: " << starts
        << " totals: " << totals
        << endl;

    // As tuples
    FixedList<std::pair<Type, Type>, N> tuples;
    for (unsigned i = 0; i < N; ++i)
    {
        tuples[i].first = starts[i];
        tuples[i].second = totals[i];
    }

    Pout<< "tuples: " << flatOutput(tuples) << endl;
}


template<bool exclusive, class Type>
void testScan(const Type& localValue)
{
    const auto comm = UPstream::worldComm;

    // sum (single value)
    {
        Type result;

        if constexpr (exclusive)
        {
            result = UPstream::mpiExscan_sum(localValue, comm);
        }
        else
        {
            result = UPstream::mpiScan_sum(localValue, comm);
        }

        Type total(result);
        if constexpr (exclusive)
        {
            total += localValue;
        }

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, (UPstream::nProcs(comm)-1));
        }

        if constexpr (exclusive)
        {
            Pout<< "input: " << localValue
                << " beg-sum: " << result
                << " total: " << total << endl;
        }
        else
        {
            Pout<< "input: " << localValue
                << " end-sum: " << result
                << " total: " << total << endl;
        }
    }

    // min (single value)
    {
        Type result;

        if constexpr (exclusive)
        {
            result = UPstream::mpiExscan_min(localValue, comm);
        }
        else
        {
            result = UPstream::mpiScan_min(localValue, comm);
        }

        if constexpr (exclusive)
        {
            Pout<< "input: " << localValue
                << " (exclusive)min : " << result
                << endl;
        }
        else
        {
            Pout<< "input: " << localValue
                << " (inclusive)min : " << result
                << endl;
        }
    }

    // max (this time with templated version)
    {
        constexpr auto opCode = UPstream::opCodes::op_max;

        Type result = UPstream::mpiScan<opCode, exclusive>(localValue, comm);

        if constexpr (exclusive)
        {
            Pout<< "input: " << localValue
                << " (exclusive)max : " << result
                << endl;
        }
        else
        {
            Pout<< "input: " << localValue
                << " (inclusive)max : " << result
                << endl;
        }
    }


    // Same as calculate OffsetRange
    if constexpr (std::is_integral_v<Type>)
    {
        const auto& size = localValue;

        // Exscan only!
        Type start = UPstream::find_offset(size, comm);

        Type total(start + size);

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, (UPstream::nProcs(comm)-1));
        }

        OffsetRange<Type> range(start, size, total);

        Pout<< "input: " << size
            << " range: " << range << endl;
    }
}


template<class Type>
void testScan(const Type& localValue, bool exclusive)
{
    if (exclusive)
    {
        testScan<true>(localValue);
    }
    else
    {
        testScan<false>(localValue);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();

    argList::addBoolOption("exclusive", "Perform exclusive scan");
    argList::addOption("factor", "Multiplier for rank values");

    #include "setRootCase.H"

    if (!UPstream::parRun())
    {
        Info<< "parallel only!" << nl;
        return 1;
    }

    const auto myProci = UPstream::myProcNo();

    const bool exclusive = args.found("exclusive");
    const scalar scaleFactor = args.getOrDefault<scalar>("factor", 1);


    // Note: regular call will cast the parameter so that the following
    // is not flagged as a problem
    {
        float value(3.14159);
        auto range0 = globalOffset::calculate(value);

        Info<< "offset of " << value << " : " << range0 << nl;

        // This one will refuse to compile...
        // auto range1 = globalOffset::calculate<float>(value);
    }

    Info<< nl;
    {
        typedef int32_t Type;

        Type localValue = (myProci+1)*scaleFactor;
        List<Type> allSizes = UPstream::listGatherValues(localValue);

        Info<< "all sizes: " << flatOutput(allSizes) << nl;

        testScan(localValue, exclusive);

        testGlobalIndex(localValue);
    }

    Info<< nl;
    {
        typedef int64_t Type;

        Type localValue = (myProci+1)*scaleFactor;
        List<Type> allSizes = UPstream::listGatherValues(localValue);

        Info<< "all sizes: " << flatOutput(allSizes) << nl;

        testScan(localValue, exclusive);

        testGlobalIndex(localValue);
    }

    Info<< nl;
    {
        typedef vector Type;

        Type localValue = vector::uniform(myProci+1);
        localValue *= scaleFactor;

        testScan(localValue, exclusive);
    }

    Info<< nl;
    {
        typedef int64_t Type;

        Type localValue = (myProci+1)*10;
        List<Type> allSizes = UPstream::listGatherValues(localValue);

        Info<< "all sizes: " << flatOutput(allSizes) << nl;
        Info<< "recv-sizes: " << globalIndex::calcRecvSizes(localValue) << nl;
    }

    Info<< nl;
    {
        typedef int64_t Type;
        typedef FixedList<Type, 4> list_type;

        FixedList<Type, 4> localSizes;
        std::iota(localSizes.begin(), localSizes.end(), (myProci+1)*10);

        typedef int64_t Type;

        List<list_type> allSizes = UPstream::listGatherValues(localSizes);
        Info<< "all  sizes: " << flatOutput(allSizes) << nl;

        testOffsetTotals(localSizes);
    }

    Info<< nl;
    {
        typedef double Type;

        Type localValue = (myProci+1)*scaleFactor;
        List<Type> allSizes = UPstream::listGatherValues(localValue);

        Info<< "all sizes: " << flatOutput(allSizes) << nl;

        testScan(localValue, exclusive);

        // This should refuse to compile
        // testGlobalIndex(localValue);
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
