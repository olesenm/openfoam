/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Test-parallel-minloc

Description
    Basic test for MPI min_loc, max_loc operations

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "Tuple2.H"
#include "Random.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"
#include "clockTime.H"
#include "memInfo.H"

namespace Foam
{

// An ad hoc tuple marked as contiguous:
template<class T> struct contigTuple
{
    T first; int second;

    bool operator==(const contigTuple& t) const
    {
        return (first == t.first && second == t.second);
    }

    bool operator!=(const contigTuple& t) const
    {
        return !(*this == t);
    }
};

// An ad hoc tuple marked as non-contiguous:
template<class T> struct noncontigTuple
{
    T first; int second;

    bool operator==(const noncontigTuple& t) const
    {
        return (first == t.first && second == t.second);
    }

    bool operator!=(const noncontigTuple& t) const
    {
        return !(*this == t);
    }

};

template<class T> struct is_contiguous<contigTuple<T>> : std::true_type {};
template<class T> struct is_contiguous<noncontigTuple<T>> : std::false_type {};


template<class T>
Istream& operator>>(Istream& is, contigTuple<T>& t)
{
    is.readBegin("tuple");
    is >> t.first >> t.second;
    is.readEnd("tuple");
    is.check(FUNCTION_NAME);
    return is;
}

template<class T>
Istream& operator>>(Istream& is, noncontigTuple<T>& t)
{
    is.readBegin("tuple");
    is >> t.first >> t.second;
    is.readEnd("tuple");
    is.check(FUNCTION_NAME);
    return is;
}

template<class T>
Ostream& operator<<(Ostream& os, const contigTuple<T>& t)
{
    os << '(' << t.first << ' ' << t.second << ')';
    return os;
}

template<class T>
Ostream& operator<<(Ostream& os, const noncontigTuple<T>& t)
{
    os << '(' << t.first << ' ' << t.second << ')';
    return os;
}


// Local wrapper for UPstream (to allow protected access)
class UPstreamTester : public UPstream
{
private:

    template<class Type, UPstream::opCodes OpCode>
    static void mpiAllReduce_location
    (
        Tuple2<Type, int> values[],
        int count,
        int communicator
    )
    {
        // Restricted to basic location types
        UPstream::mpi_allreduce
        (
            values,
            count,
            UPstream_location_dataType<Type>::datatype_id,
            OpCode,
            communicator
        );
    }

public:

    using UPstream::mpiAllReduce_min;
    using UPstream::mpiAllReduce_max;

    #undef  Pstream_declareMethods
    #define Pstream_declareMethods(OpName, OpCode, Type)                      \
                                                                              \
        static void mpiAllReduce_##OpName                                     \
        (                                                                     \
            Tuple2<Type, int>& value,                                         \
            const int communicator = UPstream::worldComm                      \
        )                                                                     \
        {                                                                     \
            mpiAllReduce_location<Type, OpCode>(&value, 1, communicator);     \
        }

        Pstream_declareMethods(min, UPstream::opCodes::op_min, int32_t);
        Pstream_declareMethods(max, UPstream::opCodes::op_max, int32_t);

        Pstream_declareMethods(min, UPstream::opCodes::op_min, int64_t);
        Pstream_declareMethods(max, UPstream::opCodes::op_max, int64_t);

        Pstream_declareMethods(min, UPstream::opCodes::op_min, float);
        Pstream_declareMethods(max, UPstream::opCodes::op_max, float);

        Pstream_declareMethods(min, UPstream::opCodes::op_min, double);
        Pstream_declareMethods(max, UPstream::opCodes::op_max, double);

        #undef Pstream_declareMethods
};

} // End namespace Foam

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();

    argList::addBoolOption("basic", "basic tests");
    argList::addBoolOption("unsupported", "test 'unsupported' operations");

    argList::addBoolOption("all", "same as -mpi, -contiguous, -non-contiguous");
    argList::addBoolOption("mpi", "array reduction with native mpi op");
    argList::addBoolOption
    (
        "contiguous",
        "listCombineReduce contiguous"
    );
    argList::addBoolOption
    (
        "non-contiguous",
        "listCombineReduce non-contiguous"
    );

    // For list tests
    argList::addOption("length", "num", "array length (default: 1000)");
    argList::addOption("repeat", "count", "repeat (default: 10)");

    #include "setRootCase.H"

    Info<< nl << "mpi vendor=" << int(UPstream::vendorCode()) << nl;

    const label repeat = args.getOrDefault<label>("repeat", 10);
    const label length = args.getOrDefault<label>("length", 1000);

    const bool testAll = args.found("all");
    const bool testBasic = args.found("basic");
    const bool testNative = args.found("mpi");
    const bool testContig = args.found("contiguous");
    const bool testNoncontig = args.found("non-contiguous");

    // basic/trivial tests
    if (testBasic)
    {
        using Type = scalar;
        using tupleType = Tuple2<Type, int>;

        Random rndGen(10 + UPstream::myProcNo());

        tupleType value
        (
            rndGen.position<Type>(0, 100),
            UPstream::myProcNo()
        );


        List<Type> values = Pstream::listGatherValues(value.first());
        labelList indices = Foam::identity(UPstream::nProcs());

        Info<< "pre values:  "; values.writeList(Info) << endl;
        Info<< "pre indices: "; indices.writeList(Info) << endl;

        {
            tupleType copy(value);
            int best = -1;

            best = UPstream::minLocation(copy, UPstream::worldComm);

            Info<< "minloc: " << copy << " best:" << best << endl;

            copy = value;

            best = UPstream::maxLocation(copy, UPstream::worldComm);
            Info<< "maxloc: " << copy << " best:" << best << endl;
        }

        Info<< "find_minLocation: " << value.first() << " => "
            << UPstream::find_minLocation(value.first(), UPstream::worldComm)
            << endl;
        Info<< "find_maxLocation: " << value.first() << " => "
            << UPstream::find_maxLocation(value.first(), UPstream::worldComm)
            << endl;
    }

    // min/max/sum with tuple
    // These will fail at runtime

    // An error occurred in MPI_Allreduce:
    // the reduction operation MPI_MIN is not defined on
    // the MPI_DOUBLE_INT datatype

    if (args.found("unsupported"))
    {
        using Type = scalar;
        using tupleType = Tuple2<Type, int>;

        Random rndGen(10 + UPstream::myProcNo());

        tupleType value
        (
            rndGen.position<Type>(0, 100),
            UPstream::myProcNo()
        );


        List<Type> values = Pstream::listGatherValues(value.first());
        labelList indices = Foam::identity(UPstream::nProcs());

        Info<< "pre values:  "; values.writeList(Info) << endl;
        Info<< "pre indices: "; indices.writeList(Info) << endl;

        {
            tupleType copy(value);
            UPstreamTester::mpiAllReduce_min(copy, UPstream::worldComm);
            Info<< "min: " << copy << nl;

            copy = value;
            UPstreamTester::mpiAllReduce_max(copy, UPstream::worldComm);
            Info<< "max: " << copy << nl;
        }
    }

    Info<< nl << endl;

    if (testAll || testNative)
    {
        using Type = scalar;
        using tupleType = Tuple2<Type, int>;

        Random rndGen(10 + UPstream::myProcNo());

        List<tupleType> list(length);

        for (auto& tup : list)
        {
            tup = tupleType
            (
                rndGen.position<Type>(0, 100),
                UPstream::myProcNo()
            );
        }

        Info<< "mpi reduce: " << flush;

        clockTime timing;
        for (label count = 0; count < repeat; ++count)
        {
            UPstream::minLocation
            (
                list.data(),
                list.size(),
                UPstream::worldComm
            );
        }

        double speed = timing.timeIncrement();

        int64_t memPeak = memInfo::get_hwm();
        Foam::reduce(memPeak, maxOp<int64_t>());

        Info<< speed << " sec, mem:" << memPeak << endl;
    }

    if (testAll || testContig)
    {
        using Type = scalar;
        using tupleType = contigTuple<Type>;
        Random rndGen(10 + UPstream::myProcNo());

        List<tupleType> list(length);

        for (auto& tup : list)
        {
            tup = { rndGen.position<Type>(0, 100), UPstream::myProcNo() };
        }

        // minloc op
        const auto bop = [](tupleType& x, const tupleType& y)
        {
            if (y.first < x.first)
            {
                x = y;
            }
            else if (y.first == x.first)
            {
                x.second = std::min(x.second, y.second);
            }
        };

        Info<< "contiguous reduce: " << flush;

        clockTime timing;
        for (label count = 0; count < repeat; ++count)
        {
            Pstream::listCombineReduce
            (
                list,
                bop,
                UPstream::msgType(),
                UPstream::worldComm
            );
        }

        double speed = timing.timeIncrement();

        int64_t memPeak = memInfo::get_hwm();
        Foam::reduce(memPeak, maxOp<int64_t>());

        Info<< speed << " sec, mem:" << memPeak << endl;
    }

    if (testAll || testNoncontig)
    {
        using Type = scalar;
        using tupleType = noncontigTuple<Type>;

        Random rndGen(10 + UPstream::myProcNo());

        List<tupleType> list(length);

        for (auto& tup : list)
        {
            tup = { rndGen.position<Type>(0, 100), UPstream::myProcNo() };
        }

        // minloc op
        const auto bop = [](tupleType& x, const tupleType& y)
        {
            if (y.first < x.first)
            {
                x = y;
            }
            else if (y.first == x.first)
            {
                x.second = std::min(x.second, y.second);
            }
        };

        Info<< "non-contiguous reduce: " << flush;

        clockTime timing;
        for (label count = 0; count < repeat; ++count)
        {
            Pstream::listCombineReduce
            (
                list,
                bop,
                UPstream::msgType(),
                UPstream::worldComm
            );
        }

        double speed = timing.timeIncrement();

        int64_t memPeak = memInfo::get_hwm();
        Foam::reduce(memPeak, maxOp<int64_t>());

        Info<< speed << " sec, mem:" << memPeak << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
