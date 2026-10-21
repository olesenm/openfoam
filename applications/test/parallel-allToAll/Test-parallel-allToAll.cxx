/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application
    Test-parallel-allToAll

Description
    Test allToAll variants with single and multiple elements per rank

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "FixedList.H"
#include "List.H"
#include "Pair.H"
#include "Pstream.H"
#include "IOstreams.H"

using namespace Foam;

static label nFailed = 0;

// The (non-zero) value sent from src to dst as element k
template<class T>
T sendValue(label src, label dst, label k)
{
    return T(1 + 1000*src + 10*dst + k);
}

// For the consensus tests: some pairs send nothing (all zero)
bool isSilent(label src, label dst)
{
    return ((src + dst) % 3 == 2);
}


void check(const char* name, bool ok)
{
    if (!ok)
    {
        ++nFailed;
        Pout<< "FAILED: " << name << endl;
    }
}


template<class T>
void testFlat(const label count, const char* name)
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    const int msgTag = __LINE__;

    List<T> sends(count*numProc);
    List<T> recvs(count*numProc, T(-1));

    forAll(sends, i)
    {
        sends[i] = sendValue<T>(myProci, i / count, i % count);
    }

    // Exchange
    UPstream::allToAll_n(sends, recvs, count, communicator);

    bool ok = true;
    forAll(recvs, i)
    {
        ok = ok && (recvs[i] == sendValue<T>(i / count, myProci, i % count));
    }
    check(name, ok);

    // Consensus (non-zero values only)
    forAll(sends, i)
    {
        const label proci = i / count;
        sends[i] =
        (
            isSilent(myProci, proci)
          ? T(0)
          : sendValue<T>(myProci, proci, i % count)
        );
    }
    recvs = T(-1);

    UPstream::allToAllConsensus_n(sends, recvs, count, msgTag, communicator);

    ok = true;
    forAll(recvs, i)
    {
        const label proci = i / count;
        const T expected =
        (
            isSilent(proci, myProci) ? T(0) : sendValue<T>(proci, myProci, i % count)
        );
        ok = ok && (recvs[i] == expected);
    }
    check(name, ok);
}


// Compile-time and run-time test of the FixedList overloads
template<class T, unsigned N>
void testFixedList(const char* name)
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    const int msgTag = __LINE__;

    List<FixedList<T, N>> sends(numProc);
    List<FixedList<T, N>> recvs(numProc);

    // Plain exchange
    forAll(sends, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            sends[proci][k] = sendValue<T>(myProci, proci, k);
        }
        recvs[proci] = T(-1);
    }

    UPstream::allToAll(sends, recvs, communicator);

    bool ok = true;
    forAll(recvs, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            ok = ok && (recvs[proci][k] == sendValue<T>(proci, myProci, k));
        }
    }
    check(name, ok);

    // Consensus
    forAll(sends, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            sends[proci][k] =
            (
                isSilent(myProci, proci)
              ? T(0)
              : sendValue<T>(myProci, proci, k)
            );
        }
        recvs[proci] = T(-1);
    }

    UPstream::allToAllConsensus(sends, recvs, msgTag, communicator);

    ok = true;
    forAll(recvs, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            const T expected =
            (
                isSilent(proci, myProci)
              ? T(0)
              : sendValue<T>(proci, myProci, k)
            );
            ok = ok && (recvs[proci][k] == expected);
        }
    }
    check(name, ok);
}


// Compile-time and run-time test of the Pair overloads
template<class T>
void testPair(const char* name)
{
    constexpr unsigned N = 2;

    constexpr auto sentinel = Pair<T>(-1, -1);

    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);
    const int msgTag = __LINE__;

    List<Pair<T>> sends(numProc);
    List<Pair<T>> recvs(numProc);

    // Plain exchange
    forAll(sends, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            sends[proci][k] = sendValue<T>(myProci, proci, k);
        }
        recvs[proci] = sentinel;
    }

    UPstream::allToAll(sends, recvs, communicator);

    bool ok = true;
    forAll(recvs, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            ok = ok && (recvs[proci][k] == sendValue<T>(proci, myProci, k));
        }
    }
    check(name, ok);

    // Consensus
    forAll(sends, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            sends[proci][k] =
            (
                isSilent(myProci, proci)
              ? T(0)
              : sendValue<T>(myProci, proci, k)
            );
        }
        recvs[proci] = sentinel;
    }

    UPstream::allToAllConsensus(sends, recvs, msgTag, communicator);

    ok = true;
    forAll(recvs, proci)
    {
        for (label k = 0; k < label(N); ++k)
        {
            const T expected =
            (
                isSilent(proci, myProci)
              ? T(0)
              : sendValue<T>(proci, myProci, k)
            );
            ok = ok && (recvs[proci][k] == expected);
        }
    }
    check(name, ok);
}


// Single-element wrappers
template<class T = Foam::label>
void testSingle()
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);
    const int msgTag = __LINE__;

    List<T> sends(numProc);
    List<T> recvs(numProc, T(-1));

    forAll(sends, proci)
    {
        sends[proci] = sendValue<T>(myProci, proci, 0);
    }

    UPstream::allToAll(sends, recvs, communicator);

    bool ok = true;
    forAll(recvs, proci)
    {
        ok = ok && (recvs[proci] == sendValue<T>(proci, myProci, 0));
    }
    check("allToAll single", ok);

    forAll(sends, proci)
    {
        sends[proci] =
        (
            isSilent(myProci, proci) ? 0 : sendValue<T>(myProci, proci, 0)
        );
    }
    recvs = -1;

    UPstream::allToAllConsensus(sends, recvs, msgTag, communicator);

    ok = true;
    forAll(recvs, proci)
    {
        const T expected =
        (
            isSilent(proci, myProci) ? 0 : sendValue<T>(proci, myProci, 0)
        );
        ok = ok && (recvs[proci] == expected);
    }
    check("allToAllConsensus single", ok);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    testSingle<label>();

    testFlat<int32_t>(1, "int32 count=1");
    testFlat<int32_t>(3, "int32 count=3");
    testFlat<int64_t>(1, "int64 count=1");
    testFlat<int64_t>(4, "int64 count=4");

    testFixedList<int32_t, 1>("FixedList<int32,1>");
    testFixedList<int32_t, 3>("FixedList<int32,3>");
    testFixedList<int64_t, 2>("FixedList<int64,2>");

    testPair<int32_t>("Pair<int32>");
    testPair<int64_t>("Pair<int64>");

    // These do not compile:
    // testFlat<uint64_t>(4, "int64 count=4");
    // testFixedList<double, 4>("FixedList<double,4>");
    // testPair<word>("Pair<word>");


    Foam::reduce(nFailed, sumOp<label>());

    Info<< nl << (nFailed ? "FAILED" : "PASSED")
        << " (" << nFailed << " failures)" << nl << endl;

    return nFailed ? 1 : 0;
}

// ************************************************************************* //
