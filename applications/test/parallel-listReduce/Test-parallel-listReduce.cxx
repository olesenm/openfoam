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
    Test-parallel-listReduce

Description
    Test for simple listCombine reductions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "mapDistribute.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Random.H"
#include "Tuple2.H"
#include "mappedPatchBase.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("Increase Pstream::debug level");

    argList args(argc, argv);

    Pstream::debug += args.verbose();

    constexpr label listLen = 100;

    const auto myProci = UPstream::myProcNo();

    // Per-rank seed so that each process generates different values
    Random rndGen(Random::defaultSeed + myProci);


    // Simple listCombineReduce
    {
        typedef Tuple2<scalar, int32_t> Type;

        Pout<< nl << "type: " << error::demangle<Type>()
            << ", contiguous:" << is_contiguous_v<Type> << endl;

        List<Type> nearest(listLen, { GREAT, myProci });

        for (auto& near : nearest)
        {
            near.first() = rndGen.position<scalar>(50, 100);
        }

        Pout<< nl << "input: " << flatOutput(nearest) << endl;

        Pstream::listCombineReduce(nearest, minFirstEqOp{});

        Pout<< nl << "reduced: " << flatOutput(nearest) << endl;
    }

    // Simple listCombineReduce
    {
        typedef Tuple2<scalar, int64_t> Type;

        Pout<< nl << "type: " << error::demangle<Type>()
            << ", contiguous:" << is_contiguous_v<Type> << endl;

        List<Type> nearest(listLen, { GREAT, myProci });

        for (auto& near : nearest)
        {
            near.first() = rndGen.position<scalar>(50, 100);
        }

        Pout<< nl << "input: " << flatOutput(nearest) << endl;

        Pstream::listCombineReduce(nearest, minFirstEqOp{});

        Pout<< nl << "reduced: " << flatOutput(nearest) << endl;
    }

    // Simple listCombineReduce
    {
        typedef Tuple2<scalar, labelPair> Type;

        Pout<< nl << "type: " << error::demangle<Type>()
            << ", contiguous:" << is_contiguous_v<Type> << endl;

        List<Type> nearest(listLen, { GREAT, { myProci, -1 }} );

        for (auto& near : nearest)
        {
            near.first() = rndGen.position<scalar>(50, 100);
        }

        Pout<< nl << "input: " << flatOutput(nearest) << endl;

        Pstream::listCombineReduce(nearest, minFirstEqOp{});

        Pout<< nl << "reduced: " << flatOutput(nearest) << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
