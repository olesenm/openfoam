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
    Test-globalIndex4

Description
    Simple tests for the globalIndex class.

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"
#include "globalOffset.H"
#include "argList.H"
#include "IndirectList.H"
#include "IOstreams.H"
#include "Random.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();

    #include "setRootCase.H"


    Info<< "Dispatch tags: "
        "gatherOnly/gatherTotal/gatherOffset/gatherNonLocal\n";

    const auto check = [](bool ok, const char* what)
    {
        if (!ok)
        {
            FatalErrorInFunction
                << "Problem with dispatch tag: " << what
                << abort(FatalError);
        }
    };

    {
        const label nLocal = 100*UPstream::myProcNo();
        const bool is_master = UPstream::master();

        // Reference (all ranks have everything)
        const globalIndex reference(nLocal);
        const label total = reference.totalSize();
        const label start = reference.localStart();

        // gatherOnly : master has everything, sub-ranks nothing
        {
            const globalIndex gi(globalIndex::gatherOnly{}, nLocal);
            Pout<< "gatherOnly: "; gi.offsets().writeList(Pout) << nl;

            if (is_master)
            {
                check(gi.offsets() == reference.offsets(), "gatherOnly master");
            }
            else
            {
                check(gi.offsets().empty(), "gatherOnly sub-rank");
            }
        }

        // gatherNonLocal : as gatherOnly, but master size treated as 0
        {
            const globalIndex gi(globalIndex::gatherNonLocal{}, nLocal);
            Pout<< "gatherNonLocal: "; gi.offsets().writeList(Pout) << nl;

            if (is_master)
            {
                check(gi.localSize(0) == 0, "gatherNonLocal master size");
                check
                (
                    gi.totalSize() == total - nLocal || !UPstream::parRun(),
                    "gatherNonLocal total"
                );
            }
            else
            {
                check(gi.offsets().empty(), "gatherNonLocal sub-rank");
            }
        }

        // gatherTotal : master has everything, sub-ranks a single total
        {
            const globalIndex gi(globalIndex::gatherTotal{}, nLocal);
            Pout<< "gatherTotal: "; gi.offsets().writeList(Pout) << nl;

            if (is_master)
            {
                check(gi.offsets() == reference.offsets(), "gatherTotal master");
            }
            else
            {
                check(gi.offsets().size() == 1, "gatherTotal sub-rank size");
                check(gi.offsets()[0] == total, "gatherTotal sub-rank total");
            }
        }

        // gatherOffset : master has everything,
        //   sub-ranks (begin,end,total) or (begin,end)
        {
            const globalIndex gi(globalIndex::gatherOffset{}, nLocal);
            Pout<< "gatherOffset: "; gi.offsets().writeList(Pout) << nl;

            if (is_master)
            {
                check(gi.offsets() == reference.offsets(), "gatherOffset master");
            }
            else
            {
                check
                (
                    (gi.offsets().size() == 2 || gi.offsets().size() == 3),
                    "gatherOffset sub-rank size"
                );

                check(gi.offsets()[0] == start, "gatherOffset begin");
                check(gi.offsets()[1] == start + nLocal, "gatherOffset end");
                check(gi.offsets().back() == total, "gatherOffset total");
            }
        }

        // reset() with the same tags, re-using an existing object
        {
            globalIndex gi;
            gi.reset(globalIndex::gatherTotal{}, nLocal);
            check
            (
                gi.offsets().size() == (is_master ? reference.offsets().size() : 1),
                "reset gatherTotal"
            );

            gi.reset(globalIndex::gatherOffset{}, nLocal);

            if (is_master)
            {
                check(gi.offsets() == reference.offsets(), "reset gatherOffset master");
            }
            else
            {
                check
                (
                    (gi.offsets().size() == 2 || gi.offsets().size() == 3),
                    "reset gatherOffset sub-rank"
                );
            }
        }

        // Non-parallel communicator: purely local content
        {
            globalIndex gi;
            gi.reset(globalIndex::gatherOffset{}, nLocal, UPstream::commSelf());
            check(gi.totalSize() == nLocal, "commSelf gatherOffset");

            gi.reset(globalIndex::gatherTotal{}, nLocal, UPstream::commSelf());
            check(gi.totalSize() == nLocal, "commSelf gatherTotal");
        }
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
