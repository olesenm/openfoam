/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "subTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(subTriSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        subTriSurfaceMesh,
        dict
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::triSurface Foam::subTriSurfaceMesh::create_subset
(
    const IOobject& io,
    const dictionary& dict
)
{
    const word subGeomName(dict.get<word>("surface"));

    const wordRes includeNames(dict.get<wordRes>("patches"));

    const triSurfaceMesh& tsm =
        io.db().lookupObject<triSurfaceMesh>(subGeomName);

    if (tsm.type() == "distributedTriSurfaceMesh")
    {
        IOWarningInFunction(dict)
            << "Found surface " << subGeomName << " of type " << tsm.type()
            << ". subTriSurfaceMesh only works with non-distributed surfaces."
            << endl;
    }

    const auto& s = static_cast<const triSurface&>(tsm);

    // This is all just for extra error reporting
    {
        const auto& surfPatches = s.patches();

        // Verify that something matches
        bool anyMatches = false;

        for (const auto& p : surfPatches)
        {
            if (includeNames.match(p.name()))
            {
                anyMatches = true;
                break;
            }
        }

        if (!anyMatches)
        {
            auto& err = FatalIOErrorInFunction(dict);

            err << "Found no regions in triSurface matching "
                << includeNames << "\n    Regions: (";

            forAll(surfPatches, patchi)
            {
                if (patchi) err << ' ';
                err << surfPatches[patchi].name();
            }

            err << ")\n" << exit(FatalIOError);
        }
    }

    return s.subsetMesh(includeNames);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subTriSurfaceMesh::subTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, create_subset(io, dict))
{}


// ************************************************************************* //
