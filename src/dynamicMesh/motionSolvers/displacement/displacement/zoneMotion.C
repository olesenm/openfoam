/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "zoneMotion.H"
#include "syncTools.H"
#include "bitSet.H"
#include "cellSet.H"
#include "cellZoneMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneMotion::selectZonePoints(const polyMesh& mesh)
{
    const auto& allCellZones = mesh.cellZones();

    if (!cellZoneIDs_.empty())
    {
        bitSet movePts(mesh.nPoints());

        // Markup points associated with cell zone(s)
        for (const label zonei : cellZoneIDs_)
        {
            for (const label celli : allCellZones[zonei])
            {
                for (const label facei : mesh.cells()[celli])
                {
                    movePts.set(mesh.faces()[facei]);
                }
            }
        }

        syncTools::syncPointList(mesh, movePts, orEqOp<unsigned int>(), 0u);

        pointIDs_ = movePts.sortedToc();

        // No cell points selected => move all points
        if (returnReduceAnd(pointIDs_.empty()))
        {
            cellZoneIDs_.clear();  // consistency
        }
    }

    selectionMode_ =
    (
        cellZoneIDs_.empty()
      ? selectionModes::smAll
      : selectionModes::smCellZone
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneMotion::zoneMotion()
:
    pointIDs_(),
    cellZoneIDs_(),
    selectionMode_(selectionModes::smAll)
{}



Foam::zoneMotion::zoneMotion
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneMotion()
{
    resetZone(mesh, dict);
}


Foam::zoneMotion::zoneMotion(const dictionary& dict, const polyMesh& mesh)
:
    zoneMotion()
{
    resetZone(mesh, dict);
}


Foam::zoneMotion::zoneMotion
(
    const polyMesh& mesh,
    const labelUList& cellZoneIds
)
:
    zoneMotion()
{
    cellZoneIDs_ = cellZoneIds;
    selectZonePoints(mesh);
}


Foam::zoneMotion::zoneMotion
(
    const polyMesh& mesh,
    const wordRe& cellZoneSelection
)
:
    zoneMotion()
{
    resetZone(mesh, cellZoneSelection);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneMotion::resetZone()
{
    pointIDs_.clear();
    cellZoneIDs_.clear();
    selectionMode_ = selectionModes::smAll;
}


void Foam::zoneMotion::resetZone
(
    const polyMesh& mesh,
    const wordRe& cellZoneSelection
)
{
    pointIDs_.clear();
    cellZoneIDs_.clear();
    selectionMode_ = selectionModes::smAll;

    if (!cellZoneSelection.empty())
    {
        const auto& allCellZones = mesh.cellZones();

        // Also handles groups, multiple zones (as wordRe match) ...
        cellZoneIDs_ = allCellZones.indices(cellZoneSelection);

        if (cellZoneIDs_.empty())
        {
            WarningInFunction
                << "No matching cellZones: " << cellZoneSelection << nl
                << "    Valid zones : "
                << flatOutput(allCellZones.names()) << nl
                << "    Valid groups: "
                << flatOutput(allCellZones.groupNames()) << nl;
        }
        else
        {
            selectZonePoints(mesh);
        }
    }
}


void Foam::zoneMotion::resetZone
(
    const polyMesh& mesh,
    const labelUList& cellZoneIds
)
{
    pointIDs_.clear();
    cellZoneIDs_ = cellZoneIds;

    if (cellZoneIDs_.empty())
    {
        selectionMode_ = selectionModes::smAll;
    }
    else
    {
        selectZonePoints(mesh);
    }
}


void Foam::zoneMotion::resetZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    pointIDs_.clear();
    cellZoneIDs_.clear();
    selectionMode_ = selectionModes::smAll;

    int useSubset(0);
    const auto& allCellZones = mesh.cellZones();

    constexpr bool verbose = true;

    // Specified cellSet?
    labelHashSet cellIDs;
    if
    (
        word cellSetName;
        (
            dict.readIfPresent("cellSet", cellSetName)
         && !cellSetName.empty()
         && cellSetName != "none"  // Compat: ignore 'none' placeholder
        )
    )
    {
        if (verbose)
        {
            Info<< "Applying motion to cellSet: " << cellSetName << endl;
        }
        useSubset |= int(selectionModes::smCellSet);
        cellIDs = cellSet::readContents(mesh, cellSetName);
    }

    // Specified cellZone(s) ?
    if
    (
        wordRe cellZoneSelection;
        (
            dict.readIfPresent("cellZone", cellZoneSelection)
         && !cellZoneSelection.empty()
         && cellZoneSelection != "none"  // Compat: ignore 'none' placeholder
        )
    )
    {
        if (verbose)
        {
            Info<< "Applying motion to cellZone: " << cellZoneSelection << endl;
        }
        useSubset |= int(selectionModes::smCellZone);

        // Also handles groups, multiple zones (as wordRe match) ...
        cellZoneIDs_ = allCellZones.indices(cellZoneSelection);

        if (cellZoneIDs_.empty())
        {
            FatalIOErrorInFunction(dict)
                << "No matching cellZones: " << cellZoneSelection << nl
                << "    Valid zones : "
                << flatOutput(allCellZones.names()) << nl
                << "    Valid groups: "
                << flatOutput(allCellZones.groupNames()) << nl
                << exit(FatalIOError);
        }
    }

    if (useSubset)
    {
        bitSet movePts(mesh.nPoints());

        // Markup points associated with cell zones
        for (const label zonei : cellZoneIDs_)
        {
            for (const label celli : mesh.cellZones()[zonei])
            {
                for (const label facei : mesh.cells()[celli])
                {
                    movePts.set(mesh.faces()[facei]);
                }
            }
        }

        // Markup points associated with cellSet
        for (const label celli : cellIDs)
        {
            for (const label facei : mesh.cells()[celli])
            {
                movePts.set(mesh.faces()[facei]);
            }
        }

        syncTools::syncPointList(mesh, movePts, orEqOp<unsigned int>(), 0u);

        pointIDs_ = movePts.sortedToc();
    }


    // No cell points selected => move all points
    if (returnReduceAnd(pointIDs_.empty()))
    {
        if (verbose)
        {
            Info<< "Applying motion to entire mesh" << endl;
        }
        useSubset = 0;
        cellZoneIDs_.clear();  // consistency
    }

    selectionMode_ = selectionModes(useSubset);
}


void Foam::zoneMotion::rebuildPointIDs
(
    const dictionary& dict,
    const polyMesh& mesh
)
{
    resetZone(mesh, dict);
}


// ************************************************************************* //
