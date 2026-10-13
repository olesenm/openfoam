/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "multiSolidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "cellZoneMesh.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        multiSolidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::multiSolidBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    points0MotionSolver(mesh, dict, typeName)
{
    motionControls_.resize(coeffDict().size());
    zoneMotions_.resize(coeffDict().size());

    const auto& allCellZones = mesh.cellZones();

    label zonei = 0;

    for (const entry& e : coeffDict())
    {
        if (const auto* dictptr = e.dictPtr())
        {
            const auto& subDict = *dictptr;
            const auto& cellZoneName = e.keyword();

            // Also handles groups, multiple zones (as wordRe match) ...
            labelList zoneIDs = mesh.cellZones().indices(cellZoneName);

            if (zoneIDs.empty())
            {
                FatalIOErrorInFunction(coeffDict())
                    << "No matching cellZones: " << cellZoneName << nl
                    << "    Valid zones : "
                    << flatOutput(allCellZones.names()) << nl
                    << "    Valid groups: "
                    << flatOutput(allCellZones.groupNames()) << nl
                    << exit(FatalIOError);
            }

            motionControls_.set
            (
                zonei,
                solidBodyMotionFunction::New(subDict, mesh.time())
            );

            // The points associated with cell zone(s)
            auto& zoneMove = zoneMotions_.emplace_set(zonei, mesh, zoneIDs);

            Info<< "Applying solid body motion "
                << motionControls_[zonei].type()
                << " to "
                << returnReduce(zoneMove.pointIDs().size(), sumOp<label>())
                << " points of cellZone " << cellZoneName << endl;

            ++zonei;
        }
    }

    motionControls_.resize(zonei);
    zoneMotions_.resize(zonei);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::multiSolidBodyMotionSolver::curPoints() const
{
    auto ttransformedPts = tmp<pointField>::New(mesh().points());
    auto& transformedPts = ttransformedPts.ref();

    forAll(motionControls_, zonei)
    {
        const labelUList& zonePoints = zoneMotions_[zonei].pointIDs();

        UIndirectList<point>(transformedPts, zonePoints) = transformPoints
        (
            motionControls_[zonei].transformation(),
            pointField(points0_, zonePoints)
        );
    }

    return ttransformedPts;
}


void Foam::multiSolidBodyMotionSolver::distribute
(
    const mapDistributePolyMesh& map
)
{
    map.distributePointData(points0_);

    const polyMesh& thisMesh = this->mesh();
    label zonei = 0;

    for (const entry& dEntry : coeffDict())
    {
        if (dEntry.isDict())
        {
            if (zonei < zoneMotions_.size())
            {
                zoneMotions_[zonei].resetZone(thisMesh, dEntry.keyword());
            }

            ++zonei;
        }
    }
}


// ************************************************************************* //
