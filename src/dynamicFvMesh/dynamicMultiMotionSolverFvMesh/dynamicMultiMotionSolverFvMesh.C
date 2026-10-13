/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "dynamicMultiMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "mapDistributePolyMesh.H"
#include "motionSolver.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMultiMotionSolverFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMultiMotionSolverFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMultiMotionSolverFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMultiMotionSolverFvMesh::dynamicMultiMotionSolverFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicMultiMotionSolverFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    IOobject dynMeshDictIO
    (
        "dynamicMeshDict",
        time().constant(),
       *this,
        IOobjectOption::MUST_READ,  //<- MUST_READ for initial setup
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    dictionary dynDict(IOdictionary::readContents(dynMeshDictIO));
    const auto& dynamicMeshCoeffs = dynDict.subDict(typeName + "Coeffs");

    // NO_READ for further construction
    dynMeshDictIO.readOpt(IOobjectOption::NO_READ);

    motionSolvers_.resize(dynamicMeshCoeffs.size());
    zoneMotions_.resize(dynamicMeshCoeffs.size());

    const auto& allCellZones = this->cellZones();

    label zonei = 0;

    for (const entry& e : dynamicMeshCoeffs)
    {
        if (const auto* dictptr = e.dictPtr())
        {
            const auto& subDict = *dictptr;

            wordRe cellZoneName;
            subDict.readEntry("cellZone", cellZoneName);

            // Also handles groups, multiple zones (as wordRe match) ...
            labelList zoneIDs = allCellZones.indices(cellZoneName);

            if (zoneIDs.empty())
            {
                FatalIOErrorInFunction(dynamicMeshCoeffs)
                    << "No matching cellZones: " << cellZoneName << nl
                    << "    Valid zones : "
                    << flatOutput(allCellZones.names()) << nl
                    << "    Valid groups: "
                    << flatOutput(allCellZones.groupNames()) << nl
                    << exit(FatalIOError);
            }

            motionSolvers_.set
            (
                zonei,
                motionSolver::New
                (
                    *this,
                    IOdictionary(dynMeshDictIO, subDict)
                )
            );

            // The points associated with cell zone(s)
            auto& zoneMove = zoneMotions_.emplace_set(zonei, *this, zoneIDs);

            Info<< "Applying motionSolver " << motionSolvers_[zonei].type()
                << " to "
                << returnReduce(zoneMove.pointIDs().size(), sumOp<label>())
                << " points of cellZone " << cellZoneName << endl;

            ++zonei;
        }
    }

    motionSolvers_.resize(zonei);
    zoneMotions_.resize(zonei);

    // Assume changed ...
    return true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Non-default in header (incomplete types)
Foam::dynamicMultiMotionSolverFvMesh::
~dynamicMultiMotionSolverFvMesh() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicMultiMotionSolverFvMesh::distribute
(
    const mapDistributePolyMesh& map
)
{
    forAll(motionSolvers_, zonei)
    {
        motionSolvers_[zonei].distribute(map);

        zoneMotions_[zonei].resetZone
        (
            *this,
            zoneMotions_[zonei].cellZoneIDs()
        );
    }
}


bool Foam::dynamicMultiMotionSolverFvMesh::update()
{
    pointField transformedPts(this->points());

    forAll(motionSolvers_, zonei)
    {
        const labelUList& zonePoints = zoneMotions_[zonei].pointIDs();

        const pointField newPoints(motionSolvers_[zonei].newPoints());

        for (const label pointi : zonePoints)
        {
            transformedPts[pointi] = newPoints[pointi];
        }
    }

    fvMesh::movePoints(transformedPts);

    static bool hasWarned = false;

    if (auto* Uptr = getObjectPtr<volVectorField>("U"))
    {
        Uptr->correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningInFunction
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //
