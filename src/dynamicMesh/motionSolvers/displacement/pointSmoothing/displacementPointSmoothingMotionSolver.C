/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenCFD Ltd.
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

#include "displacementPointSmoothingMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "mapDistributePolyMesh.H"
#include "syncTools.H"
#include "pointConstraints.H"
#include "motionSmootherAlgo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementPointSmoothingMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementPointSmoothingMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        displacementPointSmoothingMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::displacementPointSmoothingMotionSolver::markAffectedFaces
(
    const labelHashSet& changedFaces,
    labelHashSet& affectedFaces
)
{
    bitSet affectedPoints(mesh().nPoints());

    for (const label facei : changedFaces)
    {
        affectedPoints.set(mesh().faces()[facei]);
    }

    syncTools::syncPointList
    (
        mesh(),
        affectedPoints,
        orEqOp<unsigned int>(),
        0U
    );

    for (const label pointi : affectedPoints)
    {
        for (const label celli : mesh().pointCells()[pointi])
        {
            affectedFaces.insert(mesh().cells()[celli]);
        }
    }
}


bool Foam::displacementPointSmoothingMotionSolver::relax()
{
    if
    (
        relaxationFactors_.empty()
     || (relaxationFactors_.size() == 1 && relaxationFactors_[0] == 1.0)
    )
    {
        relaxedPoints_ = points0() + pointDisplacement().internalField();
        return true;
    }


    const pointField oldRelaxedPoints(relaxedPoints_);

    labelHashSet affectedFaces(facesToMove_);

    // Create a list of relaxation levels
    // -1 indicates a point which is not to be moved
    //  0 is the starting value for a moving point
    labelList relaxationLevel(mesh().nPoints(), -1);

    for (const label facei : affectedFaces)
    {
        for (const label pointi : mesh().faces()[facei])
        {
            relaxationLevel[pointi] = 0;
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        relaxationLevel,
        maxEqOp<label>(),
        label(-1)
    );

    // Loop whilst relaxation levels are being incremented
    bool complete(false);
    while (!complete)
    {
        //scalar nAffectedFaces(affectedFaces.size());
        //reduce(nAffectedFaces, sumOp<scalar>());
        //Info << "    Moving " << nAffectedFaces << " faces" << endl;

        // Move the points
        forAll(relaxationLevel, pointI)
        {
            if (relaxationLevel[pointI] >= 0)
            {
                const scalar x
                (
                    relaxationFactors_[relaxationLevel[pointI]]
                );

                relaxedPoints_[pointI] =
                    (1 - x)*oldRelaxedPoints[pointI]
                  + x*(points0()[pointI] + pointDisplacement()[pointI]);
            }
        }

        // Get a list of changed faces
        labelHashSet markedFaces;
        markAffectedFaces(affectedFaces, markedFaces);
        labelList markedFacesList(markedFaces.toc());

        // Update the geometry
        meshGeometry_.correct(relaxedPoints_, markedFacesList);

        // Check the modified face quality
        markedFaces.clear();
        motionSmootherAlgo::checkMesh
        (
            false,
            meshQualityDict_,
            meshGeometry_,
            relaxedPoints_,
            markedFacesList,
            markedFaces
        );

        // Mark the affected faces
        affectedFaces.clear();
        markAffectedFaces(markedFaces, affectedFaces);

        // Increase relaxation and check convergence
        bitSet pointsToRelax(mesh().nPoints());
        complete = true;
        for (const label facei : affectedFaces)
        {
            pointsToRelax.set(mesh().faces()[facei]);
        }

        for (const label pointi : pointsToRelax)
        {
            if (relaxationLevel[pointi] < relaxationFactors_.size() - 1)
            {
                ++ relaxationLevel[pointi];
                complete = false;
            }
        }

        // Synchronise relaxation levels
        syncTools::syncPointList
        (
            mesh(),
            relaxationLevel,
            maxEqOp<label>(),
            label(0)
        );

        // Synchronise completion
        UPstream::reduceAnd(complete);
    }

    // Check for convergence
    bool converged(true);
    forAll(mesh().faces(), facei)
    {
        const face& fPoints(mesh().faces()[facei]);

        for (const label pointi :  fPoints)
        {
            if (relaxationLevel[pointi] > 0)
            {
                facesToMove_.insert(facei);

                converged = false;

                break;
            }
        }
    }

    // Syncronise convergence
    UPstream::reduceAnd(converged);

    //if (converged)
    //{
    //    Info<< "... Converged" << endl << endl;
    //}
    //else
    //{
    //    Info<< "... Not converged" << endl << endl;
    //}

    return converged;
}


void Foam::displacementPointSmoothingMotionSolver::setFacesToMove
(
    const dictionary& dict
)
{
    if (dict.getOrDefault<bool>("moveInternalFaces", true))
    {
        facesToMove_.resize(2*mesh().nFaces());
        forAll(mesh().faces(), faceI)
        {
            facesToMove_.insert(faceI);
        }
    }
    else
    {
        facesToMove_.resize(2*(mesh().nBoundaryFaces()));
        for
        (
            label faceI = mesh().nInternalFaces();
            faceI < mesh().nFaces();
            ++ faceI
        )
        {
            facesToMove_.insert(faceI);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementPointSmoothingMotionSolver::
displacementPointSmoothingMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    meshGeometry_(mesh),
    pointSmoother_(pointSmoother::New(mesh, coeffDict())),
    nPointSmootherIter_
    (
        coeffDict().get<label>("nPointSmootherIter")
    ),
    relaxedPoints_(mesh.points())
{
    if (coeffDict().readIfPresent("relaxationFactors", relaxationFactors_))
    {
        meshQualityDict_ = coeffDict().subDict("meshQuality");
    }
    setFacesToMove(coeffDict());
}


Foam::displacementPointSmoothingMotionSolver::
displacementPointSmoothingMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    meshGeometry_(mesh),
    pointSmoother_
    (
        pointSmoother::New
        (
            mesh,
            coeffDict()
        )
    ),
    nPointSmootherIter_
    (
        coeffDict().get<label>("nPointSmootherIter")
    ),
    relaxedPoints_(mesh.points())
{
    if (coeffDict().readIfPresent("relaxationFactors", relaxationFactors_))
    {
        meshQualityDict_ = coeffDict().subDict("meshQuality");
    }
    setFacesToMove(coeffDict());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Non-default in header (incomplete types)
Foam::displacementPointSmoothingMotionSolver::
~displacementPointSmoothingMotionSolver() = default;


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementPointSmoothingMotionSolver::curPoints() const
{
    //Note: twoDCorrect already done by ::solve

    return relaxedPoints_;
}


void Foam::displacementPointSmoothingMotionSolver::distribute
(
    const mapDistributePolyMesh& map
)
{
    displacementMotionSolver::distribute(map);
    map.distributePointData(relaxedPoints_);

    pointSmoother_ = pointSmoother::New(mesh(), coeffDict());
    facesToMove_.clear();
    setFacesToMove(coeffDict());
    meshGeometry_.correct();
}


void Foam::displacementPointSmoothingMotionSolver::solve()
{
    //Pout<< "time:" << mesh().time().timeName()
    //    << " mesh faceCentres:" << gAverage(mesh().faceCentres())
    //    << " mesh cellCentres:" << gAverage(mesh().cellCentres())
    //    << endl;

    movePoints(curPoints());

    // Update values on pointDisplacement
    pointDisplacement().boundaryFieldRef().updateCoeffs();

    // Extend: face-to-point-to-cell-to-faces
    labelHashSet affectedFaces;
    markAffectedFaces(facesToMove_, affectedFaces);

    for(label i = 0; i < nPointSmootherIter_; i ++)
    {
        meshGeometry_.correct
        (
            points0() + pointDisplacement().internalField(),
            affectedFaces.toc()
        );
        //Pout<< "iter:" << i
        //    << " faceCentres:" << gAverage(meshGeometry_.faceCentres())
        //    << " cellCentres:" << gAverage(meshGeometry_.cellCentres())
        //    << endl;

        pointSmoother_->update
        (
            affectedFaces.toc(),
            points0(),
            points0() + pointDisplacement().internalField(),
            meshGeometry_,
            pointDisplacement()
        );
    }

    relax();

    twoDCorrectPoints(relaxedPoints_);

    // Update pointDisplacement for actual relaxedPoints. Keep fixed-value
    // bcs.
    pointDisplacement().primitiveFieldRef() = relaxedPoints_-points0();

    // Adhere to multi-point constraints
    const pointConstraints& pcs =
         pointConstraints::New(pointDisplacement().mesh());
    pcs.constrainDisplacement(pointDisplacement(), false);

    // Update relaxedPoints to take constraints into account
    relaxedPoints_ = points0() + pointDisplacement().internalField();
}


// ************************************************************************* //
