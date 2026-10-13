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

#include "displacementSmartPointSmoothingMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "mapDistributePolyMesh.H"
#include "syncTools.H"
#include "pointConstraints.H"
#include "motionSmootherAlgo.H"

//#include "fvMesh.H"
//#include "fvGeometryScheme.H"
#include "OBJstream.H"
#include "emptyPointPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementSmartPointSmoothingMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementSmartPointSmoothingMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        displacementSmartPointSmoothingMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::displacementSmartPointSmoothingMotionSolver::markAffectedFaces
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


bool Foam::displacementSmartPointSmoothingMotionSolver::relax()
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
        if (false)
        {
            // Use snappyHexMesh compatible checks
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
        }
        else
        {
            // Use pointSmoother specific
            tmp<scalarField> tfaceQ
            (
                pointUntangler_->faceQuality
                (
                    relaxedPoints_,
                    meshGeometry_.faceCentres(),
                    meshGeometry_.faceAreas(),
                    meshGeometry_.cellCentres(),
                    meshGeometry_.cellVolumes()
                )
            );

            if (debug)
            {
                MinMax<scalar> range(gMinMax(tfaceQ()));
                Pout<< "    min:" << range.min() << nl
                    << "    max:" << range.max() << endl;
            }

            labelList order;
            Foam::sortedOrder(tfaceQ(), order);

            label nUntangle = 0;
            forAll(order, i)
            {
                if (tfaceQ()[order[i]] > untangleQ_)
                {
                    nUntangle = i;
                    break;
                }
            }

            affectedFaces = labelList(SubList<label>(order, nUntangle));
        }

        // Increase relaxation and check convergence
        bitSet pointsToRelax(mesh().nPoints());
        complete = true;
        for (const label facei : affectedFaces)
        {
            pointsToRelax.set(mesh().faces()[facei]);
        }

        for (const label pointI : pointsToRelax)
        {
            if (relaxationLevel[pointI] < relaxationFactors_.size() - 1)
            {
                ++ relaxationLevel[pointI];

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
        for (const label pointi : mesh().faces()[facei])
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


void Foam::displacementSmartPointSmoothingMotionSolver::setFacesToMove
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


void Foam::displacementSmartPointSmoothingMotionSolver::emptyCorrectPoints
(
    pointVectorField& pointDisplacement
)
{
    // Assume empty point patches are already in correct location
    // so knock out any off-plane displacement.
    auto& fld = pointDisplacement.primitiveFieldRef();
    for (const auto& ppf : pointDisplacement.boundaryField())
    {
        if (isA<emptyPointPatchVectorField>(ppf))
        {
            const auto& mp = ppf.patch().meshPoints();
            forAll(mp, i)
            {
                pointConstraint pc;
                ppf.patch().applyConstraint(i, pc);
                fld[mp[i]] = pc.constrainDisplacement(fld[mp[i]]);
            }
        }
    }

    pointField wantedPoints(points0() + fld);
    twoDCorrectPoints(wantedPoints);
    fld = wantedPoints-points0();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementSmartPointSmoothingMotionSolver::
displacementSmartPointSmoothingMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    meshGeometry_(mesh),
    pointUntangler_
    (
        pointSmoother::New
        (
            coeffDict().get<word>("untangler"),
            mesh,
            coeffDict()
        )
    ),
    untangleQ_(coeffDict().get<scalar>("untangleQ")),
    minQ_(coeffDict().get<scalar>("minQ")),
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


Foam::displacementSmartPointSmoothingMotionSolver::
displacementSmartPointSmoothingMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    meshGeometry_(mesh),
    pointUntangler_
    (
        pointSmoother::New
        (
            coeffDict().get<word>("untangler"),
            mesh,
            coeffDict()
        )
    ),
    untangleQ_(coeffDict().get<scalar>("untangleQ")),
    minQ_(coeffDict().get<scalar>("minQ")),
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
Foam::displacementSmartPointSmoothingMotionSolver::
~displacementSmartPointSmoothingMotionSolver() = default;


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementSmartPointSmoothingMotionSolver::curPoints() const
{
    //Note: twoDCorrect already done by ::solve

    return relaxedPoints_;
}


void Foam::displacementSmartPointSmoothingMotionSolver::distribute
(
    const mapDistributePolyMesh& map
)
{
    displacementMotionSolver::distribute(map);
    map.distributePointData(relaxedPoints_);

    pointUntangler_ =
        pointSmoother::New
        (
            coeffDict().get<word>("untangler"),
            mesh(),
            coeffDict()
        );
    pointSmoother_ = pointSmoother::New(mesh(), coeffDict());
    facesToMove_.clear();
    setFacesToMove(coeffDict());
    meshGeometry_.correct();
}


void Foam::displacementSmartPointSmoothingMotionSolver::solve()
{
    movePoints(curPoints());

    // Update values on pointDisplacement. Note: should also evaluate? Since
    // e.g. cellMotionBC uses pointDisplacement value.
    pointDisplacement().boundaryFieldRef().updateCoeffs();


    fileName debugDir;
    if (debug & 2)
    {
        debugDir = mesh().time().timePath();
        mkDir(debugDir);
        OBJstream os(debugDir/"bc.obj");

        const pointField wantedPoints
        (
            points0()
          + pointDisplacement().internalField()
        );
        const auto& pbm = pointDisplacement().mesh().boundary();
        for (const auto& ppp : pbm)
        {
            if (!isA<emptyPointPatch>(ppp))
            {
                const auto& mp = ppp.meshPoints();
                for (const label pointi : mp)
                {
                    os.write
                    (
                        linePointRef
                        (
                            points0()[pointi],
                            wantedPoints[pointi]
                        )
                    );
                }
            }
        }
        Pout<< "Written " << os.nVertices() << " initial displacements to "
            << os.name() << endl;
    }


    // Extend: face-to-point-to-cell-to-faces
    labelHashSet affectedFaces;
    markAffectedFaces(facesToMove_, affectedFaces);


    for (label nIter = 0; nIter < nPointSmootherIter_; ++nIter)
    {
        const pointField wantedPoints
        (
            points0()
          + pointDisplacement().internalField()
        );

        meshGeometry_.correct
        (
            wantedPoints,
            affectedFaces.toc()
        );

        //{
        //    // Debugging: check meshGeometry consistent with fvGeometryScheme
        //    const auto& geom =
        //        reinterpret_cast<const fvMesh&>(mesh()).geometry();
        //    pointField faceCentres(mesh().nFaces());
        //    vectorField faceAreas(mesh().nFaces());
        //    pointField cellCentres(mesh().nCells());
        //    scalarField cellVolumes(mesh().nCells());
        //    geom.updateGeom
        //    (
        //        wantedPoints,
        //        mesh().points(),    // old points
        //        faceCentres,
        //        faceAreas,
        //        cellCentres,
        //        cellVolumes
        //    );
        //    forAll(faceCentres, facei)
        //    {
        //        const point& meshFc = mesh().faceCentres()[facei];
        //        const point& meshGeomFc = meshGeometry_.faceCentres()[facei];
        //        const point& updatedFc = faceCentres[facei];
        //
        //        if (updatedFc != meshGeomFc)
        //        {
        //            const face& f = mesh().faces()[facei];
        //
        //            Pout<< "At face:" << facei << nl
        //                << "    old         :" << meshFc << nl
        //                << "    new         :" << updatedFc << nl
        //                << "    polyMeshGeom:" << meshGeomFc << nl
        //                << "    oldPoints   :"
        //                << UIndirectList<point>(mesh().points(), f) << nl
        //                << "    wantedPoints:"
        //                << UIndirectList<point>(wantedPoints, f) << nl
        //                << endl;
        //        }
        //    }
        //}

        // Get measure of face quality
        tmp<scalarField> tfaceQ
        (
            pointUntangler_->faceQuality
            (
                wantedPoints,
                meshGeometry_.faceCentres(),
                meshGeometry_.faceAreas(),
                meshGeometry_.cellCentres(),
                meshGeometry_.cellVolumes()
            )
        );


        if (debug)
        {
            MinMax<scalar> range(gMinMax(tfaceQ()));
            Pout<< "    min:" << range.min() << nl
                << "    max:" << range.max() << endl;
        }

        labelList order;
        Foam::sortedOrder(tfaceQ(), order);

        label nUntangle = 0;
        forAll(order, i)
        {
            if (tfaceQ()[order[i]] > untangleQ_)
            {
                nUntangle = i;
                break;
            }
        }
        label nLow = 0;
        forAll(order, i)
        {
            if (tfaceQ()[order[i]] > minQ_)
            {
                nLow = i;
                break;
            }
        }


        if (debug)
        {
            Pout<< "    nUntangle:" << returnReduce(nUntangle, sumOp<label>())
                << nl
                << "    nLow     :" << returnReduce(nLow, sumOp<label>())
                << nl;
        }


        if (returnReduce(nUntangle, sumOp<label>()))
        {
            // Start untangling
            labelList lowQFaces(SubList<label>(order, nUntangle));
            //{
            //    // Grow set (non parallel)
            //    bitSet isMarkedFace(mesh().nFaces());
            //    for (const label facei : lowQFaces)
            //    {
            //        for (const label pointi : mesh().faces()[facei])
            //        {
            //            isMarkedFace.set(mesh().pointFaces()[pointi]);
            //        }
            //    }
            //    lowQFaces = isMarkedFace.sortedToc();
            //}

            //Pout<< "    untangling "
            //    << returnReduce(lowQFaces.size(), sumOp<label>())
            //    << " faces" << endl;
            pointUntangler_->update
            (
                lowQFaces,
                points0(),
                wantedPoints,
                meshGeometry_,
                pointDisplacement()
                //false                       // ! do NOT apply bcs, constraints
            );

            // Keep points on empty patches. Note: since pointConstraints
            // does not implement constraints on emptyPointPatches and
            // emptyPointPatchField does not either.
            emptyCorrectPoints(pointDisplacement());

            if (debug & 2)
            {
                OBJstream os(debugDir/"untangle_" + Foam::name(nIter) + ".obj");

                const pointField wantedPoints
                (
                    points0()
                  + pointDisplacement().internalField()
                );
                forAll(wantedPoints, pointi)
                {
                    os.write
                    (
                        linePointRef
                        (
                            points0()[pointi],
                            wantedPoints[pointi]
                        )
                    );
                }
                Pout<< "Written " << os.nVertices() << " wanted untangle to "
                    << os.name() << endl;
            }
        }
        else if (returnReduce(nLow, sumOp<label>()))
        {
            labelList lowQFaces(SubList<label>(order, nLow));
            //{
            //    // Grow set (non parallel)
            //    bitSet isMarkedFace(mesh().nFaces());
            //    for (const label facei : lowQFaces)
            //    {
            //        for (const label pointi : mesh().faces()[facei])
            //        {
            //            isMarkedFace.set(mesh().pointFaces()[pointi]);
            //        }
            //    }
            //    lowQFaces = isMarkedFace.sortedToc();
            //}

            //Pout<< "    smoothing "
            //    << returnReduce(lowQFaces.size(), sumOp<label>())
            //    << " faces" << endl;

            pointSmoother_->update
            (
                lowQFaces,
                points0(),
                wantedPoints,
                meshGeometry_,
                pointDisplacement()
            );
            // Keep points on empty patches
            emptyCorrectPoints(pointDisplacement());
        }
        else
        {
            //Pout<< "** converged" << endl;
            break;
        }
    }


    relax();
    //relaxedPoints_ = points0() + pointDisplacement().internalField();

    twoDCorrectPoints(relaxedPoints_);

    // Update pointDisplacement for actual relaxedPoints. Keep fixed-value
    // bcs.
    pointDisplacement().primitiveFieldRef() = relaxedPoints_-points0();

    // Adhere to multi-point constraints. Does correctBoundaryConditions +
    // multi-patch issues.
    const pointConstraints& pcs =
         pointConstraints::New(pointDisplacement().mesh());
    pcs.constrainDisplacement(pointDisplacement(), false);
}


// ************************************************************************* //
