/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "patchSeedSet.H"
#include "polyMesh.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "mappedPatchBase.H"
#include "triangulatedPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchSeedSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchSeedSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchSeedSet::calcPatchSamples
(
    const label nAvailable,
    const label nPatchPoints,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    if (nAvailable < 1)
    {
        return;
    }

    Random& rndGen = *rndGenPtr_;

    // Globally consistent addressing
    const label globalPointStart = UPstream::find_offset(nAvailable);

    const label nGlobalPatchPoints = returnReduce(nPatchPoints, sumOp<label>());

    point pt;
    label facei;
    label celli;

    const bool perturb = true;

    for (auto patchi : patchIDs_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        triangulatedPatch tp(pp, perturb);

        const label np = nAvailable*pp.size()/scalar(nGlobalPatchPoints);
        for (label i = 0; i < np; ++i)
        {
            tp.randomLocalPoint(rndGen, pt, facei, celli);

            samplingPts.append(pt);
            samplingCells.append(celli);
            samplingFaces.append(pp.start() + facei);
            samplingSegments.append(0);
            samplingCurveDist.append(globalPointStart + i);
        }
    }
}


void Foam::patchSeedSet::calcSelectedLocations
(
    const label nAvailable,
    const label nPatchPoints,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);

    if (nAvailable < 1)
    {
        return;
    }

    Random& rndGen = *rndGenPtr_;

    // Boundary faces only
    const label bndFaceOffset = mesh().nInternalFaces();
    bitSet bndFaces;
    treeBoundBox overallBb;

    // Collect faces and bounding box
    for (auto patchi : patchIDs_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        bndFaces.set(pp.offset_range());

        overallBb.add(pp.box());
    }

    // 1. Find nearest local patch face for all selectedLocations
    labelList patchFaces(nAvailable, -1);

    // The nearest selection key. Construct as missed sample
    // == (large distance, proc = max_ranks)

    List<Tuple2<scalar, int>> nearest
    (
        nAvailable,
        { VGREAT, UPstream::max_ranks() }
    );

    const bool hasLocal = bndFaces.any();

    if (hasLocal)
    {
        overallBb.inflate(rndGen, 1e-4, ROOTVSMALL);

        const indexedOctree<treeDataFace> boundaryTree
        (
            // Boundary faces only
            treeDataFace(mesh(), std::move(bndFaces), bndFaceOffset),

            overallBb,      // overall search domain
            8,              // maxLevel
            10,             // leafsize
            3.0             // duplicity
        );

        // Get some global dimension so all points are equally likely
        // to be found
        const scalar globalDistSqr
        (
            //overallBb.magSqr()
            GREAT
        );

        const auto& treeData = boundaryTree.shapes();
        const auto& samples = selectedLocations_;

        for (label samplei = 0; samplei < nAvailable; ++samplei)
        {
            const point& sample = samples[samplei];

            pointIndexHit pHit = boundaryTree.findNearest
            (
                sample,
                globalDistSqr
            );

            if (const auto& p = pHit.point(); pHit.hit())
            {
                label facei = treeData.objectIndex(pHit.index());

                patchFaces[samplei] = facei;
                nearest[samplei] = { sample.distSqr(p), myProci };
            }
        }
    }

    // Find rank with nearest hit - globally consistent
    if (UPstream::vendorCode() == UPstream::vendorCodes::intelmpi)
    {
        // Avoid slow MINLOC on intelmpi (OCT-2026)
        Pstream::listCombineReduce
        (
            nearest, minFirstEqOp{}, UPstream::msgType(), communicator
        );
    }
    else
    {
        UPstream::minLocation(nearest, communicator);
    }

    // Extract/filter local selections
    {
        label nLocal = 0;
        forAll(nearest, samplei)
        {
            if (myProci == nearest[samplei].second())
            {
                patchFaces[nLocal] = patchFaces[samplei];
                ++nLocal;
            }
        }
        patchFaces.resize(nLocal);
    }

    if (debug)
    {
        Pout<< "Found " << patchFaces.size()
            << " out of " << nAvailable
            << " on local processor" << endl;
    }


    // Shuffle and truncate if in random mode
    const label totalSize = returnReduce(patchFaces.size(), sumOp<label>());

    if (totalSize > nAvailable)
    {
        // Check what fraction of maxPoints_ I need to generate locally.
        label myMaxPoints = scalar(patchFaces.size())/totalSize*nAvailable;

        labelList subset = identity(patchFaces.size());
        for (label iter = 0; iter < 4; ++iter)
        {
            forAll(subset, i)
            {
                label j = rndGen.position<label>(0, subset.size()-1);
                std::swap(subset[i], subset[j]);
            }
        }

        // Truncate
        subset.setSize(myMaxPoints);

        // Subset patchFaces

        if (debug)
        {
            Pout<< "In random mode : selected " << subset.size()
                << " faces out of " << patchFaces.size() << endl;
        }

        patchFaces = labelUIndList(patchFaces, subset).list();
    }


    // Get points on patchFaces.

    // Globally consistent addressing
    const label globalFaceStart = UPstream::find_offset(patchFaces.size());

    samplingPts.setCapacity(patchFaces.size());
    samplingCells.setCapacity(patchFaces.size());
    samplingFaces.setCapacity(patchFaces.size());
    samplingSegments.setCapacity(patchFaces.size());
    samplingCurveDist.setCapacity(patchFaces.size());

    // For calculation of min-decomp tet base points
    (void)mesh().tetBasePtIs();

    forAll(patchFaces, i)
    {
        const label facei = patchFaces[i];
        const label celli = mesh().faceOwner()[facei];

        // Slightly shift point in since on warped face face-diagonal
        // decomposition might be outside cell for face-centre decomposition!
        pointIndexHit info = mappedPatchBase::facePoint
        (
            mesh(),
            facei,
            polyMesh::FACE_DIAG_TRIS
        );

        if (info.hit())
        {
            // Move the point into the cell
            const point& cc = mesh().cellCentres()[celli];
            samplingPts.append
            (
                info.point() + 1e-1*(cc-info.point())
            );
        }
        else
        {
            samplingPts.append(info.point());
        }
        samplingCells.append(celli);
        samplingFaces.append(facei);
        samplingSegments.append(0);
        samplingCurveDist.append(globalFaceStart+i);
    }
}


void Foam::patchSeedSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug)
    {
        write(Info);
    }
}


void Foam::patchSeedSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    DebugInfo << "patchSeedSet : sampling on patches :" << endl;

    if (!rndGenPtr_)
    {
        rndGenPtr_.reset(new Random(0));
    }

    label nPatchPoints = 0;
    for (auto patchi : patchIDs_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        nPatchPoints += pp.size();

        DebugInfo << "    " << pp.name() << " size " << pp.size() << endl;
    }

    label nAvailable = min(maxPoints_, selectedLocations_.size());

    calcSelectedLocations
    (
        nAvailable,
        nPatchPoints,
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    nAvailable = maxPoints_ - nAvailable;

    calcPatchSamples
    (
        nAvailable,
        nPatchPoints,
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchSeedSet::patchSeedSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    patchIDs_
    (
        mesh.boundaryMesh().patchSet(dict.get<wordRes>("patches")).sortedToc()
    ),
    maxPoints_(dict.get<label>("maxPoints")),
    selectedLocations_
    (
        dict.getOrDefault<pointField>
        (
            "points",
            pointField(0)
        )
    )
{
    genSamples();
}


// ************************************************************************* //
