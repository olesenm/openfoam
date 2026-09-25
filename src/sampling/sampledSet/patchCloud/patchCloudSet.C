/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "patchCloudSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "Time.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchCloudSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchCloudSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchCloudSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    const auto myRank = UPstream::myProcNo();

    DebugInfo << "patchCloudSet : sampling on patches :" << endl;

    // Boundary faces only
    const label bndFaceOffset = mesh().nInternalFaces();
    bitSet bndFaces;
    treeBoundBox overallBb;

    // Collect faces and bounding box
    for (label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        DebugInfo << "    " << pp.name() << " size " << pp.size() << endl;

        bndFaces.set(pp.offset_range());

        overallBb.add(pp.box());
    }

    // Force calculation of face-diagonal decomposition
    (void)mesh().tetBasePtIs();

    const auto myProci = UPstream::myProcNo();

    // Nearest selection: (sqr(distance), (rank, index))
    List<Tuple2<scalar, labelPair>> nearest
    (
        sampleCoords_.size(),
        { Foam::sqr(GREAT), {-1, -1} }
    );

    pointField hitPoints(sampleCoords_.size(), point::max);

    // Find the nearest locally
    if (bndFaces.any())
    {
        // Not very random
        constexpr label rand_seed = 123456;

        // Make bb asymmetric just to avoid problems on symmetric meshes
        // Make sure bb is 3D.
        overallBb.inflate_random(rand_seed, 1e-4, ROOTVSMALL);

        const indexedOctree<treeDataFace> patchTree
        (
            // Boundary faces only
            treeDataFace(mesh(), std::move(bndFaces), bndFaceOffset),
            overallBb,      // overall search domain
            8,              // maxLevel
            10,             // leafsize
            3.0             // duplicity
        );

        // Search limit for findNearest
        const scalar maxDistSqr = Foam::sqr(searchDist_);

        const auto& treeData = patchTree.shapes();

        forAll(sampleCoords_, samplei)
        {
            const point& sample = sampleCoords_[samplei];
            auto& near = nearest[samplei];
            auto& hitAt = hitPoints[samplei];

            pointIndexHit pHit = patchTree.findNearest(sample, maxDistSqr);

            // Set nearest index to mesh face label
            if (const auto& p = pHit.point(); pHit.hit())
            {
                const label facei = treeData.objectIndex(pHit.index());

                hitAt = p;

                near = { sample.distSqr(p), {myProci, facei} };
            }
        }
    }

    // Find nearest - globally consistent
    Pstream::listCombineReduce(nearest, minFirstEqOp{});

    forAll(nearest, i)
    {
        const auto& rank = nearest[i].second().first();

        if (rank != myRank)
        {
            // My rank provides nothing
            hitPoints[i] = point::max;
        }
    }

    // Selects the best point or point::max
    Foam::reduce
    (
        hitPoints.data(),
        hitPoints.size(),
        minOp<point>(),
        UPstream::msgType(),
        UPstream::worldComm
    );

    if (debug && UPstream::master())
    {
        OFstream str
        (
            mesh().time().path()
          / name()
          + "_nearest.obj"
        );
        Info<< "Dumping mapping as lines from supplied points to"
            << " nearest patch face to file " << str.name() << endl;

        label vertexIndex(0);
        forAll(nearest, i)
        {
            if (const auto rank = nearest[i].second().first(); rank >= 0)
            {
                meshTools::writeOBJ
                (
                    str,
                    sampleCoords_[i],
                    hitPoints[i],
                    vertexIndex
                );
            }
        }
    }


    // Store the sampling locations on the nearest processor
    forAll(nearest, samplei)
    {
        const auto& [rank, facei] = nearest[samplei].second();
        const auto& sample = sampleCoords_[samplei];
        const auto& p = hitPoints[samplei];

        if (rank >= 0)  // hit
        {
            if (rank == myProci)
            {
                samplingPts.push_back(p);
                samplingCells.push_back(mesh().faceOwner()[facei]);
                samplingFaces.push_back(facei);
                samplingSegments.push_back(0);
                samplingCurveDist.push_back(1.0 * samplei);
            }
        }
        else
        {
            // No processor found point near enough. Mark with special value
            // which is intercepted when interpolating
            if (UPstream::master())
            {
                samplingPts.push_back(sample);
                samplingCells.push_back(-1);
                samplingFaces.push_back(-1);
                samplingSegments.push_back(0);
                samplingCurveDist.push_back(1.0 * samplei);
            }
        }
    }
}


void Foam::patchCloudSet::genSamples()
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

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchCloudSet::patchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords,
    const labelHashSet& patchSet,
    const scalar searchDist
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords),
    patchSet_(patchSet),
    searchDist_(searchDist)
{
    genSamples();
}


Foam::patchCloudSet::patchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.get<pointField>("points")),
    patchSet_
    (
        mesh.boundaryMesh().patchSet(dict.get<wordRes>("patches"))
    ),
    searchDist_(dict.get<scalar>("maxDistance"))
{
    genSamples();
}


// ************************************************************************* //
