/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

#include "patchFieldProbe.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchFieldProbe, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchFieldProbe::findElements(const fvMesh& mesh)
{
    const auto myProci = UPstream::myProcNo();

    (void)mesh.tetBasePtIs();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    // Nearest selection: (sqr(distance), (rank, facei))
    List<Tuple2<scalar, labelPair>> nearest
    (
        this->size(),
        { Foam::sqr(GREAT), {-1, -1} }
    );

    // Initialize all hit points as point::max
    pointField hitPoints(this->size(), point::max);

    patchIDs_ = bm.patchSet(patchNames_).sortedToc();

    // Boundary faces only
    const label bndFaceOffset = mesh.nInternalFaces();
    bitSet bndFaces;
    treeBoundBox overallBb;

    // Collect faces and bounding box
    for (label patchi : patchIDs_)
    {
        if (const polyPatch& pp = bm[patchi]; pp.size())
        {
            bndFaces.set(pp.offset_range());
            overallBb.add(pp.box());
        }
    }

    if (bndFaces.any())
    {
        constexpr label rand_seed = 123456;

        overallBb.inflate_random(rand_seed, 1e-4, ROOTVSMALL);

        const indexedOctree<treeDataFace> boundaryTree
        (
            // Boundary faces only
            treeDataFace(mesh, std::move(bndFaces), bndFaceOffset),

            overallBb,      // overall search domain
            8,              // maxLevel
            10,             // leafsize
            3.0             // duplicity
        );

        const auto& treeData = boundaryTree.shapes();

        const auto& samples = probeLocations();

        forAll(samples, probei)
        {
            const point& sample = samples[probei];
            auto& near = nearest[probei];
            auto& hitAt = hitPoints[probei];

            pointIndexHit pHit =
                boundaryTree.findNearest(sample, boundaryTree.bb().magSqr());

            if (!pHit.hit())
            {
                pHit = boundaryTree.findNearest(sample, Foam::sqr(GREAT));
            }

            if (const point& p = pHit.point(); pHit.hit())
            {
                const label facei = treeData.objectIndex(pHit.index());

                const label patchi = bm.whichPatch(facei);

                if (isA<emptyPolyPatch>(bm[patchi]))
                {
                    WarningInFunction
                        << " The sample point: " << sample
                        << " belongs to " << patchi
                        << " which is an empty patch. This is not permitted. "
                        << " This sample will not be included "
                        << endl;
                }
                else
                {
                    hitAt = p;
                    near = { sample.distSqr(p), {myProci, facei} };
                }
            }
        }
    }


    // Find nearest - globally consistent
    Pstream::listCombineReduce(nearest, minFirstEqOp{});

    // Update hitPoints, with extra sanity
    forAll(nearest, i)
    {
        const auto proci = nearest[i].second().first();
        const point& sample = probeLocations()[i];

        if (proci < 0)
        {
            // No hit anywhere, preserve original sample point
            hitPoints[i] = sample;
        }
        else if (proci != myProci)
        {
            // No hit on my proc
            hitPoints[i] = point::max;
        }
    }

    // Selects the best point, the original sample or point::max
    Foam::reduce
    (
        hitPoints.data(),
        hitPoints.size(),
        minOp<point>(),
        UPstream::msgType(),
        UPstream::worldComm
    );

    if (debug)
    {
        InfoInFunction << nl;

        const auto& samples = probeLocations();

        forAll(nearest, samplei)
        {
            const auto& [proci, facei] = nearest[samplei].second();
            const point& sample = samples[samplei];

            Info<< "    " << samplei << " coord:" << sample
                << " found on processor:" << proci
                << " in local face:" << facei
                << " with location:" << hitPoints[samplei]
                << endl;
        }
    }

    // Update actual probe locations and store old ones
    {
        pointField& locations = this->probeLocations();

        oldPoints_ = locations;
        locations = std::move(hitPoints);
    }

    // Extract any local faces to sample:
    // - operator[] : actual point to sample (=nearest point on patch)
    // - oldPoints_ : original provided point (might be anywhere in the mesh)
    // - cellIds_   : cells, not used
    // - faceIds_   : mesh face id
    // - patchIds_  : patch corresponding to faceList
    // - procIds_   : processor

    cellIds_.resize_fill(nearest.size(), -1);
    faceIds_.resize_fill(nearest.size(), -1);
    procIds_.resize_fill(nearest.size(), -1);
    patchIds_.resize_fill(nearest.size(), -1);

    forAll(nearest, samplei)
    {
        const auto& [proci, facei] = nearest[samplei].second();

        // Valid (local) face to sample?
        if (proci == myProci && facei >= 0)
        {
            faceIds_[samplei] = facei;
            procIds_[samplei] = proci;
            patchIds_[samplei] = bm.whichPatch(facei);
        }
    }

    // NB: do not reduce faceIds_ : they are processor-local

    Foam::reduce
    (
        procIds_.data(),
        procIds_.size(),
        maxOp<label>(),
        UPstream::msgType(),
        UPstream::worldComm
    );

    Foam::reduce
    (
        patchIds_.data(),
        patchIds_.size(),
        maxOp<label>(),
        UPstream::msgType(),
        UPstream::worldComm
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchFieldProbe::patchFieldProbe
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    probeModel(mesh, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchFieldProbe::read(const dictionary& dict)
{
    if (!probeModel::read(dict))
    {
        return false;
    }

    if (!dict.readIfPresent("patches", patchNames_))
    {
        patchNames_.resize(1);
        patchNames_.first() = dict.get<word>("patch");
    }

    // Initialise cells to sample from supplied locations
    findElements(thisMesh_);

    return true;
}


// ************************************************************************* //
