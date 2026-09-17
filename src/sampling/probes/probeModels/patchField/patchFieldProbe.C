/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
#include "mappedPatchBase.H"
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
    (void)mesh.tetBasePtIs();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    // All the info for nearest. Construct to miss
    List<mappedPatchBase::nearInfo> nearest(this->size());

    patchIDs_ = bm.patchSet(patchNames_).sortedToc();

    label nFaces = 0;
    forAll(patchIDs_, i)
    {
        nFaces += bm[patchIDs_[i]].size();
    }

    if (nFaces > 0)
    {
        // Collect mesh faces and bounding box
        labelList bndFaces(nFaces);
        treeBoundBox overallBb;

        nFaces = 0;
        forAll(patchIDs_, i)
        {
            const polyPatch& pp = bm[patchIDs_[i]];
            forAll(pp, i)
            {
                bndFaces[nFaces++] = pp.start()+i;
                const face& f = pp[i];

                // Without reduction.
                overallBb.add(pp.points(), f);
            }
        }

        constexpr label rand_seed = 123456;

        overallBb.inflate_random(rand_seed, 1e-4, ROOTVSMALL);


        const indexedOctree<treeDataFace> boundaryTree
        (
            treeDataFace(mesh, bndFaces),  // patch faces only

            overallBb,      // overall search domain
            8,              // maxLevel
            10,             // leafsize
            3.0             // duplicity
        );

        forAll(probeLocations(), probei)
        {
            const auto& treeData = boundaryTree.shapes();
            const point sample = probeLocations()[probei];

            pointIndexHit info = boundaryTree.findNearest
            (
                sample,
                Foam::sqr(boundaryTree.bb().mag())
            );

            if (!info.hit())
            {
                info = boundaryTree.findNearest(sample, Foam::sqr(GREAT));
            }

            const label facei = treeData.objectIndex(info.index());

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
            else if (info.hit())
            {
                // Note: do we store the face centre or the actual nearest?
                // We interpolate using the faceI only though (no
                // interpolation) so it does not actually matter much, just for
                // the location written to the header.

                //const point& facePt = mesh.faceCentres()[faceI];
                const point& facePt = info.point();

                mappedPatchBase::nearInfo sampleInfo;

                sampleInfo.first() = pointIndexHit(true, facePt, facei);

                sampleInfo.second().first() = facePt.distSqr(sample);
                sampleInfo.second().second() = Pstream::myProcNo();

                nearest[probei] = sampleInfo;
            }
        }
    }


    // Find nearest - globally consistent
    Pstream::listCombineReduce(nearest, mappedPatchBase::nearestEqOp());

    oldPoints_.resize(this->size());

    pointField& probeLocations = this->probeLocations();
    // Update actual probe locations and store old ones
    forAll(nearest, samplei)
    {
        oldPoints_[samplei] = probeLocations[samplei];
        probeLocations[samplei] = nearest[samplei].first().point();
    }

    if (debug)
    {
        InfoInFunction << nl;
        forAll(nearest, samplei)
        {
            label proci = nearest[samplei].second().second();
            label locali = nearest[samplei].first().index();

            Info<< "    " << samplei << " coord:"<< probeLocations[samplei]
                << " found on processor:" << proci
                << " in local face:" << locali
                << " with location:" << nearest[samplei].first().point()
                << endl;
        }
    }

    // Extract any local faces to sample:
    // - operator[] : actual point to sample (=nearest point on patch)
    // - oldPoints_ : original provided point (might be anywhere in the mesh)
    // - cellIds_   : cells, not used
    // - faceIds_      : faces (now patch faces)
    // - patchIds_   : patch corresponding to faceList
    // - procIds_     : processor
    cellIds_.resize_nocopy(nearest.size());
    cellIds_ = -1;

    faceIds_.resize_nocopy(nearest.size());
    faceIds_ = -1;

    procIds_.resize_nocopy(nearest.size());
    procIds_ = -1;

    patchIds_.resize_nocopy(nearest.size());
    patchIds_ = -1;

    forAll(nearest, sampleI)
    {
        procIds_[sampleI] = nearest[sampleI].second().second();

        if (nearest[sampleI].second().second() == Pstream::myProcNo())
        {
            // Store the face to sample
            faceIds_[sampleI] = nearest[sampleI].first().index();
            const label facei = faceIds_[sampleI];
            if (facei != -1)
            {
                procIds_[sampleI] = Pstream::myProcNo();
                patchIds_[sampleI] = bm.whichPatch(facei);
            }
        }
        reduce(procIds_[sampleI], maxOp<label>());
        reduce(patchIds_[sampleI], maxOp<label>());
    }
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
