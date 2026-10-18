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
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // Test for valid proc number
    const auto valid_proc = [numProc](int proci) noexcept
    {
        return (proci >= 0 && proci < numProc);
    };

    (void)mesh.tetBasePtIs();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    // The extracted sampling information (output):
    // - cellIds_   : cells, not used
    // - faceIds_   : mesh face id
    // - patchIds_  : patch corresponding to faceIds_
    // - procIds_   : processor
    //
    // NOTE: faceIds_ and patchIds_ are filled with local
    // information and re-filtered later

    cellIds_.resize_fill(this->size(), -1);
    faceIds_.resize_fill(this->size(), -1);
    patchIds_.resize_fill(this->size(), -1);
    procIds_.resize_fill(this->size(), -1);

    // Nearest selection key: (sqr(distance), rank)
    // Construct as missed sample
    // == (large distance, proc = max_ranks)
    List<Tuple2<scalar, int>> nearest
    (
        this->size(),
        { VGREAT, UPstream::max_ranks() }
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

            pointIndexHit pHit =
                boundaryTree.findNearest(sample, boundaryTree.bb().magSqr());

            if (pHit.not_hit())
            {
                pHit = boundaryTree.findNearest(sample, Foam::sqr(GREAT));
            }

            if (const point& p = pHit.point(); pHit.hit())
            {
                label facei = treeData.objectIndex(pHit.index());
                label patchi = bm.whichPatch(facei);

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
                    faceIds_[probei] = facei;
                    patchIds_[probei] = patchi;
                    hitPoints[probei] = p;
                    nearest[probei] = { sample.distSqr(p), myProci };
                }
            }
        }
    }


    // Find nearest - globally consistent
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

    // Update hitPoints, with extra sanity
    forAll(nearest, i)
    {
        const auto proci = nearest[i].second();
        const point& sample = probeLocations()[i];

        // No hit on my proc
        // - set as point::max (for later reduction)
        //   or preserve original sample point if no hit anywhere

        if (proci != myProci)
        {
            hitPoints[i] = valid_proc(proci) ? point::max : sample;
        }
    }

    // Selects the best point, the original sample or point::max
    Foam::reduce
    (
        hitPoints.data(),
        hitPoints.size(),
        minOp<point>(),
        UPstream::msgType(),
        communicator
    );

    if (debug)
    {
        InfoInFunction << nl;

        const auto& samples = probeLocations();

        forAll(nearest, samplei)
        {
            const auto proci = nearest[samplei].second();
            const auto facei = faceIds_[samplei];
            const point& sample = samples[samplei];

            Info<< "    " << samplei << " coord:" << sample
                << " found on processor:"
                << (proci < UPstream::nProcs() ? proci : -1)
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

    // Re-filter the local faceIds_/patchIds_
    // to only retain the globally nearest sample point.

    forAll(nearest, samplei)
    {
        const auto proci = nearest[samplei].second();

        // Valid (local) face to sample?
        if (proci == myProci && faceIds_[samplei] >= 0)
        {
            procIds_[samplei] = proci;
        }
        else
        {
            faceIds_[samplei] = -1;
            patchIds_[samplei] = -1;
            procIds_[samplei] = -1;
        }
    }

    // NB: never reduce (cellIds_, faceIds_) : they are processor-local

    Foam::reduce
    (
        procIds_.data(),
        procIds_.size(),
        maxOp<label>(),
        UPstream::msgType(),
        communicator
    );

    Foam::reduce
    (
        patchIds_.data(),
        patchIds_.size(),
        maxOp<label>(),
        UPstream::msgType(),
        communicator
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
