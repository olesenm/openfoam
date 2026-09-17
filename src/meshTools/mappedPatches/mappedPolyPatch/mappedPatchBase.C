/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearchMeshObject.H"
#include "multiWorldConnectionsObject.H"
#include "meshTools.H"
#include "OFstream.H"
#include "treeDataFace.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "Time.H"
#include "mapDistribute.H"
#include "SubField.H"
#include "triangle.H"
#include "syncTools.H"
#include "treeDataCell.H"
#include "DynamicField.H"
#include "faceAreaWeightAMI.H"
#include "OTstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBase, 0);
}


const Foam::Enum
<
    Foam::mappedPatchBase::sampleMode
>
Foam::mappedPatchBase::sampleModeNames_
({
    { sampleMode::NEARESTCELL, "nearestCell" },
    { sampleMode::NEARESTPATCHFACE, "nearestPatchFace" },
    { sampleMode::NEARESTPATCHFACEAMI, "nearestPatchFaceAMI" },
    { sampleMode::NEARESTPATCHPOINT, "nearestPatchPoint" },
    { sampleMode::NEARESTFACE, "nearestFace" },
    { sampleMode::NEARESTONLYCELL, "nearestOnlyCell" },
});


const Foam::Enum
<
    Foam::mappedPatchBase::offsetMode
>
Foam::mappedPatchBase::offsetModeNames_
({
    { offsetMode::UNIFORM, "uniform" },
    { offsetMode::NONUNIFORM, "nonuniform" },
    { offsetMode::NORMAL, "normal" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mappedPatchBase::updateMesh(PstreamBuffers& pBufs)
{
    clearOut();
}


Foam::autoPtr<Foam::fileName> Foam::mappedPatchBase::readDatabase
(
    const dictionary& dict
)
{
    if (dict.found("sampleDatabase"))
    {
        if (dict.get<bool>("sampleDatabase"))
        {
            return autoPtr<fileName>::New
            (
                dict.getOrDefault<fileName>
                (
                    "sampleDatabasePath",
                    fileName::null
                )
            );
        }
    }
    else if (dict.found("sampleDatabasePath"))
    {
        return autoPtr<fileName>::New(dict.get<fileName>("sampleDatabasePath"));
    }

    return nullptr;
}


bool Foam::mappedPatchBase::addWorldConnection()
{
    if (sameWorld())
    {
        return true;
    }

    const Time& runTime = patch_.boundaryMesh().mesh().time();
    return const_cast<multiWorldConnections&>
    (
        multiWorldConnections::New(runTime)
    ).addConnectionByName(sampleWorld_);
}


Foam::label Foam::mappedPatchBase::getWorldCommunicator() const
{
    if (sameWorld())
    {
        return UPstream::commWorld();
    }

    const Time& runTime = patch_.boundaryMesh().mesh().time();
    return
        multiWorldConnections::New(runTime).getCommByName(sampleWorld_);
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::facePoints
(
    const polyPatch& pp
) const
{
    const polyMesh& mesh = pp.boundaryMesh().mesh();

    // Force construction of min-tet decomp
    (void)mesh.tetBasePtIs();

    // Initialise to face-centre
    auto tfacePoints = tmp<pointField>::New(patch_.size());
    auto& facePoints = tfacePoints.ref();

    forAll(pp, facei)
    {
        facePoints[facei] = facePoint
        (
            mesh,
            pp.start()+facei,
            polyMesh::FACE_DIAG_TRIS
        ).point();
    }

    return tfacePoints;
}


void Foam::mappedPatchBase::collectSamples
(
    const label mySampleWorld,      // Wanted world
    const pointField& facePoints,
    pointField& samples,            // All samples
    labelList& patchFaceWorlds,     // Per sample: wanted world
    labelList& patchFaceProcs,      // Per sample: originating processor
    labelList& patchFaces,          // Per sample: originating patchFace index
    pointField& patchFc             // Per sample: originating centre
) const
{
    DebugInFunction << nl;

    const auto myComm = getCommunicator();  // Get or create
    const auto myRank = UPstream::myProcNo(myComm);
    const auto nProcs = UPstream::nProcs(myComm);

    const auto oldWarnComm = UPstream::commWarn(myComm);

    constexpr auto commsType = UPstream::commsTypes::nonBlocking;
    const int uniqueTag = (UPstream::msgType() + 1492);

    if (debug & 2)
    {
        Perr<< "patch: " << patch_.name()
            << "[rank=" << myRank << " procs=" << nProcs
            << " comm=" << myComm << "] collect samples" << endl;
    }

    // Collect all sample points and the faces they come from.

    // Older code (v2606 and earlier) used
    // allGatherList and ListListOps::combine, but that gives
    // serialization problems (buffer overrun) with very large cases.
    //
    // Instead use combination globalIndex and broadcast which
    // enables send/recv/broadcast of contiguous data, albeit with
    // a data bottleneck on the master rank.

    do
    {
        const label nLocalPoints = facePoints.size();

        const globalIndex addr
        (
            globalIndex::gatherOnly{},
            nLocalPoints,
            myComm
        );

        label totalSize = addr.totalSize();
        Pstream::broadcast(totalSize, myComm);

        if (!totalSize)
        {
            patchFc.clear();
            samples.clear();
            break;
        }

        // Concatenate facePoints into a straight list
        addr.gather(facePoints, patchFc, uniqueTag, commsType, myComm);
        if (UPstream::is_subrank(myComm))
        {
            patchFc.resize_nocopy(totalSize);
        }
        UPstream::broadcast(patchFc.data(), patchFc.size(), myComm);

        auto tsampledValues = samplePoints(facePoints);

        // Concatenate sampled values into a straight list
        addr.gather(tsampledValues(), samples, uniqueTag, commsType, myComm);
        tsampledValues.clear();

        if (UPstream::is_subrank(myComm))
        {
            samples.resize_nocopy(totalSize);
        }
        UPstream::broadcast(samples.data(), samples.size(), myComm);

    } while (false);


    // Now handle world and processor addressing
    {
        const labelList procToWorldIndex
        (
            UPstream::allGatherValues<label>(mySampleWorld, myComm)
        );
        const labelList nFacesPerProc
        (
            UPstream::allGatherValues<label>(patch_.size(), myComm)
        );

        label nTotalFaces = 0;
        for (auto n : nFacesPerProc)
        {
            nTotalFaces += n;
        }

        patchFaceWorlds.resize_nocopy(nTotalFaces);
        patchFaceProcs.resize_nocopy(nTotalFaces);
        patchFaces.resize_nocopy(nTotalFaces);

        label start = 0;

        forAll(nFacesPerProc, proci)
        {
            const auto len = nFacesPerProc[proci];
            const auto worldIndex = procToWorldIndex[proci];

            patchFaceProcs.slice(start, len) = proci;
            patchFaceWorlds.slice(start, len) = worldIndex;

            // Local one-to-one identity mapping for patchFaces
            {
                auto slice = patchFaces.slice(start, len);
                Foam::identity(slice);
            }

            start += len;
        }
    }

    // Restore communicator settings
    UPstream::commWarn(oldWarnComm);
}


Foam::List<Foam::mappedPatchBase::nearKey>
Foam::mappedPatchBase::findLocalSamples
(
    const sampleMode mode,

    const word& sampleRegion,
    const word& samplePatch,

    const pointField& samples,

    // Indices within the locally sampled entity
    labelList& sampleIndices,
    // Locations within the locally sampled entity
    pointField& sampleLocations,
    // optional: subset filter (zero-size == no subsetting)
    const bitSet& sampleFilter
) const
{
    DebugInFunction << nl;

    const auto myComm = getCommunicator();  // Get or create
    const auto myRank = UPstream::myProcNo(myComm);

    // Lookup the correct region
    const polyMesh& mesh = lookupMesh(sampleRegion);

    // Initialize as out-of-range
    sampleIndices.resize_fill(samples.size(), -1);
    sampleLocations.resize_fill(samples.size(), point::max);

    // The nearest selection key. Construct as missed sample
    // == (large distance, proc = -1)
    List<nearKey> nearest(samples.size(), nearKey(Foam::sqr(GREAT), -1));

    // For multi-world, accept or ignore this sample?
    const auto accept = [&sampleFilter](label i) noexcept
    {
        return sampleFilter.empty_or_test(i);
    };

    // Update the indices and locations.
    // Caller already checked that pointIndexHit is hit()
    const auto setHit =
        [&sampleIndices,&sampleLocations,&nearest,myRank]
        (
            label i,
            const pointIndexHit& pHit,
            scalar distSqr
        )
        {
            sampleIndices[i] = pHit.index();
            sampleLocations[i] = pHit.point();
            nearest[i].first() = distSqr;
            nearest[i].second() = myRank;
        };


    switch (mode)
    {
        case NEARESTCELL:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const auto& tree = mesh.cellTree();

            forAll(samples, samplei)
            {
                if (!accept(samplei)) continue;

                const point& sample = samples[samplei];

                if (label celli = tree.findInside(sample); celli > 0)
                {
                    const point& cc = mesh.cellCentres()[celli];

                    pointIndexHit pHit(true, cc, celli);

                    setHit(samplei, pHit, sample.distSqr(cc));
                }
            }
            break;
        }

        case NEARESTONLYCELL:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const auto& tree = mesh.cellTree();

            forAll(samples, samplei)
            {
                if (!accept(samplei)) continue;

                const point& sample = samples[samplei];

                pointIndexHit pHit = tree.findNearest(sample, Foam::sqr(GREAT));

                if (const auto& p = pHit.point(); pHit.hit())
                {
                    setHit(samplei, pHit, sample.distSqr(p));
                }
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            constexpr label rand_seed = 123456;

            const polyPatch& pp = lookupPatch(sampleRegion, samplePatch);

            // NB: an empty patch == "missed" sampling (already handled above)

            if (!pp.empty())
            {
                treeBoundBox patchBb(pp.box());

                // Inflate by 0.1% and clamp with min size
                patchBb.inflate_random(rand_seed, 1e-3, ROOTVSMALL);

                const indexedOctree<treeDataFace> boundaryTree
                (
                    treeDataFace(pp),  // Patch faces

                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );
                const auto& treeData = boundaryTree.shapes();

                // Search limit for findNearest
                const scalar maxDistSqr = patchBb.magSqr();

                forAll(samples, samplei)
                {
                    if (!accept(samplei)) continue;

                    const point& sample = samples[samplei];

                    pointIndexHit pHit =
                        boundaryTree.findNearest(sample, maxDistSqr);

                    if (pHit.hit())
                    {
                        const point& p = treeData.centre(pHit.index());

                        pHit.setPoint(p);
                        setHit(samplei, pHit, sample.distSqr(p));
                    }
                }
            }
            break;
        }

        case NEARESTPATCHPOINT:
        {
            constexpr label rand_seed = 123456;

            const polyPatch& pp = lookupPatch(sampleRegion, samplePatch);

            // NB: an empty patch == "missed" sampling (already handled above)

            if (!pp.empty())
            {
                // patch (local) points
                treeBoundBox patchBb(pp.points(), pp.meshPoints());

                // Inflate by 0.1% and clamp with min size
                patchBb.inflate_random(rand_seed, 1e-3, ROOTVSMALL);

                const indexedOctree<treeDataPoint> boundaryTree
                (
                    // Patch points
                    treeDataPoint(mesh.points(), pp.meshPoints()),

                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );

                // Search limit for findNearest
                const scalar maxDistSqr = patchBb.magSqr();

                forAll(samples, samplei)
                {
                    if (!accept(samplei)) continue;

                    const point& sample = samples[samplei];

                    pointIndexHit pHit =
                        boundaryTree.findNearest(sample, maxDistSqr);

                    if (const auto& p = pHit.point(); pHit.hit())
                    {
                        setHit(samplei, pHit, sample.distSqr(p));
                    }
                }
            }
            break;
        }

        case NEARESTFACE:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const auto& meshSearchEngine = meshSearchMeshObject::New(mesh);

            forAll(samples, samplei)
            {
                if (!accept(samplei)) continue;

                const point& sample = samples[samplei];

                if
                (
                    label facei = meshSearchEngine.findNearestFace(sample);
                    (facei >= 0)
                )
                {
                    const point& p = mesh.faceCentres()[facei];

                    pointIndexHit pHit(true, p, facei);

                    setHit(samplei, pHit, sample.distSqr(p));
                }
            }
            break;
        }

        case NEARESTPATCHFACEAMI:
        {
            // nothing to do here, treat as "missed" sampling
            // (already handled above)
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration ("
                << int(mode) << ") = programming problem!"
                << abort(FatalError);
            break;
        }
    }

    return nearest;
}


void Foam::mappedPatchBase::findSamples
(
    const sampleMode mode,
    const label myWorld,
    const pointField& samples,
    const labelUList& wantedWorlds,
    const labelUList& origProcs,

    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    DebugInFunction << nl;

    // Find the processor/cell containing the samples. Does not account
    // for samples being found in two processors.

    const auto myComm = getCommunicator();  // Get or create
    const auto myRank = UPstream::myProcNo(myComm);
    const auto nProcs = UPstream::nProcs(myComm);

    const int uniqueTag = (UPstream::msgType() + 1492);

    const auto oldWarnComm = UPstream::commWarn(myComm);

    // The final output are flat lists, mark all as invalid
    sampleProcs.resize_fill(samples.size(), -1);
    sampleIndices.resize_fill(samples.size(), -1);
    sampleLocations.resize_fill(samples.size(), point::max);


    // Tuple of patch/region used for sampling:
    List<wordPair> sampledObjects(nProcs);
    {
        sampledObjects[myRank].first() = samplePatch_;
        sampledObjects[myRank].second() = sampleRegion_;
        Pstream::allGatherList(sampledObjects, uniqueTag, myComm);
    }


    // Subset of worlds to sample. Empty == no subset needed

    bitSet worldSelection(samples.size());

    forAll(wantedWorlds, i)
    {
        if (wantedWorlds[i] == myWorld)
        {
            worldSelection.set(i);
        }
    }

    const auto leadProci = worldSelection.find_first();

    // If sampling is a single world only, don't bother with the extra filter.
    // This will be most of the cases.
    if (worldSelection.all())
    {
        worldSelection.clearStorage();
    }

    // Assume single patch to sample for now
    const word origPatch =
    (
        (leadProci >= 0)
      ? sampledObjects[origProcs[leadProci]].first()
      : word::null
    );

    const word origRegion =
    (
        (leadProci >= 0)
      ? sampledObjects[origProcs[leadProci]].second()
      : word::null
    );

    if (debug)
    {
        const bool localOnly = worldSelection.empty();

        label nLocalOnly = (localOnly ? 1 : 0);
        Foam::reduce(nLocalOnly, sumOp<label>(), UPstream::msgType(), myComm);

        Info<< "Searching locally (single world) on " << nLocalOnly
            << '/' << nProcs << " ranks" << endl;

        if (localOnly)
        {
            Pout<< "Searching locally (single world) for " << samples.size()
                << " samples on region:" << origRegion
                << " on patch:" << origPatch << endl;
        }
        else if (worldSelection.none())
        {
            Pout<< "No local search (world:" << myWorld << ")" << endl;
        }
        else
        {
            Pout<< "Searching locally (world:" << myWorld
                << " ) for " << worldSelection.count()
                << " samples on region:" << origRegion
                << " on patch:" << origPatch << endl;
        }
    }


    // The selection keys for MINLOC-style reduction.
    // (distance, rank)

    List<nearKey> selectionKeys
    (
        findLocalSamples
        (
            mode,
            origRegion,
            origPatch,
            samples,

            sampleIndices,
            sampleLocations,

            worldSelection
        )
    );

    if (debug)
    {
        PrintMemoryInFunction << endl;
    }

    // Reduction: Stage 1
    //
    // reduce the (distance, rank) key as a MINLOC type of operation

    Pstream::listCombineReduce
    (
        selectionKeys,
        minFirstEqOp{},
        uniqueTag,
        myComm
    );


    // Reduction: Stage 2
    //
    // use (globally consistent) selection key information to decide
    // on the originating rank. Use this to markup the corresponding
    // entries in sample procs, indices, locations.
    //
    // With (-1, -1, point::max) denoting an invalid entry,
    // can simply use max() for indices/ranks and min() for locations.

    forAll(selectionKeys, i)
    {
        const auto& rank = selectionKeys[i].second();

        if (rank == myRank)
        {
            // My rank provides values.
            // Retain sample indices and location
            sampleProcs[i] = myRank;
        }
        else
        {
            // My rank provides nothing
            sampleIndices[i] = -1;
            sampleProcs[i] = -1;
            sampleLocations[i] = point::max;
        }
    }

    // Reduction: Stage 3
    //
    // let MPI_Allreduce work everything out

    // Selects the winning index or -1
    Foam::reduce
    (
        sampleIndices.data(),
        sampleIndices.size(),
        maxOp<label>(),
        UPstream::msgType(),
        myComm
    );

    // Selects the winning rank or -1
    Foam::reduce
    (
        sampleProcs.data(),
        sampleProcs.size(),
        maxOp<label>(),
        UPstream::msgType(),
        myComm
    );

    // Selects the best point or point::max
    Foam::reduce
    (
        sampleLocations.data(),
        sampleLocations.size(),
        minOp<point>(),
        UPstream::msgType(),
        myComm
    );

    if constexpr (false)  // OR: if (debug >= 4)
    {
        Pout<< "** After combining:" << endl;
        forAll(samples, i)
        {
            const auto rank  = selectionKeys[i].second();

            Pout<< "  sample: " << samples[i]
                << " from proc:" << origProcs[i]
                << " (world:" << wantedWorlds[i]
                << " ) found on proc:" << rank
                << " at index:" << sampleIndices[i]
                << " point:" << sampleLocations[i] << nl;
        }
        Pout<< endl;
    }

    // Restore communicator settings
    UPstream::commWarn(oldWarnComm);
}


void Foam::mappedPatchBase::calcMapping() const
{
    DebugInFunction;

    static bool hasWarned = false;
    if (mapPtr_)
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    DebugInFunction << nl;


    // Early construction of tetBasePtIs since does parallel
    // communication which might conflict with inter-world communicator
    // (for multi-world operation)
    (void)patch_.boundaryMesh().mesh().tetBasePtIs();

    const label myComm = getCommunicator();  // Get or create

    //// Make sure if running in database that there is a syncObjects FO
    //if (sampleDatabase() && !sameWorld())
    //{
    //    const word syncName("syncObjects");
    //    const polyMesh& mesh = patch_.boundaryMesh().mesh();
    //    const Time& runTime = mesh.time();
    //    const functionObjectList& fos = runTime.functionObjects();
    //
    //    forAll(fos, i)
    //    {
    //        Pout<< "** FO:" << fos[i].name() << " tpye:" << fos[i].type()
    //            << endl;
    //    }
    //
    //    if (fos.findObjectID(syncName) == -1)
    //    {
    //        //FatalErrorInFunction
    //        WarningInFunction
    //            << "Did not detect functionObject " << syncName
    //            << ". This is used to synchronise data inbetween worlds"
    //            << endl;
    //    }
    //}


    // Get points on face (since cannot use face-centres - might be off
    // face-diagonal decomposed tets.
    tmp<pointField> patchPoints(facePoints(patch_));

    // Get offsetted points
    const pointField offsettedPoints(samplePoints(patchPoints()));

    // Do a sanity check - am I sampling my own patch?
    // This only makes sense for a non-zero offset.
    const bool sampleMyself =
    (
        mode_ == NEARESTPATCHFACE
     && sampleWorld() == UPstream::myWorld()
     && sampleRegion() == patch_.boundaryMesh().mesh().name()
     && samplePatch() == patch_.name()
    );

    if (sampleMyself)
    {
        // Check offset
        vectorField d(offsettedPoints-patchPoints());
        bool coincident = (gAverage(mag(d)) <= ROOTVSMALL);

        if (sampleMyself && coincident)
        {
            WarningInFunction
                << "Invalid offset " << d << endl
                << "Offset is the vector added to the patch face centres to"
                << " find the patch face supplying the data." << endl
                << "Setting it to " << d
                << " on the same patch, on the same region"
                << " will find the faces themselves which does not make sense"
                << " for anything but testing." << endl
                << "patch:" << patch_.name() << endl
                << "sampleRegion:" << sampleRegion() << endl
                << "mode:" << sampleModeNames_[mode_] << endl
                << "samplePatch:" << samplePatch() << endl
                << "offsetMode:" << offsetModeNames_[offsetMode_] << endl;
        }
    }

    // Get local world and world-to-sample in index form
    const label myWorld = Pstream::myWorldID();
    const label mySampleWorld = Pstream::allWorlds().find(sampleWorld_);


    // Get global list of all samples and the processor and face they come from.
    pointField samples;         // coordinates
    labelList patchFaceWorlds;  // world to sample
    labelList patchFaceProcs;   // originating processor
    labelList patchFaces;       // originating face on processor patch
    pointField patchFc;         // originating face centre
    collectSamples
    (
        mySampleWorld,          // world I want to sample
        patchPoints,

        samples,
        patchFaceWorlds,
        patchFaceProcs,
        patchFaces,
        patchFc
    );

    if constexpr (false)  // OR: if (debug >= 4)
    {
        forAll(samples, i)
        {
            Pout<< "    sample:" << samples[i]
                << " from proc:" << patchFaceProcs[i]
                << " face:" << patchFaces[i]
                << " on world:" << patchFaceWorlds[i] << nl;
        }
    }

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    pointField sampleLocations;
    findSamples
    (
        mode_,
        myWorld,        // my world (in index form)
        samples,
        patchFaceWorlds,
        patchFaceProcs,

        sampleProcs,
        sampleIndices,
        sampleLocations
    );

    // Check for samples that were not found. This will only happen for
    // NEARESTCELL since finds cell containing a location
    if (mode_ == NEARESTCELL)
    {
        // Collect the samples that cannot be found.
        // Note that sampleProcs is already parallel consistent, so no futher
        // reduction or coordination is required here.
        DynamicList<label> subMap;
        forAll(sampleProcs, samplei)
        {
            if (sampleProcs[samplei] == -1)
            {
                subMap.push_back(samplei);
            }
        }

        if (label nNotFound = subMap.size(); nNotFound > 0)
        {
            if (!hasWarned)
            {
                hasWarned = true;
                WarningInFunction
                    << "Did not find " << nNotFound
                    << " out of " << sampleProcs.size() << " total samples."
                    << " Sampling these on owner cell centre instead." << nl
                    << "On patch " << patch_.name()
                    << " on region " << sampleRegion()
                    << " in mode " << sampleModeNames_[mode_] << nl
                    << "with offset mode " << offsetModeNames_[offsetMode_]
                    << ". Suppressing further warnings from " << type() << endl;
            }

            const pointField samples_subset(samples, subMap);
            const labelList worlds_subset(patchFaceWorlds, subMap);
            const labelList procs_subset(patchFaceProcs, subMap);

            // And re-search for pure nearest (should not fail)
            labelList subSampleProcs;
            labelList subSampleIndices;
            pointField subSampleLocations;
            findSamples
            (
                NEARESTONLYCELL,
                myWorld,        // my world (in index form)

                samples_subset,
                worlds_subset,
                procs_subset,

                subSampleProcs,
                subSampleIndices,
                subSampleLocations
            );

            // Insert
            UIndirectList<label>(sampleProcs, subMap) = subSampleProcs;
            UIndirectList<label>(sampleIndices, subMap) = subSampleIndices;
            UIndirectList<point>(sampleLocations, subMap) = subSampleLocations;
        }
    }

    // Now we have all the data we need:
    // - where sample originates from (so destination when mapping):
    //   patchFaces, patchFaceProcs.
    // - cell/face sample is in (so source when mapping)
    //   sampleIndices, sampleProcs.

    if (UPstream::master(myComm))
    {
        forAll(samples, i)
        {
            if (sampleProcs[i] == -1)
            {
                FatalErrorInFunction << "did not find sample "
                    << patchFc[i] << " on patch " << patch_.name()
                    << " on region "
                    << patch_.boundaryMesh().mesh().name()
                    << " on processor " << patchFaceProcs[i]
                    << exit(FatalError);
            }
        }
    }


    if (debug && UPstream::master(myComm))
    {
        if constexpr (false)  // OR: if (debug >= 4)
        {
            const auto& meshName = patch_.boundaryMesh().mesh().name();
            const auto& sampleName = sampleRegion();

            forAll(samples, i)
            {
                Pout<< i << " need data in region " << meshName
                    << " for proc:" << patchFaceProcs[i]
                    << " face:" << patchFaces[i]
                    << " at:" << patchFc[i] << endl
                    << "Found data in region " << sampleName
                    << " at proc:" << sampleProcs[i]
                    << " face:" << sampleIndices[i]
                    << " at:" << sampleLocations[i]
                    << nl << endl;
            }
        }

        OFstream str
        (
            patch_.boundaryMesh().mesh().time().path()
          / patch_.name()
          + "_mapped.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cell/faceCentres/points to file " << str.name()
            << endl;

        label vertIndex(0);
        forAll(patchFc, i)
        {
            meshTools::writeOBJ(str, patchFc[i], sampleLocations[i], vertIndex);
        }
    }

    // Determine schedule.
    mapPtr_.reset(new mapDistribute(sampleProcs, patchFaceProcs, myComm));

    // Rework the schedule from indices into samples to cell data to send,
    // face data to receive.

    labelListList& subMap = mapPtr_().subMap();
    labelListList& constructMap = mapPtr_().constructMap();

    forAll(subMap, proci)
    {
        subMap[proci] = labelUIndList(sampleIndices, subMap[proci]);
        constructMap[proci] = labelUIndList(patchFaces, constructMap[proci]);

        if (debug)
        {
            Pout<< "To proc:" << proci << " sending values of cells/faces:"
                << subMap[proci] << endl;
            Pout<< "From proc:" << proci
                << " receiving values of patch faces:"
                << constructMap[proci] << endl;
        }
    }


    // Redo constructSize
    mapPtr_().constructSize() = patch_.size();

    if (debug)
    {
        // Check that all elements get a value.
        bitSet used(patch_.size());
        forAll(constructMap, proci)
        {
            const labelList& map = constructMap[proci];

            forAll(map, i)
            {
                label facei = map[i];

                if (used.test(facei))
                {
                    FatalErrorInFunction
                        << "On patch " << patch_.name()
                        << " patchface " << facei
                        << " is assigned to more than once."
                        << abort(FatalError);
                }
                else
                {
                    used.set(facei);
                }
            }
        }
        forAll(used, facei)
        {
            if (!used.test(facei))
            {
                FatalErrorInFunction
                    << "On patch " << patch_.name()
                    << " patchface " << facei
                    << " is never assigned to."
                    << abort(FatalError);
            }
        }
    }

    updateMeshTime().setUpToDate();
    if (sameWorld())
    {
        updateSampleMeshTime().setUpToDate();
    }
}


const Foam::autoPtr<Foam::searchableSurface>& Foam::mappedPatchBase::surfPtr()
const
{
    const word surfType(surfDict_.getOrDefault<word>("type", "none"));

    if (!surfPtr_ && surfType != "none")
    {
        word surfName(surfDict_.getOrDefault("name", patch_.name()));

        const polyMesh& mesh = patch_.boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


void Foam::mappedPatchBase::calcAMI() const
{
    if (AMIPtr_->upToDate())
    {
        DebugInFunction
            << "AMI already up-to-date"
            << endl;

        return;
    }

    DebugInFunction << nl;

    const label myComm = getCommunicator();  // Get or create

    // Pre-calculate surface (if any)
    const auto& surf = surfPtr();

    // Check if running locally
    if (sampleWorld_.empty() || sameWorld())
    {
        const polyPatch& nbr = samplePolyPatch();

        // Transform neighbour patch to local system
        const pointField nbrPoints(samplePoints(nbr.localPoints()));

        const primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );


        if (debug)
        {
            OFstream os(patch_.name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, samplePolyPatch().localFaces(), nbrPoints);

            OFstream osN(patch_.name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0, nbrPoints);

            OFstream osO(patch_.name() + "_ownerPatch.obj");
            meshTools::writeOBJ(osO, patch_.localFaces(), patch_.localPoints());
        }

        // Construct/apply AMI interpolation to determine addressing and
        // weights.

        // Change to use inter-world communicator
        const auto oldWarnComm = UPstream::warnComm;
        const auto oldWorldComm = UPstream::worldComm;
        if (!sameWorld())
        {
            UPstream::commWarn(myComm);
            UPstream::commWorld(myComm);
        }

        AMIPtr_->calculate(patch_, nbrPatch0, surf);

        // Restore communicator settings
        UPstream::commWarn(oldWarnComm);
        UPstream::commWorld(oldWorldComm);
    }
    else
    {
        faceList dummyFaces;
        pointField dummyPoints;
        const primitivePatch dummyPatch
        (
            SubList<face>(dummyFaces),
            dummyPoints
        );

        // Change to use inter-world communicator
        AMIPtr_->comm(myComm);

        const auto oldWarnComm = UPstream::commWarn(myComm);
        const auto oldWorldComm = UPstream::commWorld(myComm);

        if (masterWorld())
        {
            // Construct/apply AMI interpolation to determine addressing
            // and weights. Have patch_ for src faces, 0 faces for the
            // target side
            AMIPtr_->calculate(patch_, dummyPatch, surf);
        }
        else
        {
            // Construct/apply AMI interpolation to determine addressing
            // and weights. Have 0 faces for src side, patch_ for the tgt
            // side
            AMIPtr_->calculate(dummyPatch, patch_, surf);
        }
        // Now the AMI addressing/weights will be from src side (on masterWorld
        // processors) to tgt side (on other processors)

        // Restore communicator settings
        UPstream::commWarn(oldWarnComm);
        UPstream::commWorld(oldWorldComm);
    }
}


const Foam::objectRegistry& Foam::mappedPatchBase::subRegistry
(
    const objectRegistry& obr,
    const wordList& names,
    const label index
)
{
    const objectRegistry& sub = obr.subRegistry(names[index], true);
    if (index == names.size()-1)
    {
        return sub;
    }
    else
    {
        return subRegistry(sub, names, index+1);
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::mappedPatchBase(const polyPatch& pp)
:
    patch_(pp),
    sampleWorld_(),
    sampleRegion_(patch_.boundaryMesh().mesh().name()),
    mode_(NEARESTPATCHFACE),
    samplePatch_(),
    coupleGroup_(),
    sampleDatabasePtr_(),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(pp.size(), offset_),
    distance_(0),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_(true),
    mapPtr_(nullptr),
    AMIReverse_(false),
    AMIPtr_(new faceAreaWeightAMI(true, AMIReverse_)),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{
    // NOTE: same region, no sample-world. Thus no world-world communication
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vectorField& offsets
)
:
    mappedPatchBase
    (
        pp,
        sampleRegion,
        mode,
        samplePatch,
        scalar(0)
    )
{
    mappedPatchBase::setOffset(offsets);
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vector& uniformOffset
)
:
    mappedPatchBase
    (
        pp,
        sampleRegion,
        mode,
        samplePatch,
        scalar(0)
    )
{
    mappedPatchBase::setOffset(uniformOffset);
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const scalar normalDistance
)
:
    patch_(pp),
    sampleWorld_(),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    sampleDatabasePtr_(),
    offsetMode_(NORMAL),
    offset_(Zero),
    offsets_(0),
    distance_(normalDistance),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIReverse_(false),
    AMIPtr_(new faceAreaWeightAMI(true, AMIReverse_)),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{
    addWorldConnection();
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    sampleWorld_(dict.getOrDefault<word>("sampleWorld", word::null)),
    sampleRegion_(dict.getOrDefault<word>("sampleRegion", word::null)),
    mode_(sampleModeNames_.get("sampleMode", dict)),
    samplePatch_(dict.getOrDefault<word>("samplePatch", word::null)),
    coupleGroup_(dict),
    sampleDatabasePtr_(readDatabase(dict)),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(0),
    distance_(0),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIReverse_
    (
        dict.getOrDefault
        (
            "reverseTarget",        // AMIInterpolation uses this keyword
            dict.getOrDefault
            (
                "flipNormals",
                false
            )
        )
    ),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.getOrDefault("AMIMethod", faceAreaWeightAMI::typeName),
            dict,
            AMIReverse_
        )
    ),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    addWorldConnection();

    if (!coupleGroup_.good())
    {
        if (sampleWorld_.empty() && sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
            sameRegion_ = true;
        }
    }

    if (offsetModeNames_.readIfPresent("offsetMode", dict, offsetMode_))
    {
        switch (offsetMode_)
        {
            case UNIFORM:
            {
                dict.readEntry("offset", offset_);
            }
            break;

            case NONUNIFORM:
            {
                offsets_ = pointField("offsets", dict, patch_.size());
            }
            break;

            case NORMAL:
            {
                dict.readEntry("distance", distance_);
            }
            break;
        }
    }
    else if (dict.readIfPresent("offset", offset_))
    {
        offsetMode_ = UNIFORM;
    }
    else if (dict.found("offsets"))
    {
        offsetMode_ = NONUNIFORM;
        offsets_ = pointField("offsets", dict, patch_.size());
    }
    else if (mode_ != NEARESTPATCHFACE && mode_ != NEARESTPATCHFACEAMI)
    {
        FatalIOErrorInFunction(dict)
            << "Please supply the offsetMode as one of "
            << offsetModeNames_
            << exit(FatalIOError);
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const sampleMode mode,
    const dictionary& dict
)
:
    patch_(pp),
    sampleWorld_(dict.getOrDefault<word>("sampleWorld", word::null)),
    sampleRegion_(dict.getOrDefault<word>("sampleRegion", word::null)),
    mode_(mode),
    samplePatch_(dict.getOrDefault<word>("samplePatch", word::null)),
    coupleGroup_(dict), //dict.getOrDefault<word>("coupleGroup", word::null)),
    sampleDatabasePtr_(readDatabase(dict)),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(0),
    distance_(0),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIReverse_
    (
        dict.getOrDefault
        (
            "reverseTarget",        // AMIInterpolation uses this keyword
            dict.getOrDefault
            (
                "flipNormals",
                false
            )
        )
    ),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.getOrDefault("AMIMethod", faceAreaWeightAMI::typeName),
            dict,
            AMIReverse_
        )
    ),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    addWorldConnection();

    if (mode != NEARESTPATCHFACE && mode != NEARESTPATCHFACEAMI)
    {
        FatalIOErrorInFunction(dict)
            << "Construct from sampleMode and dictionary only applicable for "
            << " collocated patches in modes "
            << sampleModeNames_[NEARESTPATCHFACE] << ','
            << sampleModeNames_[NEARESTPATCHFACEAMI]
            << exit(FatalIOError);
    }

    if (!coupleGroup_.good())
    {
        if (sampleWorld_.empty() && sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
            sameRegion_ = true;
        }
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb
)
:
    patch_(pp),
    sampleWorld_(mpb.sampleWorld_),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    sampleDatabasePtr_
    (
        mpb.sampleDatabasePtr_
      ? new fileName(mpb.sampleDatabasePtr_())
      : nullptr
    ),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_(mpb.offsets_),
    distance_(mpb.distance_),
    communicator_(mpb.communicator_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    AMIPtr_(mpb.AMIPtr_->clone()),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb,
    const labelUList& mapAddressing
)
:
    patch_(pp),
    sampleWorld_(mpb.sampleWorld_),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    sampleDatabasePtr_
    (
        mpb.sampleDatabasePtr_
      ? new fileName(mpb.sampleDatabasePtr_())
      : nullptr
    ),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_
    (
        offsetMode_ == NONUNIFORM
      ? vectorField(mpb.offsets_, mapAddressing)
      : vectorField()
    ),
    distance_(mpb.distance_),
    communicator_(mpb.communicator_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    AMIPtr_(mpb.AMIPtr_->clone()),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::~mappedPatchBase()
{
    clearOut();
}


void Foam::mappedPatchBase::clearOut()
{
    mapPtr_.reset(nullptr);
    AMIPtr_->upToDate(false);

    //Note: no need to clear out surface since not mesh related
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedPatchBase::setOffset(const scalar normalDist)
{
    clearOut();
    offsetMode_ = offsetMode::NORMAL;
    offset_ = Zero;
    offsets_.clear();
    distance_ = normalDist;
}


void Foam::mappedPatchBase::setOffset(const vector& uniformOffset)
{
    clearOut();
    offsetMode_ = offsetMode::UNIFORM;
    offset_ = uniformOffset;
    offsets_.clear();
    distance_ = Zero;
}


void Foam::mappedPatchBase::setOffset(const vectorField& offsets)
{
    clearOut();
    offsetMode_ = offsetMode::NONUNIFORM;
    offset_ = Zero;
    offsets_ = offsets;
    distance_ = Zero;
}


const Foam::polyMesh& Foam::mappedPatchBase::lookupMesh
(
    const word& sampleRegion
) const
{
    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();
    return
    (
        sampleRegion.empty() || sampleRegion == thisMesh.name()
      ? thisMesh
      : thisMesh.time().lookupObject<polyMesh>(sampleRegion)
    );
}


const Foam::polyPatch& Foam::mappedPatchBase::lookupPatch
(
    const word& sampleRegion,
    const word& samplePatch
) const
{
    const polyMesh& nbrMesh = lookupMesh(sampleRegion);

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch
            << " in region " << sampleRegion_ << endl
            << exit(FatalError);
    }
    return nbrMesh.boundaryMesh()[patchi];
}


const Foam::polyMesh& Foam::mappedPatchBase::sampleMesh() const
{
    if (UPstream::myWorld() != sampleWorld_)
    {
        FatalErrorInFunction
            << "sampleWorld : " << sampleWorld_
            << " is not the current world : " << UPstream::myWorld()
            << exit(FatalError);
    }
    return lookupMesh(sampleRegion());
}


const Foam::polyPatch& Foam::mappedPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch()
            << " in region " << sampleRegion_ << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchi];
}


bool Foam::mappedPatchBase::upToDate() const
{
    const auto myComm = getCommunicator();  // Get or create

    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();

    bool thisUpToDate = thisMesh.upToDatePoints(updateMeshTime());
    bool sampleUpToDate =
    (
        sameWorld()
      ? sampleMesh().upToDatePoints(updateSampleMeshTime())
      : true
    );

    // Ensure that states are globally consistent
    if (thisMesh.moving())
    {
        UPstream::reduceAnd(thisUpToDate, myComm);
    }
    if (sameWorld() && sampleMesh().moving() && isPatchSampling())
    {
        UPstream::reduceAnd(sampleUpToDate, myComm);
    }

    // Lambda to check for points on the patch being the same
    // based on the distance change being within a tolerance value

    const scalar maxDistSqr = Foam::sqr(SMALL);

    auto checkPointMovement = [maxDistSqr,myComm]
    (
        const polyMesh& mesh,
        const polyPatch& patch,
        regIOobject& state
    ) -> bool
    {
        bool upToDate = true;
        const auto& oldPoints = mesh.oldPoints();
        const auto& points = mesh.points();

        for (const label pointi : patch.meshPoints())
        {
            if (maxDistSqr < oldPoints[pointi].distSqr(points[pointi]))
            {
                upToDate = false;
                break;
            }
        }

        UPstream::reduceAnd(upToDate, myComm);

        if (upToDate)
        {
            state.setUpToDate();
        }

        return upToDate;
    };


    if (!thisUpToDate && thisMesh.moving())
    {
        // Moving (but not topoChanging mesh) : do more accurate check:
        // compare actual patch point position

        thisUpToDate = checkPointMovement(thisMesh, patch_, updateMeshTime());
    }

    if (!sampleUpToDate && sampleMesh().moving() && isPatchSampling())
    {
        // implicitly (sameWorld() == true) - see above
        sampleUpToDate = checkPointMovement
        (
            sampleMesh(),
            samplePolyPatch(),
            updateSampleMeshTime()
        );
    }

    return (thisUpToDate && sampleUpToDate);
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints
(
    const pointField& fc
) const
{
    auto tfld = tmp<pointField>::New(fc);
    auto& fld = tfld.ref();

    switch (offsetMode_)
    {
        case UNIFORM:
        {
            fld += offset_;
            break;
        }
        case NONUNIFORM:
        {
            fld += offsets_;
            break;
        }
        case NORMAL:
        {
            // Get outwards pointing normal
            vectorField n(patch_.faceAreas());
            n /= mag(n);

            fld += distance_*n;
            break;
        }
    }

    return tfld;
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints() const
{
    return samplePoints(facePoints(patch_));
}


Foam::pointIndexHit Foam::mappedPatchBase::facePoint
(
    const polyMesh& mesh,
    const label facei,
    const polyMesh::cellDecomposition decompMode
)
{
    const point& fc = mesh.faceCentres()[facei];

    switch (decompMode)
    {
        case polyMesh::FACE_PLANES:
        case polyMesh::FACE_CENTRE_TRIS:
        {
            // For both decompositions the face centre is guaranteed to be
            // on the face
            return pointIndexHit(true, fc, facei);
        }
        break;

        case polyMesh::FACE_DIAG_TRIS:
        case polyMesh::CELL_TETS:
        {
            // Find the intersection of a ray from face centre to cell centre
            // Find intersection of (face-centre-decomposition) centre to
            // cell-centre with face-diagonal-decomposition triangles.

            const pointField& p = mesh.points();
            const face& f = mesh.faces()[facei];

            if (f.size() <= 3)
            {
                // Return centre of triangle.
                return pointIndexHit(true, fc, 0);
            }

            const label celli = mesh.faceOwner()[facei];
            const point& cc = mesh.cellCentres()[celli];
            vector d = fc-cc;

            const label fp0 = mesh.tetBasePtIs()[facei];
            const point& basePoint = p[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = p[f[fp]];
                label nextFp = f.fcIndex(fp);
                const point& nextPoint = p[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);
                pointHit hitInfo = tri.intersection
                (
                    cc,
                    d,
                    intersection::HALF_RAY
                );

                if (hitInfo.hit() && hitInfo.distance() > 0)
                {
                    return pointIndexHit(true, hitInfo.point(), i-2);
                }

                fp = nextFp;
            }

            // Fall-back
            return pointIndexHit(false, fc, -1);
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "problem" << abort(FatalError);
            return pointIndexHit();
        }
    }
}


const Foam::objectRegistry& Foam::mappedPatchBase::subRegistry
(
    const objectRegistry& obr,
    const fileName& path
)
{
    // Lookup (and create if non-existing) a registry using
    // '/' separated path. Like 'mkdir -p'

    fileName cleanedPath(path);
    cleanedPath.clean();  // Remove unneeded ".."
    const wordList names(cleanedPath.components());

    if (names.empty())
    {
        return obr;
    }
    else
    {
        return subRegistry(obr, names, 0);
    }
}


Foam::fileName Foam::mappedPatchBase::sendPath
(
    const fileName& root,
    const label proci
)
{
    const word processorName("processor"+Foam::name(proci));
    return root/"send"/processorName;
}


Foam::fileName Foam::mappedPatchBase::sendPath(const label proci) const
{
    return sendPath(sampleDatabasePath(), proci);
}


Foam::fileName Foam::mappedPatchBase::receivePath
(
    const fileName& root,
    const label proci
)
{
    const word processorName("processor"+Foam::name(proci));
    return root/"receive"/processorName;
}


Foam::fileName Foam::mappedPatchBase::receivePath(const label proci) const
{
    return receivePath(sampleDatabasePath(), proci);
}


void Foam::mappedPatchBase::writeDict
(
    const objectRegistry& obr,
    dictionary& dict
)
{
    forAllIters(obr, iter)
    {
        regIOobject* objPtr = iter.val();
        const regIOobject& obj = *objPtr;

        if (isA<objectRegistry>(obj))
        {
            dictionary& d = dict.subDictOrAdd(obj.name());
            writeDict(dynamic_cast<const objectRegistry&>(obj), d);
        }
        else if
        (
            writeIOField<scalar>(obj, dict)
         || writeIOField<vector>(obj, dict)
         || writeIOField<sphericalTensor>(obj, dict)
         || writeIOField<symmTensor>(obj, dict)
         || writeIOField<tensor>(obj, dict)
        )
        {
            // IOField specialisation
        }
        else
        {
            // Not tested. No way of retrieving data anyway ...
            OTstream os;
            obj.writeData(os);

            primitiveEntry* pePtr = new primitiveEntry(obj.name(), os.tokens());
            dict.add(pePtr);
        }
    }
}


void Foam::mappedPatchBase::readDict(const dictionary& d, objectRegistry& obr)
{
    // Reverse of writeDict - reads fields from dictionary into objectRegistry
    for (const entry& e : d)
    {
        if (e.isDict())
        {
            // Add sub registry
            objectRegistry& sub = const_cast<objectRegistry&>
            (
                obr.subRegistry(e.keyword(), true)
            );

            readDict(e.dict(), sub);
        }
        else
        {
            ITstream& is = e.stream();
            token tok(is);

            if
            (
                constructIOField<scalar>(e.keyword(), tok, is, obr)
             || constructIOField<vector>(e.keyword(), tok, is, obr)
             || constructIOField<sphericalTensor>(e.keyword(), tok, is, obr)
             || constructIOField<symmTensor>(e.keyword(), tok, is, obr)
             || constructIOField<tensor>(e.keyword(), tok, is, obr)
            )
            {
                // Done storing field
            }
            else
            {
                FatalErrorInFunction << "Unsupported type " << e.keyword()
                    << exit(FatalError);
            }
        }
    }
}


void Foam::mappedPatchBase::write(Ostream& os) const
{
    os.writeEntry("sampleMode", sampleModeNames_[mode_]);
    os.writeEntryIfDifferent<word>("sampleWorld", word::null, sampleWorld_);
    os.writeEntryIfDifferent<word>("sampleRegion", word::null, sampleRegion_);
    os.writeEntryIfDifferent<word>("samplePatch", word::null, samplePatch_);

    if (sampleDatabasePtr_)
    {
        os.writeEntry("sampleDatabase", Switch::name(true));
        // Write database path if differing
        os.writeEntryIfDifferent<fileName>
        (
            "sampleDatabasePath",
            fileName::null,
            sampleDatabasePtr_()
        );
    }
    coupleGroup_.write(os);

    if
    (
        offsetMode_ == UNIFORM
     && offset_ == vector::zero
     && (mode_ == NEARESTPATCHFACE || mode_ == NEARESTPATCHFACEAMI)
    )
    {
        // Collocated mode. No need to write offset data
    }
    else
    {
        os.writeEntry("offsetMode", offsetModeNames_[offsetMode_]);

        switch (offsetMode_)
        {
            case UNIFORM:
            {
                os.writeEntry("offset", offset_);
                break;
            }
            case NONUNIFORM:
            {
                offsets_.writeEntry("offsets", os);
                break;
            }
            case NORMAL:
            {
                os.writeEntry("distance", distance_);
                break;
            }
        }
    }

    if (mode_ == NEARESTPATCHFACEAMI)
    {
        if (AMIPtr_)
        {
            // Use AMI to write itself. Problem: outputs:
            // - restartUncoveredSourceFace
            // - reverseTarget (instead of flipNormals)
            AMIPtr_->write(os);
        }
        if (!surfDict_.empty())
        {
            surfDict_.writeEntry(surfDict_.dictName(), os);
        }
    }
}


// ************************************************************************* //
