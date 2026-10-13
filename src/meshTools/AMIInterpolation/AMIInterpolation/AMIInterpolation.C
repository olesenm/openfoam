/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

#include "AMIInterpolation.H"
#include "meshTools.H"
#include "mapDistribute.H"
#include "flipOp.H"
#include "profiling.H"
#include "triangle.H"
#include "OFstream.H"
#include "registerSwitch.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AMIInterpolation, 0);
    defineRunTimeSelectionTable(AMIInterpolation, dict);
    defineRunTimeSelectionTable(AMIInterpolation, component);
}

bool Foam::AMIInterpolation::cacheIntersections_ = false;

int Foam::AMIInterpolation::useLocalComm_
(
    debug::optimisationSwitch("localAMIComm", 1)
);
registerOptSwitch
(
    "localAMIComm",
    int,
    Foam::AMIInterpolation::useLocalComm_
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::AMIInterpolation::addToCache(const point& refPt)
{
    DebugInfo<< "-- addToCache" << endl;

    cache_.addToCache(*this, refPt);
}


bool Foam::AMIInterpolation::restoreCache
(
    const point& refPt,
    const Time& runTime
)
{
    DebugInfo<< "-- restoreCache" << endl;

    upToDate_ = cache_.restoreCache(refPt, runTime);

    return upToDate_;
}


Foam::autoPtr<Foam::indexedOctree<Foam::AMIInterpolation::treeType>>
Foam::AMIInterpolation::createTree(const primitivePatch &patch) const
{
    treeBoundBox bb(patch.points(), patch.meshPoints());
    bb.inflate(0.01);

    return autoPtr<indexedOctree<treeType>>::New
    (
        treeType
        (
            false,
            patch,
            indexedOctree<treeType>::perturbTol()
        ),
        bb,                         // overall search domain
        8,                          // maxLevel
        10,                         // leaf size
        3.0                         // duplicity
    );
}


Foam::label Foam::AMIInterpolation::calcDistribution
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const label comm,
    autoPtr<UPstream::communicator>& geomComm
) const
{
    // Either not parallel or no faces on any processor
    label proci = 0;

    if (UPstream::parRun())
    {
        // Involved in communication pattern?
        bool inCommGroup = (srcPatch.size() > 0 || tgtPatch.size() > 0);

        // Track which procs are involved
        const List<bool> hasFaces
        (
            UPstream::allGatherValues<bool>(inCommGroup, comm)
        );

        // Always include master (0) in comm-group?
        // - so messages come from master
        if (useLocalComm_ > 1 && UPstream::master(comm))
        {
            inCommGroup = true;
        }

        // Number of communicating procs (ie, they have local faces)
        label nCommProcs(0);

        // First proc with local faces (when nCommProcs == 1)
        label whichProci(-1);

        // Could use UPstream::splitCommunicator(), but that also incurs
        // global communication to determine who belongs to the same set.
        // Instead, an Allgather of the members and compare with the
        // existing communicator to decide if a new communicator is
        // required.

        DynamicList<label> subProcs(hasFaces.size());
        forAll(hasFaces, i)
        {
            if (hasFaces.test(i))
            {
                whichProci = i;
                ++nCommProcs;
                subProcs.push_back(i);
            }
            else if
            (
                inCommGroup
             && (useLocalComm_ > 1) && (i == UPstream::masterNo())
            )
            {
                // Also include master (0) in comm-group?
                // - so messages come from master
                subProcs.push_back(UPstream::masterNo());
            }
        }

        //
        // Define the AMI communication style
        //

        if (nCommProcs == 0)
        {
            // Probably does not happen. No AMI faces? => no communicator
            geomComm.reset();
        }
        else if (nCommProcs == 1)
        {
            proci = whichProci;

            // No local communicator needed
            geomComm.reset();

            DebugInFunction
                << "AMI local to processor" << proci << endl;
        }
        else  // (nCommProcs > 1)
        {
            proci = -1;

            const label currComm = (geomComm.good() ? geomComm().comm() : -1);

            if (useLocalComm_ == 0)
            {
                // Backwards compatible : no local communicator
                geomComm.reset();
            }
            else if (nCommProcs == UPstream::nProcs(comm))
            {
                // Everyone is involved : no local communicator
                geomComm.reset();
            }
            else if (inCommGroup)
            {
                if (UPstream::sameProcs(currComm, subProcs))
                {
                    // Keep geomComm
                    if (debug)
                    {
                        Pout<< "Retained geomComm:" << currComm
                            << " with " << subProcs.size()
                            << " processors out of " << UPstream::nProcs(comm)
                            << endl;
                    }
                }
                else
                {
                    geomComm.reset(new UPstream::communicator(comm, subProcs));
                    if (debug)
                    {
                        Pout<< "Allocated geomComm:" << geomComm().comm()
                            << " from " << subProcs.size()
                            << " processors out of " << UPstream::nProcs(comm)
                            << endl;
                    }
                }
            }
            else
            {
                // Not inCommGroup, but with local communicator elsewhere
                geomComm.reset(new UPstream::communicator());
                if (debug & 2)
                {
                    Pout<< "Allocated dummy geomComm:" << geomComm().comm()
                        << " src-size:" << srcPatch.size()
                        << " tgt-size:" << tgtPatch.size() << endl;
                }
            }
        }

        DebugInFunction
            << "AMI split across multiple processors "
            << flatOutput(subProcs) << endl;
    }

    return proci;
}


void Foam::AMIInterpolation::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    addProfiling(ami, "AMIInterpolation::projectPointsToSurface");

    DebugInfo<< "AMI: projecting points to surface" << endl;

    List<pointIndexHit> nearInfo;

    surf.findNearest(pts, scalarField(pts.size(), GREAT), nearInfo);

    label nMiss = 0;
    forAll(nearInfo, i)
    {
        const pointIndexHit& pi = nearInfo[i];

        if (pi.hit())
        {
            pts[i] = pi.point();
        }
        else
        {
            // Point remains unchanged
            ++nMiss;
        }
    }

    if (nMiss > 0)
    {
        FatalErrorInFunction
            << "Error projecting points to surface: "
            << nMiss << " faces could not be determined"
            << abort(FatalError);
    }
}


void Foam::AMIInterpolation::normaliseWeights
(
    const scalarList& patchAreas,
    const word& patchName,
    const labelListList& addr,
    scalarListList& wght,
    scalarField& wghtSum,
    const bool conformal,
    const bool output,
    const scalar lowWeightTol,
    const label comm
)
{
    addProfiling(ami, "AMIInterpolation::normaliseWeights");

    // Normalise the weights
    wghtSum.resize_nocopy(wght.size());
    label nLowWeight = 0;

    forAll(wght, facei)
    {
        scalarList& w = wght[facei];

        if (w.size())
        {
            scalar denom = patchAreas[facei];

            scalar s = sum(w);
            scalar t = s/denom;
            if (conformal)
            {
                denom = s;
            }

            forAll(w, i)
            {
                w[i] /= denom;
            }

            wghtSum[facei] = t;
            if (t < lowWeightTol)
            {
                ++nLowWeight;
            }
        }
        else
        {
            wghtSum[facei] = 0;
        }
    }

    if (output && comm != -1 && returnReduceOr(wght.size(), comm))
    {
        auto limits = gMinMax(wghtSum, comm);
        auto avg = gAverage(wghtSum, comm);

        label nLow =
            returnReduce(nLowWeight, sumOp<label>(), UPstream::msgType(), comm);

        Info.masterStream(comm)
            << indent
            << "AMI: Patch " << patchName
            << " sum(weights)"
            << " min:" << limits.min()
            << " max:" << limits.max()
            << " average:" << avg << nl;

        if (nLow)
        {
            Info.masterStream(comm)
                << indent
                << "AMI: Patch " << patchName
                << " identified " << nLow
                << " faces with weights less than " << lowWeightTol
                << endl;
        }
    }
}


void Foam::AMIInterpolation::agglomerate
(
    const autoPtr<mapDistribute>& targetMapPtr,
    const scalarList& fineSrcMagSf,
    const labelListList& fineSrcAddress,
    const scalarListList& fineSrcWeights,

    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing,

    scalarList& srcMagSf,
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarField& srcWeightsSum,
    autoPtr<mapDistribute>& tgtMap,
    const label comm
)
{
    addProfiling(ami, "AMIInterpolation::agglomerate");

    const label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    const label targetCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    // Agglomerate face areas
    {
        //srcMagSf.setSize(sourceRestrictAddressing.size(), 0.0);
        srcMagSf.setSize(sourceCoarseSize, 0.0);

        forAll(sourceRestrictAddressing, facei)
        {
            label coarseFacei = sourceRestrictAddressing[facei];
            srcMagSf[coarseFacei] += fineSrcMagSf[facei];
        }
    }

    // Agglomerate weights and indices
    if (targetMapPtr)
    {
        // We are involved in the communicator but our maps are still empty.
        // Fix 'm up so they are the same size as the communicator.
        const mapDistribute& map = *targetMapPtr;

        if (map.constructMap().empty())
        {
            auto& cMap = const_cast<labelListList&>(map.constructMap());
            cMap.resize_nocopy(UPstream::nProcs(map.comm()));
        }
        if (map.subMap().empty())
        {
            auto& cMap = const_cast<labelListList&>(map.subMap());
            cMap.resize_nocopy(UPstream::nProcs(map.comm()));
        }


        // Get all restriction addressing.
        labelList allRestrict(targetRestrictAddressing);
        map.distribute(allRestrict);

        // So now we have agglomeration of the target side in
        // allRestrict:
        //  0..size-1 : local agglomeration (= targetRestrictAddressing
        //              (but potentially permutated))
        //  size..    : agglomeration data from other processors


        // The trickiness in this algorithm is finding out the compaction
        // of the remote data (i.e. allocation of the coarse 'slots'). We could
        // either send across the slot compaction maps or just make sure
        // that we allocate the slots in exactly the same order on both sending
        // and receiving side (e.g. if the submap is set up to send 4 items,
        // the constructMap is also set up to receive 4 items.


        // Short note about the various types of indices:
        // - face indices : indices into the geometry.
        // - coarse face indices : how the faces get agglomerated
        // - transferred data : how mapDistribute sends/receives data,
        // - slots : indices into data after distribution (e.g. stencil,
        //           srcAddress/tgtAddress). Note: for fully local addressing
        //           the slots are equal to face indices.
        // A mapDistribute has:
        // - a subMap : these are face indices
        // - a constructMap : these are from 'transferred-data' to slots

        labelListList tgtSubMap(Pstream::nProcs(comm));

        // Local subMap is just identity
        {
            tgtSubMap[Pstream::myProcNo(comm)] = identity(targetCoarseSize);
        }

        forAll(map.subMap(), proci)
        {
            if (proci != Pstream::myProcNo(comm))
            {
                // Combine entries that point to the same coarse element.
                // The important bit is to loop over the data (and hand out
                // compact indices ) in 'transferred data' order. This
                // guarantees that we're doing exactly the
                // same on sending and receiving side - e.g. the fourth element
                // in the subMap is the fourth element received in the
                // constructMap

                const labelList& elems = map.subMap()[proci];
                const labelList& elemsMap =
                    map.constructMap()[Pstream::myProcNo(comm)];
                labelList& newSubMap = tgtSubMap[proci];
                newSubMap.resize_nocopy(elems.size());

                labelList oldToNew(targetCoarseSize, -1);
                label newi = 0;

                for (const label elemi : elems)
                {
                    label fineElem = elemsMap[elemi];
                    label coarseElem = allRestrict[fineElem];
                    if (oldToNew[coarseElem] == -1)
                    {
                        oldToNew[coarseElem] = newi;
                        newSubMap[newi] = coarseElem;
                        ++newi;
                    }
                }
                newSubMap.resize(newi);
            }
        }

        // Reconstruct constructMap by combining entries. Note that order
        // of handing out indices should be the same as loop above to compact
        // the sending map

        labelListList tgtConstructMap(Pstream::nProcs(comm));

        // Local constructMap is just identity
        {
            tgtConstructMap[Pstream::myProcNo(comm)] =
                identity(targetCoarseSize);
        }

        labelList tgtCompactMap(map.constructSize());

        {
            // Note that in special cases (e.g. 'appending' two AMIs) the
            // local size after distributing can be longer than the number
            // of faces. I.e. it duplicates elements.
            // Since we don't know this size instead we loop over all
            // reachable elements (using the local constructMap)

            const labelList& elemsMap =
                map.constructMap()[Pstream::myProcNo(comm)];
            for (const label fineElem : elemsMap)
            {
                label coarseElem = allRestrict[fineElem];
                tgtCompactMap[fineElem] = coarseElem;
            }
        }

        label compacti = targetCoarseSize;

        // Compact data from other processors
        forAll(map.constructMap(), proci)
        {
            if (proci != Pstream::myProcNo(comm))
            {
                // Combine entries that point to the same coarse element. All
                // elements now are remote data so we cannot use any local
                // data here - use allRestrict instead.
                const labelList& elems = map.constructMap()[proci];

                labelList& newConstructMap = tgtConstructMap[proci];
                newConstructMap.resize_nocopy(elems.size());

                if (elems.size())
                {
                    // Get the maximum target coarse size for this set of
                    // received data.
                    label remoteTargetCoarseSize = labelMin;
                    for (const label elemi : elems)
                    {
                        remoteTargetCoarseSize = max
                        (
                            remoteTargetCoarseSize,
                            allRestrict[elemi]
                        );
                    }
                    remoteTargetCoarseSize += 1;

                    // Combine locally data coming from proci
                    labelList oldToNew(remoteTargetCoarseSize, -1);
                    label newi = 0;

                    for (const label fineElem : elems)
                    {
                        // fineElem now points to section from proci
                        label coarseElem = allRestrict[fineElem];
                        if (oldToNew[coarseElem] == -1)
                        {
                            oldToNew[coarseElem] = newi;
                            tgtCompactMap[fineElem] = compacti;
                            newConstructMap[newi] = compacti++;
                            ++newi;
                        }
                        else
                        {
                            // Get compact index
                            label compacti = oldToNew[coarseElem];
                            tgtCompactMap[fineElem] = newConstructMap[compacti];
                        }
                    }
                    newConstructMap.resize(newi);
                }
            }
        }

        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, facei)
        {
            // All the elements contributing to facei. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[facei];
            const scalarList& weights = fineSrcWeights[facei];
            const scalar fineArea = fineSrcMagSf[facei];

            label coarseFacei = sourceRestrictAddressing[facei];

            labelList& newElems = srcAddress[coarseFacei];
            scalarList& newWeights = srcWeights[coarseFacei];

            forAll(elems, i)
            {
                label elemi = elems[i];
                label coarseElemi = tgtCompactMap[elemi];

                label index = newElems.find(coarseElemi);
                if (index == -1)
                {
                    newElems.append(coarseElemi);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }

        tgtMap.reset
        (
            new mapDistribute
            (
                compacti,
                std::move(tgtSubMap),
                std::move(tgtConstructMap),
                false,      //subHasFlip
                false,      //constructHasFlip
                comm
            )
        );
    }
    else
    {
        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, facei)
        {
            // All the elements contributing to facei. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[facei];
            const scalarList& weights = fineSrcWeights[facei];
            const scalar fineArea = fineSrcMagSf[facei];

            label coarseFacei = sourceRestrictAddressing[facei];

            labelList& newElems = srcAddress[coarseFacei];
            scalarList& newWeights = srcWeights[coarseFacei];

            forAll(elems, i)
            {
                const label elemi = elems[i];
                const label coarseElemi = targetRestrictAddressing[elemi];

                const label index = newElems.find(coarseElemi);
                if (index == -1)
                {
                    newElems.append(coarseElemi);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }
    }

    // Weights normalisation
    normaliseWeights
    (
        srcMagSf,
        "source",
        srcAddress,
        srcWeights,
        srcWeightsSum,
        true,
        false,
        -1,
        comm
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AMIInterpolation::AMIInterpolation
(
    const dictionary& dict,
    const bool reverseTarget
)
:
    requireMatch_(dict.getOrDefault("requireMatch", true)),
    reverseTarget_(dict.getOrDefault("reverseTarget", reverseTarget)),
    lowWeightCorrection_(dict.getOrDefault<scalar>("lowWeightCorrection", -1)),
    singlePatchProc_(-999),
    comm_(UPstream::worldComm),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcMapPtr_(nullptr),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    cache_(dict)
{}


Foam::AMIInterpolation::AMIInterpolation
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection
)
:
    requireMatch_(requireMatch),
    reverseTarget_(reverseTarget),
    lowWeightCorrection_(lowWeightCorrection),
    singlePatchProc_(-999),
    comm_(UPstream::worldComm),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcPatchPts_(),
    srcMapPtr_(nullptr),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtPatchPts_(),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    cache_()
{}


Foam::AMIInterpolation::AMIInterpolation
(
    const AMIInterpolation& fineAMI,
    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing
)
:
    requireMatch_(fineAMI.requireMatch_),
    reverseTarget_(fineAMI.reverseTarget_),
    lowWeightCorrection_(-1.0), // Deactivated?
    singlePatchProc_(fineAMI.singlePatchProc_),
    comm_(fineAMI.comm()),  // use fineAMI geomComm if present, comm otherwise
    geomComm_(),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcPatchPts_(),
    srcMapPtr_(nullptr),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    tgtPatchPts_(),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    cache_
    (
        fineAMI.cache(),
        fineAMI,
        sourceRestrictAddressing,
        targetRestrictAddressing
    )
{
    const label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    const label neighbourCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    if (debug & 2)
    {
        Pout<< "AMI: Creating addressing and weights as agglomeration of AMI :"
            << " source:" << fineAMI.srcAddress().size()
            << " target:" << fineAMI.tgtAddress().size()
            << " fineComm:" << fineAMI.comm()
            << " coarse source size:" << sourceCoarseSize
            << " neighbour source size:" << neighbourCoarseSize
            << endl;
    }

    if
    (
        fineAMI.srcAddress().size() != sourceRestrictAddressing.size()
     || fineAMI.tgtAddress().size() != targetRestrictAddressing.size()
    )
    {
        FatalErrorInFunction
            << "Size mismatch." << nl
            << "Source patch size:" << fineAMI.srcAddress().size() << nl
            << "Source agglomeration size:"
            << sourceRestrictAddressing.size() << nl
            << "Target patch size:" << fineAMI.tgtAddress().size() << nl
            << "Target agglomeration size:"
            << targetRestrictAddressing.size()
            << exit(FatalError);
    }


    // Agglomerate addresses and weights

    if (comm() != -1)
    {
        //Pout<< "** agglomerating srcAddress, tgtMap" << endl;
        //if (fineAMI.tgtMapPtr_.valid())
        //{
        //    const auto& fineTgtMap = fineAMI.tgtMapPtr_();
        //    Pout<< "    fineAMI.tgtMapPtr_ comm:" << fineTgtMap.comm()
        //        << "   procs:"
        //        <<  (
        //                fineTgtMap.comm() != -1
        //              ? UPstream::procID(fineTgtMap.comm())
        //              : labelList::null()
        //            )
        //        << endl;
        //}
        //else
        //{
        //    Pout<< "    NO fineAMI.tgtMapPtr_" << endl;
        //}
        //
        agglomerate
        (
            fineAMI.tgtMapPtr_,
            fineAMI.srcMagSf(),
            fineAMI.srcAddress(),
            fineAMI.srcWeights(),

            sourceRestrictAddressing,
            targetRestrictAddressing,

            srcMagSf_,
            srcAddress_,
            srcWeights_,
            srcWeightsSum_,
            tgtMapPtr_,
            comm()
        );

        //Pout<< "** agglomerating tgtAddress, srcMap" << endl;
        //if (fineAMI.srcMapPtr_.valid())
        //{
        //    const auto& fineSrcMap = fineAMI.srcMapPtr_();
        //    Pout<< "    fineAMI.srcMapPtr_ comm:" << fineSrcMap.comm()
        //        << "   procs:"
        //        <<  (
        //                fineSrcMap.comm() != -1
        //              ? UPstream::procID(fineSrcMap.comm())
        //              : labelList::null()
        //            )
        //        << endl;
        //}
        //else
        //{
        //    Pout<< "    NO fineAMI.srcMapPtr_" << endl;
        //}
        agglomerate
        (
            fineAMI.srcMapPtr_,
            fineAMI.tgtMagSf(),
            fineAMI.tgtAddress(),
            fineAMI.tgtWeights(),

            targetRestrictAddressing,
            sourceRestrictAddressing,

            tgtMagSf_,
            tgtAddress_,
            tgtWeights_,
            tgtWeightsSum_,
            srcMapPtr_,
            comm()
        );
    }
}


Foam::AMIInterpolation::AMIInterpolation(const AMIInterpolation& ami)
:
    requireMatch_(ami.requireMatch_),
    reverseTarget_(ami.reverseTarget_),
    lowWeightCorrection_(ami.lowWeightCorrection_),
    singlePatchProc_(ami.singlePatchProc_),
    comm_(ami.comm_),
    geomComm_(ami.geomComm_),   // ? steals communicator
    srcMagSf_(ami.srcMagSf_),
    srcAddress_(ami.srcAddress_),
    srcWeights_(ami.srcWeights_),
    srcWeightsSum_(ami.srcWeightsSum_),
    srcCentroids_(ami.srcCentroids_),
    srcMapPtr_(ami.srcMapPtr_.clone()),
    tgtMagSf_(ami.tgtMagSf_),
    tgtAddress_(ami.tgtAddress_),
    tgtWeights_(ami.tgtWeights_),
    tgtWeightsSum_(ami.tgtWeightsSum_),
    tgtCentroids_(ami.tgtCentroids_),
    tgtMapPtr_(ami.tgtMapPtr_.clone()),
    upToDate_(ami.upToDate_),
    cache_(ami.cache_)
{}


Foam::AMIInterpolation::AMIInterpolation(Istream& is)
:
    requireMatch_(readBool(is)),
    reverseTarget_(readBool(is)),
    lowWeightCorrection_(readScalar(is)),
    singlePatchProc_(readLabel(is)),
    comm_(readLabel(is)),   // either geomComm_ or comm_ from sending side

    srcMagSf_(is),
    srcAddress_(is),
    srcWeights_(is),
    srcWeightsSum_(is),
    srcCentroids_(is),
    //srcPatchPts_(is),
    srcMapPtr_(nullptr),

    tgtMagSf_(is),
    tgtAddress_(is),
    tgtWeights_(is),
    tgtWeightsSum_(is),
    tgtCentroids_(is),
    //tgtPatchPts_(is),
    tgtMapPtr_(nullptr),

    upToDate_(readBool(is)),

    cache_(is)
{
    // Hopefully no need to stream geomComm_ since only used in processor
    // agglomeration?

    if (singlePatchProc_ == -1 && comm_ != -1)
    {
        srcMapPtr_.reset(new mapDistribute(is));
        tgtMapPtr_.reset(new mapDistribute(is));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::AMIInterpolation::calculate
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (upToDate_)
    {
        return false;
    }

    addProfiling(ami, "AMIInterpolation::calculate");


    // Clear storage (only needed if src/tgt become zero size)
    {
        // Reset distribution maps so that a transition from distributed
        // to non-distributed, e.g. after dynamic load balancing, does not
        // leave stale non-null maps that corrupt GAMG agglomeration.
        srcMapPtr_.reset(nullptr);
        tgtMapPtr_.reset(nullptr);
        geomComm_.reset();

        if (srcMagSf_.size())
        {
            srcMagSf_.resize_nocopy(srcPatch.size());
        }
        srcAddress_.resize_nocopy(srcPatch.size());
        srcWeights_.resize_nocopy(srcPatch.size());
        srcWeightsSum_.resize_nocopy(srcPatch.size());
        if (srcCentroids_.size())
        {
            srcCentroids_.resize_nocopy(srcPatch.size());
        }

        if (tgtMagSf_.size())
        {
            tgtMagSf_.resize_nocopy(tgtPatch.size());
        }
        tgtAddress_.resize_nocopy(tgtPatch.size());
        tgtWeights_.resize_nocopy(tgtPatch.size());
        tgtWeightsSum_.resize_nocopy(tgtPatch.size());

        if (tgtCentroids_.size())
        {
            tgtCentroids_.resize_nocopy(tgtPatch.size());
        }
    }


    if (surfPtr)
    {
        srcPatchPts_ = srcPatch.points();
        projectPointsToSurface(surfPtr(), srcPatchPts_);
        tsrcPatch0_ = refPtr<primitivePatch>::New
        (
            SubList<face>(srcPatch),
            srcPatchPts_
        );

        tgtPatchPts_ = tgtPatch.points();
        projectPointsToSurface(surfPtr(), tgtPatchPts_);
        ttgtPatch0_ = refPtr<primitivePatch>::New
        (
            SubList<face>(tgtPatch),
            tgtPatchPts_
        );
    }
    else
    {
        tsrcPatch0_.cref(srcPatch);
        ttgtPatch0_.cref(tgtPatch);
    }

    // Note: use original communicator for statistics
    const label srcTotalSize = returnReduce
    (
        srcPatch.size(),
        sumOp<label>(),
        UPstream::msgType(),
        comm_
    );

    if (srcTotalSize == 0)
    {
        DebugInfo
            << "AMI: no source faces present - no addressing constructed"
            << endl;

        singlePatchProc_ = UPstream::myProcNo(comm_);

        return false;
    }

    const label tgtTotalSize = returnReduce
    (
        tgtPatch.size(),
        sumOp<label>(),
        UPstream::msgType(),
        comm_
    );

    // Calculate:
    // - which processors have faces
    // - allocates a communicator (geomComm_) for those
    // - if it is only one processor that holds all faces
    singlePatchProc_ = calcDistribution(srcPatch, tgtPatch, comm_, geomComm_);

    Info<< indent << "AMI: Patch source faces: " << srcTotalSize << nl
        << indent << "AMI: Patch target faces: " << tgtTotalSize << nl;

    if (distributed())
    {
        Info<< indent << "AMI: distributed" << endl;
    }

    DebugInfo
        << "AMI: patch proc:" << singlePatchProc_
        << endl;

    return true;
}


void Foam::AMIInterpolation::reset
(
    autoPtr<mapDistribute>&& srcToTgtMap,
    autoPtr<mapDistribute>&& tgtToSrcMap,
    labelListList&& srcAddress,
    scalarListList&& srcWeights,
    labelListList&& tgtAddress,
    scalarListList&& tgtWeights,
    const label singlePatchProc
)
{
    DebugInFunction<< endl;

    srcAddress_.transfer(srcAddress);
    srcWeights_.transfer(srcWeights);
    tgtAddress_.transfer(tgtAddress);
    tgtWeights_.transfer(tgtWeights);

    // Reset the sums of the weights
    srcWeightsSum_.resize_nocopy(srcWeights_.size());
    forAll(srcWeights_, facei)
    {
        srcWeightsSum_[facei] = sum(srcWeights_[facei]);
    }

    tgtWeightsSum_.resize_nocopy(tgtWeights_.size());
    forAll(tgtWeights_, facei)
    {
        tgtWeightsSum_[facei] = sum(tgtWeights_[facei]);
    }

    srcMapPtr_ = std::move(srcToTgtMap);
    tgtMapPtr_ = std::move(tgtToSrcMap);

    singlePatchProc_ = singlePatchProc;

    upToDate_ = true;
}


void Foam::AMIInterpolation::append
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
)
{
    addProfiling(ami, "AMIInterpolation::append");

    // Create a new interpolation
    // Note: set upToDate to false to force the cloned AMI to recalculate
    auto newPtr = clone();
    newPtr->upToDate(false);

    newPtr->calculate(srcPatch, tgtPatch);

    // If parallel then combine the mapDistribution and re-index
    if (distributed() && comm() != -1)
    {
        labelListList& srcSubMap = srcMapPtr_->subMap();
        labelListList& srcConstructMap = srcMapPtr_->constructMap();

        labelListList& tgtSubMap = tgtMapPtr_->subMap();
        labelListList& tgtConstructMap = tgtMapPtr_->constructMap();

        labelListList& newSrcSubMap = newPtr->srcMapPtr_->subMap();
        labelListList& newSrcConstructMap = newPtr->srcMapPtr_->constructMap();

        labelListList& newTgtSubMap = newPtr->tgtMapPtr_->subMap();
        labelListList& newTgtConstructMap = newPtr->tgtMapPtr_->constructMap();

        // Re-mapping/re-indexing - use max sizing
        labelList oldMapMap
        (
            max
            (
                srcMapPtr_->constructMapTotalSize(),
                tgtMapPtr_->constructMapTotalSize()
            )
        );
        labelList newMapMap
        (
            max
            (
                newPtr->srcMapPtr_->constructMapTotalSize(),
                newPtr->tgtMapPtr_->constructMapTotalSize()
            )
        );

        // Re-calculate the source indices
        {
            label total = 0;
            auto iter1 = oldMapMap.begin();
            auto iter2 = newMapMap.begin();

            forAll(srcSubMap, proci)
            {
                const label len1 = srcConstructMap[proci].size();
                const label len2 = newSrcConstructMap[proci].size();

                std::iota(iter1, (iter1 + len1), total);
                iter1 += len1;
                total += len1;

                std::iota(iter2, (iter2 + len2), total);
                iter2 += len2;
                total += len2;
            }
        }

        // Renumber the source indices
        {
            for (labelList& list : srcConstructMap)
            {
                for (label& value : list)
                {
                    value = oldMapMap[value];
                }
            }

            for (labelList& list : newSrcConstructMap)
            {
                for (label& value : list)
                {
                    value = newMapMap[value];
                }
            }

            for (labelList& list : tgtAddress_)
            {
                for (label& value : list)
                {
                    value = oldMapMap[value];
                }
            }

            for (labelList& list : newPtr->tgtAddress_)
            {
                for (label& value : list)
                {
                    value = newMapMap[value];
                }
            }
        }


        // Re-calculate the target indices
        {
            label total = 0;
            auto iter1 = oldMapMap.begin();
            auto iter2 = newMapMap.begin();

            forAll(srcSubMap, proci)
            {
                const label len1 = tgtConstructMap[proci].size();
                const label len2 = newTgtConstructMap[proci].size();

                std::iota(iter1, (iter1 + len1), total);
                iter1 += len1;
                total += len1;

                std::iota(iter2, (iter2 + len2), total);
                iter2 += len2;
                total += len2;
            }
        }

        // Renumber the target indices
        {
            for (labelList& list : tgtConstructMap)
            {
                for (label& value : list)
                {
                    value = oldMapMap[value];
                }
            }

            for (labelList& list : newTgtConstructMap)
            {
                for (label& value : list)
                {
                    value = newMapMap[value];
                }
            }

            for (labelList& list : srcAddress_)
            {
                for (label& value : list)
                {
                    value = oldMapMap[value];
                }
            }

            for (labelList& list : newPtr->srcAddress_)
            {
                for (label& value : list)
                {
                    value = newMapMap[value];
                }
            }
        }

        // Sum the construction sizes
        srcMapPtr_->constructSize() += newPtr->srcMapPtr_->constructSize();
        tgtMapPtr_->constructSize() += newPtr->tgtMapPtr_->constructSize();

        // Combine the maps
        forAll(srcSubMap, proci)
        {
            srcSubMap[proci].push_back(newSrcSubMap[proci]);
            srcConstructMap[proci].push_back(newSrcConstructMap[proci]);

            tgtSubMap[proci].push_back(newTgtSubMap[proci]);
            tgtConstructMap[proci].push_back(newTgtConstructMap[proci]);
        }
    }

    // Combine new and current source data
    forAll(srcMagSf_, srcFacei)
    {
        srcAddress_[srcFacei].push_back(newPtr->srcAddress()[srcFacei]);
        srcWeights_[srcFacei].push_back(newPtr->srcWeights()[srcFacei]);
        srcWeightsSum_[srcFacei] += newPtr->srcWeightsSum()[srcFacei];
    }

    // Combine new and current target data
    forAll(tgtMagSf_, tgtFacei)
    {
        tgtAddress_[tgtFacei].push_back(newPtr->tgtAddress()[tgtFacei]);
        tgtWeights_[tgtFacei].push_back(newPtr->tgtWeights()[tgtFacei]);
        tgtWeightsSum_[tgtFacei] += newPtr->tgtWeightsSum()[tgtFacei];
    }
}


void Foam::AMIInterpolation::normaliseWeights
(
    const bool conformal,
    const bool output
)
{
    normaliseWeights
    (
        srcMagSf_,
        "source",
        srcAddress_,
        srcWeights_,
        srcWeightsSum_,
        conformal,
        output,
        lowWeightCorrection_,
        comm()
    );

    normaliseWeights
    (
        tgtMagSf_,
        "target",
        tgtAddress_,
        tgtWeights_,
        tgtWeightsSum_,
        conformal,
        output,
        lowWeightCorrection_,
        comm()
    );
}


Foam::label Foam::AMIInterpolation::srcPointFace
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const vector& n,
    const label tgtFacei,
    point& tgtPoint
)
const
{
    const pointField& srcPoints = srcPatch.points();

    // Source face addresses that intersect target face tgtFacei
    const labelList& addr = tgtAddress_[tgtFacei];

    pointHit nearest;
    nearest.setDistance(GREAT);
    label nearestFacei = -1;

    for (const label srcFacei : addr)
    {
        const face& f = srcPatch[srcFacei];

        pointHit ray =
            f.ray(tgtPoint, n, srcPoints, intersection::algorithm::VISIBLE);

        if (ray.hit())
        {
            tgtPoint = ray.point();
            return srcFacei;
        }
        else if (ray.distance() < nearest.distance())
        {
            nearest = ray;
            nearestFacei = srcFacei;
        }
    }

    if (nearest.hit() || nearest.eligibleMiss())
    {
        tgtPoint = nearest.point();
        return nearestFacei;
    }

    return -1;
}


Foam::label Foam::AMIInterpolation::tgtPointFace
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const vector& n,
    const label srcFacei,
    point& srcPoint
)
const
{
    const pointField& tgtPoints = tgtPatch.points();

    pointHit nearest;
    nearest.setDistance(GREAT);
    label nearestFacei = -1;

    // Target face addresses that intersect source face srcFacei
    const labelList& addr = srcAddress_[srcFacei];

    for (const label tgtFacei : addr)
    {
        const face& f = tgtPatch[tgtFacei];

        pointHit ray =
            f.ray(srcPoint, n, tgtPoints, intersection::algorithm::VISIBLE);

        if (ray.hit())
        {
            srcPoint = ray.point();
            return tgtFacei;
        }
        const pointHit near = f.nearestPoint(srcPoint, tgtPoints);

        if (near.distance() < nearest.distance())
        {
            nearest = near;
            nearestFacei = tgtFacei;
        }
    }
    if (nearest.hit() || nearest.eligibleMiss())
    {
        srcPoint = nearest.point();
        return nearestFacei;
    }

    return -1;
}


bool Foam::AMIInterpolation::checkSymmetricWeights(const bool log) const
{
    if (UPstream::parRun() && this->distributed())
    {
        Log << "Checks only valid for serial running (currently)" << endl;

        return true;
    }

    bool symmetricSrc = true;

    Log << "    Checking for missing src face in tgt lists" << nl;

    forAll(srcAddress_, srcFacei)
    {
        const labelList& tgtIds = srcAddress_[srcFacei];
        for (const label tgtFacei : tgtIds)
        {
            if (!tgtAddress_[tgtFacei].found(srcFacei))
            {
                symmetricSrc = false;

                Log << "       srcFacei:" << srcFacei
                    << " not found in tgtToSrc list for tgtFacei:"
                    << tgtFacei << nl;
            }
        }
    }

    if (symmetricSrc)
    {
        Log << "    - symmetric" << endl;
    }

    bool symmetricTgt = true;

    Log << "    Checking for missing tgt face in src lists" << nl;

    forAll(tgtAddress_, tgtFacei)
    {
        const labelList& srcIds = tgtAddress_[tgtFacei];
        for (const label srcFacei : srcIds)
        {
            if (!srcAddress_[srcFacei].found(tgtFacei))
            {
                symmetricTgt = false;

                Log << "       tgtFacei:" << tgtFacei
                    << " not found in srcToTgt list for srcFacei:"
                    << srcFacei << nl;
            }
        }
    }

    if (symmetricTgt)
    {
        Log << "    - symmetric" << endl;
    }

    return symmetricSrc && symmetricTgt;
}


void Foam::AMIInterpolation::writeFaceConnectivity
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const labelListList& srcAddress
)
const
{
    OFstream os("faceConnectivity" + Foam::name(Pstream::myProcNo()) + ".obj");

    label pti = 1;

    forAll(srcAddress, i)
    {
        const labelList& addr = srcAddress[i];
        const point& srcPt = srcPatch.faceCentres()[i];

        for (const label tgtPti : addr)
        {
            const point& tgtPt = tgtPatch.faceCentres()[tgtPti];

            meshTools::writeOBJ(os, srcPt);
            meshTools::writeOBJ(os, tgtPt);

            os  << "l " << pti << " " << pti + 1 << endl;

            pti += 2;
        }
    }
}


void Foam::AMIInterpolation::write(Ostream& os) const
{
    os.writeEntry("AMIMethod", type());

    if (!requireMatch_)
    {
        os.writeEntry("requireMatch", requireMatch_);
    }

    if (reverseTarget_)
    {
        os.writeEntry("reverseTarget", reverseTarget_);
    }

    if (lowWeightCorrection_ > 0)
    {
        os.writeEntry("lowWeightCorrection", lowWeightCorrection_);
    }

    cache_.write(os);
}


bool Foam::AMIInterpolation::writeData(Ostream& os) const
{
    os  << requireMatch()
        << token::SPACE<< reverseTarget()
        << token::SPACE<< lowWeightCorrection()
        << token::SPACE<< singlePatchProc()
        << token::SPACE<< comm()    // either geomComm_ or comm_

        << token::SPACE<< srcMagSf()
        << token::SPACE<< srcAddress()
        << token::SPACE<< srcWeights()
        << token::SPACE<< srcWeightsSum()
        << token::SPACE<< srcCentroids()

        << token::SPACE<< tgtMagSf()
        << token::SPACE<< tgtAddress()
        << token::SPACE<< tgtWeights()
        << token::SPACE<< tgtWeightsSum()
        << token::SPACE<< tgtCentroids_

        << token::SPACE<< upToDate();

    cache_.writeData(os);

    if (distributed() && comm() != -1)
    {
        os  << token::SPACE<< srcMap()
            << token::SPACE<< tgtMap();
    }

    return os.good();
}


// ************************************************************************* //
