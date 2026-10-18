/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "nearestFaceAMI.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearestFaceAMI, 0);
    addToRunTimeSelectionTable(AMIInterpolation, nearestFaceAMI, dict);
    addToRunTimeSelectionTable(AMIInterpolation, nearestFaceAMI, component);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute> Foam::nearestFaceAMI::calcFaceMap
(
    const UList<nearKey>& localInfo,
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
) const
{
    const int communicator = this->comm();
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // Generate the list of processor bounding boxes for tgtPatch
    List<treeBoundBox> procBbs(numProc);
    if (!tgtPatch.empty())
    {
        procBbs[myProci] = treeBoundBox(tgtPatch.box());
    }
    boundBox::allGather(procBbs, communicator);

    // Identify which of my local src faces intersect with each processor
    // tgtPatch bb within the current match's search distance
    const pointField& srcCcs = srcPatch.faceCentres();
    List<DynamicList<label>> dynSendMap(numProc);

    forAll(localInfo, srcFacei)
    {
        // Test local srcPatch face centres against remote processor tgtPatch bb
        // using distance from initial pass

        // Misses carry a huge distance, so are tested against all
        // remote bounding boxes
        const scalar distSqr = localInfo[srcFacei].first();
        const point& srcCc = srcCcs[srcFacei];

        forAll(procBbs, proci)
        {
            const auto& bb = procBbs[proci];

            if (proci != myProci && bb.overlaps(srcCc, distSqr))
            {
                dynSendMap[proci].push_back(srcFacei);
            }
        }
    }

    // Convert dynamicList to labelList
    labelListList sendMap(numProc);
    forAll(sendMap, proci)
    {
        sendMap[proci].transfer(dynSendMap[proci]);

        if (debug)
        {
            Pout<< "send map - to proc " << proci << " sending "
                << sendMap[proci].size() << " elements" << endl;
        }
    }

    return autoPtr<mapDistribute>::New
    (
        std::move(sendMap),
        false,  // subHasFlip = false
        false,  // constructHasFlip = false
        communicator
    );
}


Foam::autoPtr<Foam::mapDistribute> Foam::nearestFaceAMI::calcDistributed
(
    const primitivePatch& src,
    const primitivePatch& tgt,
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght
) const
{
    const int communicator = comm();
    const auto myProci = UPstream::myProcNo(communicator);
    //const auto numProc = UPstream::nProcs(communicator);

    autoPtr<indexedOctree<treeType>> tgtTreePtr;
    if (tgt.size())
    {
        tgtTreePtr = this->createTree(tgt);
    }

    // Create global indexing for tgtPatch
    const globalIndex globalNumbering(tgt.size(), communicator);


    // First pass
    // ==========
    // For each srcPatch face, determine local match on tgtPatch

    // Identify local nearest matches
    pointField srcCcs(src.faceCentres());

    List<nearKey> localInfo(src.size(), { Foam::sqr(GREAT), -1 });

    if (tgt.size())
    {
        const auto& tgtTree = tgtTreePtr();

        forAll(srcCcs, srcCelli)
        {
            const point& srcCc = srcCcs[srcCelli];

            // With a search radius2 of GREAT all cells should receive a hit
            pointIndexHit test = tgtTree.findNearest(srcCc, GREAT);

            if (const auto& p = test.point(); test.hit())
            {
                localInfo[srcCelli] =
                (
                    nearKey
                    (
                        p.distSqr(srcCc),
                        globalNumbering.toGlobal(myProci, test.index())
                    )
                );
            }
        }
    }
    // else: no local tgtPatch faces - all srcPatch faces remain as misses
    // so that they [should be] set by remote tgtPatch faces


    // Second pass
    // ===========
    // Determine remote matches

    // Map returns labels of src patch faces to send to each proc
    autoPtr<mapDistribute> mapPtr = calcFaceMap(localInfo, src, tgt);
    mapDistribute& map = mapPtr();

    List<nearKey> remoteInfo(localInfo);
    map.distribute(remoteInfo);

    // Note: re-using srcCcs
    map.distribute(srcCcs);

    if (tgt.size())
    {
        const auto& tgtTree = tgtTreePtr();

        // Test remote srcPatch faces against local tgtPatch faces
        forAll(remoteInfo, i)
        {
            const point& srcCc = srcCcs[i];

            pointIndexHit test =
                tgtTree.findNearest(srcCc, remoteInfo[i].first());

            if (const auto& p = test.point(); test.hit())
            {
                nearKey testInfo
                (
                    p.distSqr(srcCc),
                    globalNumbering.toGlobal(myProci, test.index())
                );

                // Combine local into remote
                nearestEqOp()(remoteInfo[i], testInfo);
            }
        }
    }

    // Send back to originating processor. Choose best if sent to multiple
    // processors. Note that afterwards all unused entries have the unique
    // value nearestZero (distance < 0, index < 0).
    // This is used later on to see if the sample was sent to any processor.

    const nearKey nearestZero(-GREAT, -1);

    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        UList<labelPair>::null(),
        src.size(),
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        remoteInfo,
        nearestZero,
        nearestEqOp(),
        identityOp(),       // No flipping
        UPstream::msgType(),
        communicator
    );


    // Third pass
    // ==========
    // Combine local and remote info and filter out any connections that are
    // further away than threshold distance squared

    srcToTgtAddr.setSize(src.size());
    srcToTgtWght.setSize(src.size());
    forAll(srcToTgtAddr, srcFacei)
    {
        // Combine remote into local
        nearestEqOp()(localInfo[srcFacei], remoteInfo[srcFacei]);

        const auto [distSqr, tgtFacei] = localInfo[srcFacei];

        if
        (
            nearestEqOp::good_index(tgtFacei)
         && distSqr < maxDistance2_
        )
        {
            srcToTgtAddr[srcFacei].resize_fill(1, tgtFacei);
            srcToTgtWght[srcFacei].resize_fill(1, 1.0);
        }
    }

    List<Map<label>> cMap;
    return autoPtr<mapDistribute>::New
    (
        globalNumbering,
        srcToTgtAddr,
        cMap,
        UPstream::msgType(),
        communicator
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearestFaceAMI::nearestFaceAMI
(
    const dictionary& dict,
    const bool reverseTarget
)
:
    AMIInterpolation(dict, reverseTarget),
    maxDistance2_(dict.getOrDefault<scalar>("maxDistance2", GREAT))
{}


Foam::nearestFaceAMI::nearestFaceAMI
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection
)
:
    AMIInterpolation(requireMatch, reverseTarget, lowWeightCorrection),
    maxDistance2_(GREAT)
{}


Foam::nearestFaceAMI::nearestFaceAMI(const nearestFaceAMI& ami)
:
    AMIInterpolation(ami),
    maxDistance2_(ami.maxDistance2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nearestFaceAMI::calculate
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

    AMIInterpolation::calculate(srcPatch, tgtPatch, surfPtr);

    const auto& src = this->srcPatch0();
    const auto& tgt = this->tgtPatch0();

    // Set face area magnitudes
    srcMagSf_ = mag(src.faceAreas());
    tgtMagSf_ = mag(tgt.faceAreas());

    // TODO: implement symmetric calculation controls; assume yes for now
    bool symmetric_ = true;

    if (this->distributed() && comm() != -1)
    {
        tgtMapPtr_ =
            calcDistributed
            (
                src,
                tgt,
                srcAddress_,
                srcWeights_
            );

        if (symmetric_)
        {
            srcMapPtr_ =
                calcDistributed
                (
                    tgt,
                    src,
                    tgtAddress_,
                    tgtWeights_
                );
            }
    }
    else
    {
        srcAddress_.setSize(src.size());
        srcWeights_.setSize(src.size());

        if (symmetric_)
        {
            tgtAddress_.setSize(tgt.size());
            tgtWeights_.setSize(tgt.size());
        }

        const pointField& srcCcs = src.faceCentres();
        const pointField& tgtCcs = tgt.faceCentres();

        const auto tgtTreePtr = this->createTree(tgtPatch);
        const auto& tgtTree = tgtTreePtr();

        forAll(srcCcs, srcFacei)
        {
            const point& srcCc = srcCcs[srcFacei];
            const pointIndexHit hit = tgtTree.findNearest(srcCc, GREAT);

            if
            (
                hit.hit()
             && (magSqr(srcCc - tgtCcs[hit.index()]) < maxDistance2_)
            )
            {
                label tgtFacei = hit.index();
                srcAddress_[srcFacei] = labelList(1, tgtFacei);
                srcWeights_[srcFacei] = scalarList(1, 1.0);

                if (symmetric_)
                {
                    tgtAddress_[tgtFacei] = labelList(1, srcFacei);
                    tgtWeights_[tgtFacei] = scalarList(1, 1.0);
                }
            }
            else
            {
                if (debug)
                {
                    WarningInFunction
                        << "Unable to find target face for source face "
                        << srcFacei << endl;
                }
            }
        }

        if (symmetric_)
        {
            const auto srcTreePtr = this->createTree(srcPatch);
            const auto& srcTree = srcTreePtr();

            // Check that all source cells have connections and populate any
            // missing entries
            forAll(tgtWeights_, tgtCelli)
            {
                if (tgtAddress_[tgtCelli].empty())
                {
                    const point& tgtCc = tgtCcs[tgtCelli];
                    pointIndexHit hit = srcTree.findNearest(tgtCc, GREAT);

                    if
                    (
                        hit.hit()
                     && (magSqr(tgtCc - srcCcs[hit.index()]) < maxDistance2_)
                    )
                    {
                        tgtAddress_[tgtCelli] = labelList(1, hit.index());
                        tgtWeights_[tgtCelli] = scalarList(1, 1.0);
                    }
                }
                else
                {
                    if (debug)
                    {
                        WarningInFunction
                            << "Unable to find source face for target face "
                            << tgtCelli << endl;
                    }
                }
            }
        }
    }

    srcWeightsSum_.setSize(srcWeights_.size(), 1);
    tgtWeightsSum_.setSize(tgtWeights_.size(), 1);

    upToDate_ = true;

    return upToDate_;
}


void Foam::nearestFaceAMI::write(Ostream& os) const
{
    AMIInterpolation::write(os);
    os.writeEntryIfDifferent<scalar>("maxDistance2", GREAT, maxDistance2_);
}


// ************************************************************************* //
