/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "GAMGAgglomeration.H"
#include "GAMGInterface.H"
#include "processorGAMGInterface.H"
#include "cyclicLduInterface.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef FOAM_OFFLOAD
void Foam::GAMGAgglomeration::createSort
(
    const labelList& list,
    labelList& sort
)
{
    sort.setSize(list.size());

    const label n = list.size();

    labelRange range(n);

    // Initialize:
    // sort[i] = i
    apply_policy(n, [&](auto policy)
    {
        std::copy
        (
            policy,
            range.begin(),
            range.end(),
            sort.begin()
        );
    });

    // Copy keys
    labelList listTmp(list);

    // Sort keys and permute indices together
    apply_policy(n, [&](auto policy)
    {
        std::stable_sort
        (
            policy,
            sort.begin(),
            sort.end(),
            [&](const label a, const label b)
            {
                return listTmp[a] < listTmp[b];
            }
        );
    });
}


void Foam::GAMGAgglomeration::createTarget
(
    const labelList& list,
    const labelList& sort,
    labelList& target,
    labelList& targetStart
)
{
    const label n = list.size();
    if (!n)
    {
        target.clear();
        targetStart.resize_nocopy(1);
        targetStart[0] = 0;
        return;
    }

    labelList listSort(n);
    labelList flags(n);
    labelList offsets(n);
    labelRange range(n);

    //  Sort list
    apply_policy(n, [&](auto policy)
    {
        std::transform
        (
            policy,
            range.begin(),
            range.end(),
            listSort.begin(),
            [&](const label i)
            {
                return list[sort[i]];
            }
        );
    });

    // Mark the first entry of every unique target group
    apply_policy(n, [&](auto policy)
    {
        std::transform
        (
            policy,
            range.begin(),
            range.end(),
            flags.begin(),
            [&](const label i)
            {
                return (i == 0 || listSort[i] != listSort[i - 1]) ? 1 : 0;
            }
        );
    });

    // find total flags to get target size and offsets for scatter
    label targetSize = 0;
    apply_policy(n, [&](auto policy)
    {
        targetSize = std::reduce
        (
            policy,
            flags.begin(),
            flags.end(),
            label(0)
        );
    });

    // offsets[i] gives the compressed target index for group starts
    apply_policy(n, [&](auto policy)
    {
        std::exclusive_scan
        (
            policy,
            flags.begin(),
            flags.end(),
            offsets.begin(),
            label(0)
        );
    });

    target.setSize(targetSize);
    targetStart.setSize(targetSize + 1);

    // Scatter unique target values and their start offsets
    apply_policy(n, [&](auto policy)
    {
        std::for_each
        (
            policy,
            range.begin(),
            range.end(),
            [&](const label i)
            {
                if (flags[i])
                {
                    const label compacti = offsets[i];
                    target[compacti] = listSort[i];
                    targetStart[compacti] = i;
                }
            }
        );
    });
    targetStart[targetSize] = n;
 }


void Foam::GAMGAgglomeration::buildFullRestrictAddr
(
    const labelList& addr,
    const label level
)
{
    restrictSortAddressing_.set(level, new labelField(addr.size()));
    restrictTargetAddressing_.set(level, new labelField());
    restrictTargetStartAddressing_.set(level, new labelField());

    createSort
    (
        addr,
        restrictSortAddressing_[level]
    );

    createTarget
    (
        addr,
        restrictSortAddressing_[level],
        restrictTargetAddressing_[level],
        restrictTargetStartAddressing_[level]
    );
}


void Foam::GAMGAgglomeration::buildFullFaceRestrictAddr
(
    const labelList& addr,
    const label level
)
{
    faceRestrictSortAddressing_.set(level, new labelList(addr.size()));
    faceRestrictTargetAddressing_.set(level, new labelList());
    faceRestrictTargetStartAddressing_.set(level, new labelList());

    createSort
    (
        addr,
        faceRestrictSortAddressing_[level]
    );

    createTarget
    (
        addr,
        faceRestrictSortAddressing_[level],
        faceRestrictTargetAddressing_[level],
        faceRestrictTargetStartAddressing_[level]
    );
}


void Foam::GAMGAgglomeration::buildFullPatchFaceRestrictAddr
(
    const labelListList& addr,
    const label level
)
{
    patchFaceRestrictSortAddressing_.set
    (
        level,
        new labelListList(addr.size())
    );
    patchFaceRestrictTargetAddressing_.set
    (
        level,
        new labelListList(addr.size())
    );
    patchFaceRestrictTargetStartAddressing_.set
    (
        level,
        new labelListList(addr.size())
    );

    labelListList& patchFineToCoarseSortDevice =
        patchFaceRestrictSortAddressing_[level];
    labelListList& patchFineToCoarseTargetDevice =
        patchFaceRestrictTargetAddressing_[level];
    labelListList& patchFineToCoarseTargetStartDevice =
        patchFaceRestrictTargetStartAddressing_[level];

    forAll(addr, patchi)
    {
        createSort
        (
            addr[patchi],
            patchFineToCoarseSortDevice[patchi]
        );

        createTarget
        (
            addr[patchi],
            patchFineToCoarseSortDevice[patchi],
            patchFineToCoarseTargetDevice[patchi],
            patchFineToCoarseTargetStartDevice[patchi]
        );
    }
}


void Foam::GAMGAgglomeration::buildFullAddressing()
{
    if (restrictSortAddressing_.size())
    {
        return;
    }

    forAll(restrictAddressing_, leveli)
    {
        buildFullRestrictAddr(restrictAddressing_[leveli], leveli);
    }

    forAll(faceRestrictAddressing_, leveli)
    {
        buildFullFaceRestrictAddr(faceRestrictAddressing_[leveli], leveli);
    }

    forAll(patchFaceRestrictAddressing_, leveli)
    {
        buildFullPatchFaceRestrictAddr(patchFaceRestrictAddressing_[leveli],
            leveli);
    }
}
#endif


void Foam::GAMGAgglomeration::agglomerateLduAddressing
(
    const label fineLevelIndex
)
{
    const lduMesh& fineMesh = meshLevel(fineLevelIndex);
    const lduAddressing& fineMeshAddr = fineMesh.lduAddr();

    const labelUList& upperAddr = fineMeshAddr.upperAddr();
    const labelUList& lowerAddr = fineMeshAddr.lowerAddr();

    const label nFineFaces = upperAddr.size();

    // Get restriction map for current level
    const labelField& restrictMap = restrictAddressing(fineLevelIndex);

    if (min(restrictMap) == -1)
    {
        FatalErrorInFunction
            << "min(restrictMap) == -1" << exit(FatalError);
    }

    if (restrictMap.size() != fineMeshAddr.size())
    {
        FatalErrorInFunction
            << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << restrictMap.size()
            << " nEqns: " << fineMeshAddr.size()
            << abort(FatalError);
    }


    // Get the number of coarse cells
    const label nCoarseCells = nCells_[fineLevelIndex];

    // Storage for coarse cell neighbours and coefficients

    // Guess initial maximum number of neighbours in coarse cell
    label maxNnbrs = 10;

    // Number of faces for each coarse-cell
    labelList cCellnFaces(nCoarseCells, Zero);

    // Setup initial packed storage for coarse-cell faces
    labelList cCellFaces(maxNnbrs*nCoarseCells);

    // Create face-restriction addressing
    faceRestrictAddressing_.set(fineLevelIndex, new labelList(nFineFaces));
    labelList& faceRestrictAddr = faceRestrictAddressing_[fineLevelIndex];

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineFaces);

    // Counter for coarse faces
    label& nCoarseFaces = nFaces_[fineLevelIndex];
    nCoarseFaces = 0;

    // Loop through all fine faces
    forAll(upperAddr, fineFacei)
    {
        label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
        label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

        if (rmUpperAddr == rmLowerAddr)
        {
            // For each fine face inside of a coarse cell keep the address
            // of the cell corresponding to the face in the faceRestrictAddr
            // as a negative index
            faceRestrictAddr[fineFacei] = -(rmUpperAddr + 1);
        }
        else
        {
            // this face is a part of a coarse face

            label cOwn = rmUpperAddr;
            label cNei = rmLowerAddr;

            // get coarse owner and neighbour
            if (rmUpperAddr > rmLowerAddr)
            {
                cOwn = rmLowerAddr;
                cNei = rmUpperAddr;
            }

            // check the neighbour to see if this face has already been found
            label* ccFaces = &cCellFaces[maxNnbrs*cOwn];

            bool nbrFound = false;
            label& ccnFaces = cCellnFaces[cOwn];

            for (int i=0; i<ccnFaces; i++)
            {
                if (initCoarseNeighb[ccFaces[i]] == cNei)
                {
                    nbrFound = true;
                    faceRestrictAddr[fineFacei] = ccFaces[i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnFaces >= maxNnbrs)
                {
                    label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    cCellFaces.setSize(maxNnbrs*nCoarseCells);

                    forAllReverse(cCellnFaces, i)
                    {
                        label* oldCcNbrs = &cCellFaces[oldMaxNnbrs*i];
                        label* newCcNbrs = &cCellFaces[maxNnbrs*i];

                        for (int j=0; j<cCellnFaces[i]; j++)
                        {
                            newCcNbrs[j] = oldCcNbrs[j];
                        }
                    }

                    ccFaces = &cCellFaces[maxNnbrs*cOwn];
                }

                ccFaces[ccnFaces] = nCoarseFaces;
                initCoarseNeighb[nCoarseFaces] = cNei;
                faceRestrictAddr[fineFacei] = nCoarseFaces;
                ccnFaces++;

                // new coarse face created
                nCoarseFaces++;
            }
        }
    } // end for all fine faces


    // Renumber into upper-triangular order

    // All coarse owner-neighbour storage
    labelList coarseOwner(nCoarseFaces);
    labelList coarseNeighbour(nCoarseFaces);
    labelList coarseFaceMap(nCoarseFaces);

    label coarseFacei = 0;

    forAll(cCellnFaces, cci)
    {
        label* cFaces = &cCellFaces[maxNnbrs*cci];
        label ccnFaces = cCellnFaces[cci];

        for (int i=0; i<ccnFaces; i++)
        {
            coarseOwner[coarseFacei] = cci;
            coarseNeighbour[coarseFacei] = initCoarseNeighb[cFaces[i]];
            coarseFaceMap[cFaces[i]] = coarseFacei;
            coarseFacei++;
        }
    }

    forAll(faceRestrictAddr, fineFacei)
    {
        if (faceRestrictAddr[fineFacei] >= 0)
        {
            faceRestrictAddr[fineFacei] =
                coarseFaceMap[faceRestrictAddr[fineFacei]];
        }
    }


    // Create face-flip status
    faceFlipMap_.set(fineLevelIndex, new boolList(nFineFaces, false));
    boolList& faceFlipMap = faceFlipMap_[fineLevelIndex];


    forAll(faceRestrictAddr, fineFacei)
    {
        label coarseFacei = faceRestrictAddr[fineFacei];

        if (coarseFacei >= 0)
        {
            // Maps to coarse face
            label cOwn = coarseOwner[coarseFacei];
            label cNei = coarseNeighbour[coarseFacei];

            label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
            label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

            if (cOwn == rmUpperAddr && cNei == rmLowerAddr)
            {
                faceFlipMap[fineFacei] = true;
            }
            else if (cOwn == rmLowerAddr && cNei == rmUpperAddr)
            {
                //faceFlipMap[fineFacei] = false;
            }
            else
            {
                FatalErrorInFunction
                    << "problem."
                    << " fineFacei:" << fineFacei
                    << " rmUpperAddr:" << rmUpperAddr
                    << " rmLowerAddr:" << rmLowerAddr
                    << " coarseFacei:" << coarseFacei
                    << " cOwn:" << cOwn
                    << " cNei:" << cNei
                    << exit(FatalError);
            }
        }
    }



    // Clear the temporary storage for the coarse cell data
    cCellnFaces.clear();
    cCellFaces.clear();
    initCoarseNeighb.clear();
    coarseFaceMap.clear();


    // Create coarse-level interfaces

    // Get reference to fine-level interfaces
    const lduInterfacePtrsList& fineInterfaces = interfaceLevel(fineLevelIndex);

    nPatchFaces_.set
    (
        fineLevelIndex,
        new labelList(fineInterfaces.size(), Zero)
    );
    labelList& nPatchFaces = nPatchFaces_[fineLevelIndex];

    patchFaceRestrictAddressing_.set
    (
        fineLevelIndex,
        new labelListList(fineInterfaces.size())
    );
    labelListList& patchFineToCoarse =
        patchFaceRestrictAddressing_[fineLevelIndex];


    const label startOfRequests = UPstream::nRequests();

    // Initialise transfer of restrict addressing on the interface
    // The finest mesh uses patchAddr from the original lduAdressing.
    // the coarser levels create their own adressing for faceCells
    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            if (fineLevelIndex == 0)
            {
                fineInterfaces[inti].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    restrictMap,
                    fineMeshAddr.patchAddr(inti)
                );
            }
            else
            {
                fineInterfaces[inti].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    restrictMap
                );
            }
        }
    }

    UPstream::waitRequests(startOfRequests);


    // Add the coarse level
    meshLevels_.set
    (
        fineLevelIndex,
        new lduPrimitiveMesh
        (
            nCoarseCells,
            coarseOwner,
            coarseNeighbour,
            fineMesh.comm(),
            true
        )
    );

    lduInterfacePtrsList coarseInterfaces(fineInterfaces.size());

    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            tmp<labelField> restrictMapInternalField;

            // The finest mesh uses patchAddr from the original lduAdressing.
            // the coarser levels create their own adressing for faceCells
            if (fineLevelIndex == 0)
            {
                restrictMapInternalField =
                    fineInterfaces[inti].interfaceInternalField
                    (
                        restrictMap,
                        fineMeshAddr.patchAddr(inti)
                    );
            }
            else
            {
                restrictMapInternalField =
                    fineInterfaces[inti].interfaceInternalField
                    (
                        restrictMap
                    );
            }

            tmp<labelField> nbrRestrictMapInternalField =
                fineInterfaces[inti].internalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    restrictMap
                );

            coarseInterfaces.set
            (
                inti,
                GAMGInterface::New
                (
                    inti,
                    meshLevels_[fineLevelIndex].rawInterfaces(),
                    fineInterfaces[inti],
                    restrictMapInternalField(),
                    nbrRestrictMapInternalField(),
                    fineLevelIndex,
                    fineMesh.comm()
                ).ptr()
            );

            /* Same as below:
            coarseInterfaces.set
            (
                inti,
                GAMGInterface::New
                (
                    inti,
                    meshLevels_[fineLevelIndex].rawInterfaces(),
                    fineInterfaces[inti],
                    fineInterfaces[inti].interfaceInternalField(restrictMap),
                    fineInterfaces[inti].internalFieldTransfer
                    (
                        Pstream::commsTypes::nonBlocking,
                        restrictMap
                    ),
                    fineLevelIndex,
                    fineMesh.comm()
                ).ptr()
            );
            */

            nPatchFaces[inti] = coarseInterfaces[inti].faceCells().size();
            patchFineToCoarse[inti] = refCast<const GAMGInterface>
            (
                coarseInterfaces[inti]
            ).faceRestrictAddressing();
        }
    }

    meshLevels_[fineLevelIndex].addInterfaces
    (
        coarseInterfaces,
        lduPrimitiveMesh::nonBlockingSchedule<processorGAMGInterface>
        (
            coarseInterfaces
        )
    );

    #ifdef FOAM_OFFLOAD
    // faceFlipMap_.set(fineLevelIndex, new boolList(faceFlipMap));
    labelList faceRestrictAddressingTmp(faceRestrictAddr);
    buildFullFaceRestrictAddr(faceRestrictAddressingTmp, fineLevelIndex);

    labelListList patchFaceRestrictAddressingTmp(patchFineToCoarse.size());
    forAll(patchFineToCoarse, patchi)
    {
        patchFaceRestrictAddressingTmp[patchi] = patchFineToCoarse[patchi];
    }
    buildFullPatchFaceRestrictAddr
    (
        patchFaceRestrictAddressingTmp,
        fineLevelIndex
    );
    #endif

    if (debug & 2)
    {
        const auto& coarseAddr = meshLevels_[fineLevelIndex].lduAddr();

        Pout<< "GAMGAgglomeration :"
            << " agglomerated level " << fineLevelIndex
            << " from nCells:" << fineMeshAddr.size()
            << " nFaces:" << upperAddr.size()
            << " to nCells:" << nCoarseCells
            << " nFaces:" << nCoarseFaces << nl
            << "    lower:" << flatOutput(coarseAddr.lowerAddr()) << nl
            << "    upper:" << flatOutput(coarseAddr.upperAddr()) << nl
            << endl;
    }
}


void Foam::GAMGAgglomeration::procAgglomerateLduAddressing
(
    const label meshComm,
    const labelList& procAgglomMap,
    const labelList& procIDs,
    const label allMeshComm,

    const label levelIndex
)
{
    // - Assemble all the procIDs in meshComm onto a single master
    //   (procIDs[0]). This constructs a new communicator ('comm') first.
    // - The master communicates with neighbouring masters using
    //   allMeshComm

    const lduMesh& myMesh = meshLevels_[levelIndex-1];
    const label nOldInterfaces = myMesh.interfaces().size();

    procAgglomMap_.set(levelIndex, new labelList(procAgglomMap));
    agglomProcIDs_.set(levelIndex, new labelList(procIDs));
    procCommunicator_[levelIndex] = allMeshComm;

    procAgglomCommunicator_.set
    (
        levelIndex,
        new UPstream::communicator
        (
            meshComm,
            procIDs
        )
    );
    const label comm = agglomCommunicator(levelIndex);

    // These could only be set on the master procs but it is
    // quite convenient to also have them on the slaves
    procCellOffsets_.set(levelIndex, new labelList(0));
    procFaceMap_.set(levelIndex, new labelListList(0));
    procBoundaryMap_.set(levelIndex, new labelListList(0));
    procBoundaryFaceMap_.set(levelIndex, new labelListListList(0));


    // Collect meshes. Does no renumbering
    PtrList<lduPrimitiveMesh> otherMeshes;
    lduPrimitiveMesh::gather(comm, myMesh, otherMeshes);

    if (Pstream::myProcNo(meshComm) == procIDs[0])
    {
        // Combine all addressing

        labelList procFaceOffsets;
        meshLevels_.set
        (
            levelIndex-1,
            new lduPrimitiveMesh
            (
                allMeshComm,
                procAgglomMap,

                procIDs,
                myMesh,
                otherMeshes,

                procCellOffsets_[levelIndex],
                procFaceOffsets,
                procFaceMap_[levelIndex],
                procBoundaryMap_[levelIndex],
                procBoundaryFaceMap_[levelIndex]
            )
        );
    }


    // Scatter the procBoundaryMap back to the originating processor
    // so it knows which proc boundaries are to be kept. This is used
    // so we only send over interfaceFields on kept processors (see
    // GAMGSolver::procAgglomerateMatrix)
    // TBD: using sub-communicator here (instead of explicit procIDs). Should
    //      use sub-communicators more in other places.
    {
        const CompactListList<label> data
        (
            CompactListList<label>::pack<labelList>
            (
                procBoundaryMap_[levelIndex]
            )
        );
        const labelList localSizes = data.localSizes();
        const labelList& localStarts = data.offsets();

        // Make space
        procBoundaryMap_[levelIndex].setSize(procIDs.size());
        labelList& bMap = procBoundaryMap_[levelIndex][Pstream::myProcNo(comm)];
        bMap.setSize(nOldInterfaces);

        // Scatter relevant section to originating processor
        UPstream::mpiScatterv
        (
            data.values().cdata(),

            // Pass as List<int> for MPI
            ConstPrecisionAdaptor<int, label, List>(localSizes).cref(),
            ConstPrecisionAdaptor<int, label, List>(localStarts).cref(),

            bMap.data(),
            bMap.size(),
            comm
        );
    }


    // Combine restrict addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    procAgglomerateRestrictAddressing
    (
        meshComm,
        procIDs,
        levelIndex
    );

    if (Pstream::myProcNo(meshComm) != procIDs[0])
    {
        clearLevel(levelIndex);
    }
}


void Foam::GAMGAgglomeration::procAgglomerateRestrictAddressing
(
    const label comm,
    const labelList& procIDs,
    const label levelIndex
)
{
    const bool master =
    (
        UPstream::myProcNo(comm) == (procIDs.empty() ? 0 : procIDs[0])
    );

    // Determine the fine/coarse sizes (offsets) for gathering
    labelList fineOffsets;
    labelList coarseOffsets;

    {
        List<labelPair> sizes = globalIndex::listGatherValues
        (
            comm,
            procIDs,
            labelPair
            (
                // fine
                restrictAddressing_[levelIndex].size(),
                // coarse
                nCells_[levelIndex]
            ),
            UPstream::msgType(),
            UPstream::commsTypes::scheduled
        );

        // Calculate offsets, as per globalIndex::calcOffsets()
        // but extracting from the pair
        if (master && !sizes.empty())
        {
            const label len = sizes.size();

            fineOffsets.resize(len+1);
            coarseOffsets.resize(len+1);

            label fineCount = 0;
            label coarseCount = 0;

            for (label i = 0; i < len; ++i)
            {
                fineOffsets[i] = fineCount;
                fineCount += sizes[i].first();

                coarseOffsets[i] = coarseCount;
                coarseCount += sizes[i].second();
            }

            fineOffsets[len] = fineCount;
            coarseOffsets[len] = coarseCount;
        }
    }


    // (cell)restrictAddressing
    labelList procRestrictAddressing;
    if (master)
    {
        // pre-size on master
        procRestrictAddressing.resize(fineOffsets.back());
    }
    globalIndex::gather
    (
        fineOffsets,
        comm,
        procIDs,
        restrictAddressing_[levelIndex],
        procRestrictAddressing,
        UPstream::msgType(),
        UPstream::commsTypes::nonBlocking
    );

    if (master)
    {
        nCells_[levelIndex] = coarseOffsets.back();  // ie, totalSize()

        // Renumber consecutively
        for (label proci = 1; proci < procIDs.size(); ++proci)
        {
            SubList<label> procSlot
            (
                procRestrictAddressing,
                fineOffsets[proci+1]-fineOffsets[proci],
                fineOffsets[proci]
            );

            // procSlot += coarseOffsets[proci];
            forAll(procSlot, i)
            {
                procSlot[i] += coarseOffsets[proci];
            }
        }

        restrictAddressing_[levelIndex].transfer(procRestrictAddressing);
    }
}


void Foam::GAMGAgglomeration::combineLevels(const label curLevel)
{
    label prevLevel = curLevel - 1;

    // Set the previous level nCells to the current
    nCells_[prevLevel] = nCells_[curLevel];
    nFaces_[prevLevel] = nFaces_[curLevel];

    // Map the restrictAddressing from the coarser level into the previous
    // finer level

    const labelList& curResAddr = restrictAddressing_[curLevel];
    labelList& prevResAddr = restrictAddressing_[prevLevel];

    const labelList& curFaceResAddr = faceRestrictAddressing_[curLevel];
    labelList& prevFaceResAddr = faceRestrictAddressing_[prevLevel];
    const boolList& curFaceFlipMap = faceFlipMap_[curLevel];
    boolList& prevFaceFlipMap = faceFlipMap_[prevLevel];

    forAll(prevFaceResAddr, i)
    {
        if (prevFaceResAddr[i] >= 0)
        {
            label fineFacei = prevFaceResAddr[i];
            prevFaceResAddr[i] = curFaceResAddr[fineFacei];
            prevFaceFlipMap[i] = curFaceFlipMap[fineFacei];
        }
        else
        {
            label fineFacei = -prevFaceResAddr[i] - 1;
            prevFaceResAddr[i] = -curResAddr[fineFacei] - 1;
            prevFaceFlipMap[i] = curFaceFlipMap[fineFacei];
        }
    }

    #ifdef FOAM_OFFLOAD
    faceRestrictSortAddressing_.set(curLevel, nullptr);
    faceRestrictTargetAddressing_.set(curLevel, nullptr);
    faceRestrictTargetStartAddressing_.set(curLevel, nullptr);

    labelList restrictAddressingTmp(prevResAddr);
    // faceFlipMap_[prevLevel] = prevFaceFlipMap;
    restrictSortAddressing_.set(prevLevel, new labelField(prevResAddr.size()));

    createSort
    (
        restrictAddressingTmp,
        restrictSortAddressing_[prevLevel]
    );

    createTarget
    (
        restrictAddressingTmp,
        restrictSortAddressing_[prevLevel],
        restrictTargetAddressing_[prevLevel],
        restrictTargetStartAddressing_[prevLevel]
    );
    #endif

    // Delete the restrictAddressing for the coarser level
    faceRestrictAddressing_.set(curLevel, nullptr);
    faceFlipMap_.set(curLevel, nullptr);
    #ifdef FOAM_OFFLOAD
    restrictSortAddressing_.set(curLevel, nullptr);
    restrictTargetAddressing_.set(curLevel, nullptr);
    restrictTargetStartAddressing_.set(curLevel, nullptr);
    #endif

    forAll(prevResAddr, i)
    {
        prevResAddr[i] = curResAddr[prevResAddr[i]];
    }

    const labelListList& curPatchFaceResAddr =
        patchFaceRestrictAddressing_[curLevel];
    labelListList& prevPatchFaceResAddr =
        patchFaceRestrictAddressing_[prevLevel];

    forAll(prevPatchFaceResAddr, inti)
    {
        const labelList& curResAddr = curPatchFaceResAddr[inti];
        labelList& prevResAddr = prevPatchFaceResAddr[inti];
        forAll(prevResAddr, i)
        {
            label fineFacei = prevResAddr[i];
            prevResAddr[i] = curResAddr[fineFacei];
        }

        #ifdef FOAM_OFFLOAD
        labelList patchFaceRestrictAddressingTmp(prevResAddr);

        createSort
        (
            patchFaceRestrictAddressingTmp,
            patchFaceRestrictSortAddressing_[prevLevel][inti]
        );

        createTarget
        (
            patchFaceRestrictAddressingTmp,
            patchFaceRestrictSortAddressing_[prevLevel][inti],
            patchFaceRestrictTargetAddressing_[prevLevel][inti],
            patchFaceRestrictTargetStartAddressing_[prevLevel][inti]
        );
        #endif
    }

    // Delete the restrictAddressing for the coarser level
    restrictAddressing_.set(curLevel, nullptr);

    // Patch faces
    nPatchFaces_[prevLevel] = nPatchFaces_[curLevel];



    // Adapt the restrict addressing for the patches
    const lduInterfacePtrsList& curInterLevel =
        meshLevels_[curLevel].rawInterfaces();
    const lduInterfacePtrsList& prevInterLevel =
        meshLevels_[prevLevel].rawInterfaces();

    forAll(prevInterLevel, inti)
    {
        if (prevInterLevel.set(inti))
        {
            GAMGInterface& prevInt = refCast<GAMGInterface>
            (
                const_cast<lduInterface&>
                (
                    prevInterLevel[inti]
                )
            );
            const GAMGInterface& curInt = refCast<const GAMGInterface>
            (
                curInterLevel[inti]
            );
            prevInt.combine(curInt);
        }
    }

    // Delete the matrix addressing and coefficients from the previous level
    // and replace with the corresponding entry from the coarser level
    meshLevels_.set(prevLevel, meshLevels_.set(curLevel, nullptr));
}


void Foam::GAMGAgglomeration::calculateRegionMaster
(
    const label comm,
    const labelList& procAgglomMap,
    labelList& masterProcs,
    List<label>& agglomProcIDs
)
{
    // Determine the master processors
    Map<label> agglomToMaster(procAgglomMap.size());

    forAll(procAgglomMap, proci)
    {
        const label coarsei = procAgglomMap[proci];

        auto iter = agglomToMaster.find(coarsei);
        if (iter.good())
        {
            iter.val() = min(iter.val(), proci);
        }
        else
        {
            agglomToMaster.insert(coarsei, proci);
        }
    }

    masterProcs.setSize(agglomToMaster.size());
    forAllConstIters(agglomToMaster, iter)
    {
        masterProcs[iter.key()] = iter.val();
    }


    // Collect all the processors in my agglomeration
    label myProcID = Pstream::myProcNo(comm);
    label myAgglom = procAgglomMap[myProcID];

    // Get all processors agglomerating to the same coarse
    // processor
    agglomProcIDs = findIndices(procAgglomMap, myAgglom);

    // Make sure the master is the first element.
    const label index =
        agglomProcIDs.find(agglomToMaster[myAgglom]);

    std::swap(agglomProcIDs[0], agglomProcIDs[index]);
}


// ************************************************************************* //
