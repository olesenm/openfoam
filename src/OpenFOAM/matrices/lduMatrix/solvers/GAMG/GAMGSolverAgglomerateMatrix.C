/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "GAMGSolver.H"
#include "GAMGInterfaceField.H"
#include "processorLduInterfaceField.H"
#include "processorGAMGInterfaceField.H"
#ifdef FOAM_OFFLOAD
#include "GAMGSolverAgglomerateMatrixF.H"
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::agglomerateMatrix
(
    const label fineLevelIndex,
    const lduMesh& coarseMesh,
    const lduInterfacePtrsList& coarseMeshInterfaces
)
{
    // Get fine matrix
    const lduMatrix& fineMatrix = matrixLevel(fineLevelIndex);

    if (UPstream::myProcNo(fineMatrix.mesh().comm()) != -1)
    {
        const label nCoarseFaces = agglomeration_.nFaces(fineLevelIndex);
        const label nCoarseCells = agglomeration_.nCells(fineLevelIndex);

        // Set the coarse level matrix
        matrixLevels_.set
        (
            fineLevelIndex,
            new lduMatrix(coarseMesh)
        );
        lduMatrix& coarseMatrix = matrixLevels_[fineLevelIndex];


        // Coarse matrix diagonal initialised by restricting the finer mesh
        // diagonal. Note that we size with the cached coarse nCells and not
        // the actual coarseMesh size since this might be dummy when processor
        // agglomerating.
        scalarField& coarseDiag = coarseMatrix.diag(nCoarseCells);

        // Use manual restriction to avoid resize in restrictField
        const labelList& cellRestrictAddr =
            agglomeration_.restrictAddressing(fineLevelIndex);

        coarseDiag = Zero;
        forAll(fineMatrix.diag(), i)
        {
            coarseDiag[cellRestrictAddr[i]] += fineMatrix.diag()[i];
        }

        // Get reference to fine-level interfaces
        const lduInterfaceFieldPtrsList& fineInterfaces =
            interfaceLevel(fineLevelIndex);

        // Create coarse-level interfaces
        primitiveInterfaceLevels_.set
        (
            fineLevelIndex,
            new PtrList<lduInterfaceField>(fineInterfaces.size())
        );

        PtrList<lduInterfaceField>& coarsePrimInterfaces =
            primitiveInterfaceLevels_[fineLevelIndex];

        interfaceLevels_.set
        (
            fineLevelIndex,
            new lduInterfaceFieldPtrsList(fineInterfaces.size())
        );

        lduInterfaceFieldPtrsList& coarseInterfaces =
            interfaceLevels_[fineLevelIndex];

        // Get coarse level interface sizes
        const labelList& nPatchFaces =
            agglomeration_.nPatchFaces(fineLevelIndex);

        // Calculate total size and offsets for flat storage
        // Only include SET interfaces (skip unset wall boundaries)
        interfaceLevelsOffsets_.set
        (
            fineLevelIndex,
            new labelList(fineInterfaces.size() + 1, 0)
        );
        labelList& coeffsOffsets = interfaceLevelsOffsets_[fineLevelIndex];

        label totalSize = 0;
        forAll(nPatchFaces, inti)
        {
            coeffsOffsets[inti+1] = coeffsOffsets[inti];
            if (fineInterfaces.set(inti))
            {
                coeffsOffsets[inti+1] += nPatchFaces[inti];
                totalSize += nPatchFaces[inti];
            }
        }

        // Set coarse-level boundary coefficients (flat storage)
        interfaceLevelsBouCoeffs_.set
        (
            fineLevelIndex,
            new scalarField(totalSize, Zero)
        );
        scalarField& coarseInterfaceBouCoeffs =
            interfaceLevelsBouCoeffs_[fineLevelIndex];

        // Set coarse-level internal coefficients (flat storage)
        interfaceLevelsIntCoeffs_.set
        (
            fineLevelIndex,
            new scalarField(totalSize, Zero)
        );
        scalarField& coarseInterfaceIntCoeffs =
            interfaceLevelsIntCoeffs_[fineLevelIndex];

        // Add the coarse level
        agglomerateInterfaceCoefficients
        (
            fineLevelIndex,
            coarseMeshInterfaces,
            coarsePrimInterfaces,
            coarseInterfaces,
            coarseInterfaceBouCoeffs,
            coarseInterfaceIntCoeffs,
            coeffsOffsets
        );


        // Get face restriction map for current level
        const labelList& faceRestrictAddr =
            agglomeration_.faceRestrictAddressing(fineLevelIndex);
        const boolList& faceFlipMap =
            agglomeration_.faceFlipMap(fineLevelIndex);
        #ifdef FOAM_OFFLOAD
        const labelList& faceRestrictSortAddr =
            agglomeration_.faceRestrictSortAddressing(fineLevelIndex);
        const labelList& faceRestrictTargetAddr =
            agglomeration_.faceRestrictTargetAddressing(fineLevelIndex);
        const labelList& faceRestrictTargetStartAddr =
            agglomeration_.faceRestrictTargetStartAddressing(fineLevelIndex);
        #endif

        // Check if matrix is asymmetric and if so agglomerate both upper
        // and lower coefficients ...
        if (fineMatrix.hasLower())
        {
            // Get off-diagonal matrix coefficients
            const scalarField& fineUpper = fineMatrix.upper();
            const scalarField& fineLower = fineMatrix.lower();

            // Coarse matrix upper coefficients. Note passed in size
            scalarField& coarseUpper = coarseMatrix.upper(nCoarseFaces);
            scalarField& coarseLower = coarseMatrix.lower(nCoarseFaces);

            #ifdef FOAM_OFFLOAD
            labelRange faceRange(0, faceRestrictTargetAddr.size());

            apply_policy(faceRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    faceRange.begin(),
                    faceRange.end(),
                    GAMG::asymAgglomerate
                    (
                        coarseUpper,
                        coarseLower,
                        fineUpper,
                        fineLower,
                        faceFlipMap,
                        faceRestrictSortAddr,
                        faceRestrictTargetAddr,
                        faceRestrictTargetStartAddr
                    )
                );
            });

            apply_policy(faceRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    faceRange.begin(),
                    faceRange.end(),
                    GAMG::diagAsymAgglomerate
                    (
                        coarseDiag,
                        fineUpper,
                        fineLower,
                        faceRestrictSortAddr,
                        faceRestrictTargetAddr,
                        faceRestrictTargetStartAddr
                    )
                );
            });
            #else
            forAll(faceRestrictAddr, fineFacei)
            {
                label cFace = faceRestrictAddr[fineFacei];

                if (cFace >= 0)
                {
                    // Check the orientation of the fine-face relative to the
                    // coarse face it is being agglomerated into
                    if (!faceFlipMap[fineFacei])
                    {
                        coarseUpper[cFace] += fineUpper[fineFacei];
                        coarseLower[cFace] += fineLower[fineFacei];
                    }
                    else
                    {
                        coarseUpper[cFace] += fineLower[fineFacei];
                        coarseLower[cFace] += fineUpper[fineFacei];
                    }
                }
                else
                {
                    // Add the fine face coefficients into the diagonal.
                    coarseDiag[-1 - cFace] +=
                        fineUpper[fineFacei] + fineLower[fineFacei];
                }
            }
            #endif
        }
        else // ... Otherwise it is symmetric so agglomerate just the upper
        {
            // Get off-diagonal matrix coefficients
            const scalarField& fineUpper = fineMatrix.upper();

            // Coarse matrix upper coefficients
            scalarField& coarseUpper = coarseMatrix.upper(nCoarseFaces);

            #ifdef FOAM_OFFLOAD
            labelRange faceRange(0, faceRestrictTargetAddr.size());

            apply_policy(faceRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    faceRange.begin(),
                    faceRange.end(),
                    GAMG::symAgglomerate
                    (
                        coarseUpper,
                        fineUpper,
                        faceRestrictSortAddr,
                        faceRestrictTargetAddr,
                        faceRestrictTargetStartAddr
                    )
                );
            });

            apply_policy(faceRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    faceRange.begin(),
                    faceRange.end(),
                    GAMG::diagSymAgglomerate
                    (
                        coarseDiag,
                        fineUpper,
                        faceRestrictSortAddr,
                        faceRestrictTargetAddr,
                        faceRestrictTargetStartAddr
                    )
                );
            });
            #else
            forAll(faceRestrictAddr, fineFacei)
            {
                label cFace = faceRestrictAddr[fineFacei];

                if (cFace >= 0)
                {
                    coarseUpper[cFace] += fineUpper[fineFacei];
                }
                else
                {
                    // Add the fine face coefficient into the diagonal.
                    coarseDiag[-1 - cFace] += 2*fineUpper[fineFacei];
                }
            }
            #endif
        }
    }
}


void Foam::GAMGSolver::agglomerateInterfaceCoefficients
(
    const label fineLevelIndex,
    const lduInterfacePtrsList& coarseMeshInterfaces,
    PtrList<lduInterfaceField>& coarsePrimInterfaces,
    lduInterfaceFieldPtrsList& coarseInterfaces,
    scalarField& coarseInterfaceBouCoeffs,
    scalarField& coarseInterfaceIntCoeffs,
    const labelList& coeffsOffsets
) const
{
    // Get reference to fine-level interfaces
    const lduInterfaceFieldPtrsList& fineInterfaces =
        interfaceLevel(fineLevelIndex);

    // Get reference to fine-level coefficients (flat storage)
    const scalarField& fineInterfaceBouCoeffs =
        interfaceBouCoeffsLevelData(fineLevelIndex);

    const scalarField& fineInterfaceIntCoeffs =
        interfaceIntCoeffsLevelData(fineLevelIndex);

    const labelList& fineCoeffsOffsets =
        interfaceCoeffsLevelOffsets(fineLevelIndex);

    const labelListList& patchFineToCoarse =
        agglomeration_.patchFaceRestrictAddressing(fineLevelIndex);

    // Add the coarse level
    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            const GAMGInterface& coarseInterface =
                refCast<const GAMGInterface>
                (
                    coarseMeshInterfaces[inti]
                );

            coarsePrimInterfaces.set
            (
                inti,
                GAMGInterfaceField::New
                (
                    coarseInterface,
                    fineInterfaces[inti]
                ).ptr()
            );
            coarseInterfaces.set
            (
                inti,
                &coarsePrimInterfaces[inti]
            );

            const labelList& faceRestrictAddressing = patchFineToCoarse[inti];

            // Get SubFields for this interface in flat storage
            const label fineStart = fineCoeffsOffsets[inti];
            const label fineSize = fineCoeffsOffsets[inti+1] - fineStart;
            SubField<scalar> fineBou(fineInterfaceBouCoeffs, fineSize, fineStart);
            SubField<scalar> fineInt(fineInterfaceIntCoeffs, fineSize, fineStart);

            const label coarseStart = coeffsOffsets[inti];
            const label coarseSize = coeffsOffsets[inti+1] - coarseStart;
            SubField<scalar> coarseBou
            (
                coarseInterfaceBouCoeffs,
                coarseSize,
                coarseStart
            );
            SubField<scalar> coarseInt
            (
                coarseInterfaceIntCoeffs,
                coarseSize,
                coarseStart
            );

            #ifdef FOAM_OFFLOAD
            const labelListList& fToCSort =
                agglomeration_.patchFaceRestrictSortAddressing(fineLevelIndex);
            const labelListList& fToCTarget =
                agglomeration_.patchFaceRestrictTargetAddressing
                (
                    fineLevelIndex
                );
            const labelListList& fToCTargetStart =
                agglomeration_.patchFaceRestrictTargetStartAddressing
                (
                    fineLevelIndex
                );

            const labelList& faceRestrictSortAddressing = fToCSort[inti];
            const labelList& faceRestrictTargetAddressing = fToCTarget[inti];
            const labelList& faceRestrictTargetStartAddressing =
                fToCTargetStart[inti];

            agglomeration_.restrictField
            (
                coarseBou,
                fineBou,
                faceRestrictSortAddressing,
                faceRestrictTargetAddressing,
                faceRestrictTargetStartAddressing
            );
            #else
            agglomeration_.restrictField
            (
                coarseBou,
                fineBou,
                faceRestrictAddressing
            );
            #endif

            #ifdef FOAM_OFFLOAD
            agglomeration_.restrictField
            (
                coarseInt,
                fineInt,
                faceRestrictSortAddressing,
                faceRestrictTargetAddressing,
                faceRestrictTargetStartAddressing
            );
            #else
            agglomeration_.restrictField
            (
                coarseInt,
                fineInt,
                faceRestrictAddressing
            );
            #endif
        }
    }
}


void Foam::GAMGSolver::gatherMatrices
(
    const label destLevel,
    const label comm,

    // Local matrix
    const lduMatrix& mat,
    const scalarField& interfaceBouCoeffs,
    const scalarField& interfaceIntCoeffs,
    const labelList& coeffsOffsets,
    const lduInterfaceFieldPtrsList& interfaces,

    // Remote matrices
    PtrList<lduMatrix>& otherMats,
    PtrList<FieldField<Field, scalar>>& otherBouCoeffs,
    PtrList<FieldField<Field, scalar>>& otherIntCoeffs,
    PtrList<PtrList<lduInterfaceField>>& otherInterfaces
) const
{
    if (debug & 2)
    {
        const auto& procIDs = UPstream::procID(comm);

        Pout<< "GAMGSolver::gatherMatrices :"
            << " collecting matrices from procs:" << procIDs
            << " using comm:" << comm << endl;
    }

    const auto& boundaryMap = agglomeration_.boundaryMap(destLevel);

    PstreamBuffers pBufs(comm);

    // Send to master
    if (!UPstream::master(comm))
    {
        // Mark valid interfaces
        // -1   : not set
        // >= 0 : coupled interface (might also be unmerged processor boundary)
        //
        // Note: most processor interfaces will disappear. Originally
        // we did not know which ones were kept but this is now stored
        // on the boundaryMap (even on the slave processors). So we can
        // already filter here and avoid sending across typeNames etc.

        const label proci = UPstream::myProcNo(comm);

        // All interfaceBouCoeffs need to be sent across
        bitSet validCoeffs(interfaces.size());
        forAll(interfaces, intI)
        {
            if (coeffsOffsets[intI+1] > coeffsOffsets[intI])
            {
                validCoeffs.set(intI);
            }
        }

        // Only preserved interfaces need to be sent across
        bitSet validInterface(interfaces.size());
        forAll(interfaces, intI)
        {
            const label allIntI = boundaryMap[proci][intI];
            if (interfaces.set(intI) && allIntI != -1)
            {
                validInterface.set(intI);
            }
        }

        UOPstream toMaster(UPstream::masterNo(), pBufs);

        toMaster
            << mat
            << token::SPACE << validCoeffs
            << token::SPACE << validInterface;

        for (const label intI : validCoeffs)
        {
            const label start = coeffsOffsets[intI];
            const label size = coeffsOffsets[intI+1] - start;

            toMaster
                << SubField<scalar>(interfaceBouCoeffs, size, start)
                << SubField<scalar>(interfaceIntCoeffs, size, start);
        }
        for (const label intI : validInterface)
        {
            const auto& interface = refCast<const GAMGInterfaceField>
            (
                interfaces[intI]
            );

            toMaster << interface.type();
            interface.write(toMaster);
        }
    }

    // Wait for finish
    pBufs.finishedGathers();

    // Consume
    if (UPstream::master(comm))
    {
        const label nProcs = UPstream::nProcs(comm);

        const lduMesh& destMesh = agglomeration_.meshLevel(destLevel);
        lduInterfacePtrsList destInterfaces = destMesh.interfaces();

        // Master.
        otherMats.resize(nProcs-1);
        otherBouCoeffs.resize(nProcs-1);
        otherIntCoeffs.resize(nProcs-1);
        otherInterfaces.resize(nProcs-1);

        for (const int proci : UPstream::subProcs(comm))
        {
            const label otherI = proci-1;

            UIPstream fromProc(proci, pBufs);

            otherMats.set(otherI, new lduMatrix(destMesh, fromProc));

            // Receive bit-sets of valid interfaceCoeffs/interfaces
            const bitSet validCoeffs(fromProc);
            const bitSet validInterface(fromProc);

            otherBouCoeffs.emplace_set(otherI, validCoeffs.size());
            otherIntCoeffs.emplace_set(otherI, validCoeffs.size());
            otherInterfaces.emplace_set(otherI, validInterface.size());

            // Receive individual interface contributions
            for (const label intI : validCoeffs)
            {
                otherBouCoeffs[otherI].emplace_set(intI, fromProc);
                otherIntCoeffs[otherI].emplace_set(intI, fromProc);
            }

            // Receive individual interface contributions
            for (const label intI : validInterface)
            {
                const word coupleType(fromProc);

                const label allIntI = boundaryMap[proci][intI];

                otherInterfaces[otherI].set
                (
                    intI,
                    GAMGInterfaceField::New
                    (
                        coupleType,
                        refCast<const GAMGInterface>
                        (
                            destInterfaces[allIntI]
                        ),
                        fromProc
                    ).release()
                );
            }
        }
    }
}


void Foam::GAMGSolver::procAgglomerateMatrix
(
    // Agglomeration information
    const labelList& procAgglomMap,
    const List<label>& agglomProcIDs,

    const label levelI,

    // Resulting matrix
    autoPtr<lduMatrix>& allMatrixPtr,
    autoPtr<scalarField>& allInterfaceBouCoeffs,
    autoPtr<scalarField>& allInterfaceIntCoeffs,
    autoPtr<labelList>& allCoeffsOffsets,
    PtrList<lduInterfaceField>& allPrimitiveInterfaces,
    lduInterfaceFieldPtrsList& allInterfaces
) const
{
    const lduMatrix& coarsestMatrix = matrixLevels_[levelI];
    const lduInterfaceFieldPtrsList& coarsestInterfaces =
        interfaceLevels_[levelI];

    // Non-master processors retain empty flat storage.
    allInterfaceBouCoeffs.reset(new scalarField(0));
    allInterfaceIntCoeffs.reset(new scalarField(0));
    allCoeffsOffsets.reset(new labelList(1, 0));

    // Use the existing flat storage directly for gatherMatrices.
    const scalarField& bouCoeffsData = interfaceLevelsBouCoeffs_[levelI];
    const scalarField& intCoeffsData = interfaceLevelsIntCoeffs_[levelI];
    const labelList& coeffsOffsets = interfaceLevelsOffsets_[levelI];

    // Communicator containing all processors to combine (=agglomProcIDs).
    // Result will be on master of communicator.
    const label agglomComm = agglomeration_.agglomCommunicator(levelI+1);


    // Gather all matrix coefficients onto agglomProcIDs[0]
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<lduMatrix> otherMats;
    PtrList<FieldField<Field, scalar>> otherBouCoeffs;
    PtrList<FieldField<Field, scalar>> otherIntCoeffs;
    PtrList<PtrList<lduInterfaceField>> otherInterfaces;
    gatherMatrices
    (
        levelI+1,                       // allMesh level (only on master)
        agglomComm,

        coarsestMatrix,                 // master before gathering
        bouCoeffsData,                  // master before gathering
        intCoeffsData,                  // master before gathering
        coeffsOffsets,                  // master before gathering
        coarsestInterfaces,             // master before gathering

        otherMats,                      // slave info
        otherBouCoeffs,                 // slave info
        otherIntCoeffs,                 // slave info
        otherInterfaces                 // slave info
    );


    if (UPstream::master(agglomComm))
    {
        // Agglomerate all matrix
        // ~~~~~~~~~~~~~~~~~~~~~~

        const lduMesh& allMesh = agglomeration_.meshLevel(levelI+1);
        const labelList& cellOffsets = agglomeration_.cellOffsets(levelI+1);
        const labelListList& faceMap = agglomeration_.faceMap(levelI+1);
        const labelListList& boundaryMap = agglomeration_.boundaryMap(levelI+1);
        const labelListListList& boundaryFaceMap =
            agglomeration_.boundaryFaceMap(levelI+1);

        allMatrixPtr.reset(new lduMatrix(allMesh));
        lduMatrix& allMatrix = allMatrixPtr();

        if (coarsestMatrix.hasDiag())
        {
            scalarField& allDiag = allMatrix.diag();

            SubList<scalar>
            (
                allDiag,
                coarsestMatrix.diag().size()
            ) = coarsestMatrix.diag();

            forAll(otherMats, i)
            {
                SubList<scalar>
                (
                    allDiag,
                    otherMats[i].diag().size(),
                    cellOffsets[i+1]
                ) = otherMats[i].diag();
            }
        }
        if (coarsestMatrix.hasLower())
        {
            scalarField& allLower = allMatrix.lower();
            UIndirectList<scalar>
            (
                allLower,
                faceMap[0]
            ) = coarsestMatrix.lower();
            forAll(otherMats, i)
            {
                UIndirectList<scalar>
                (
                    allLower,
                    faceMap[i+1]
                ) = otherMats[i].lower();
            }
        }
        if (coarsestMatrix.hasUpper())
        {
            scalarField& allUpper = allMatrix.upper();
            UIndirectList<scalar>
            (
                allUpper,
                faceMap[0]
            ) = coarsestMatrix.upper();
            forAll(otherMats, i)
            {
                UIndirectList<scalar>
                (
                    allUpper,
                    faceMap[i+1]
                ) = otherMats[i].upper();
            }
        }


        // Agglomerate interface fields and coefficients
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        lduInterfacePtrsList allMeshInterfaces = allMesh.interfaces();

        allCoeffsOffsets.reset
        (
            new labelList(allMeshInterfaces.size() + 1, 0)
        );

        forAll(allMeshInterfaces, intI)
        {
            (*allCoeffsOffsets)[intI+1] = (*allCoeffsOffsets)[intI]
                + allMeshInterfaces[intI].faceCells().size();
        }

        const label totalSize = allCoeffsOffsets->last();

        allInterfaceBouCoeffs.reset(new scalarField(totalSize, Zero));
        allInterfaceIntCoeffs.reset(new scalarField(totalSize, Zero));
        allPrimitiveInterfaces.setSize(allMeshInterfaces.size());
        allInterfaces.setSize(allMeshInterfaces.size());

        UPtrList<lduInterfaceField> otherFlds(0);

        forAll(boundaryMap, proci)
        {
            const labelList& bMap = boundaryMap[proci];
            forAll(bMap, procIntI)
            {
                label allIntI = bMap[procIntI];

                if (allIntI != -1)
                {
                    // So this boundary has been preserved. Copy
                    // data across.

                    if (!allInterfaces.set(allIntI))
                    {
                        const GAMGInterface& intf = refCast<const GAMGInterface>
                        (
                            allMeshInterfaces[allIntI]
                        );

                        if (proci == 0)
                        {
                            // Clone my local interfaceField. Since it is from
                            // this processor it will already exist, even if it
                            // is a processor one.

                            const auto& ffld =
                            refCast<const GAMGInterfaceField>
                            (
                                coarsestInterfaces[procIntI]
                            );

                            allPrimitiveInterfaces.set
                            (
                                allIntI,
                                ffld.clone
                                (
                                    intf,
                                    otherFlds
                                ).release()
                            );
                        }
                        else
                        {
                            // Recreate a remote interfaceField
                            if (otherInterfaces[proci-1].set(procIntI))
                            {
                                const auto& ffld =
                                refCast<const GAMGInterfaceField>
                                (
                                    otherInterfaces[proci-1][procIntI]
                                );

                                allPrimitiveInterfaces.set
                                (
                                    allIntI,
                                    ffld.clone
                                    (
                                        intf,
                                        otherFlds
                                    ).release()
                                );
                            }
                            else
                            {
                                // Recreate a default interfaceField with
                                // sensible defaults.
                                // Should not occur since all unmerged
                                // processor interfaces get transferred.
                                allPrimitiveInterfaces.set
                                (
                                    allIntI,
                                    GAMGInterfaceField::New
                                    (
                                        intf,
                                        false,  //doTransform,
                                        0       //rank
                                    ).ptr()
                                );
                            }
                        }
                        allInterfaces.set
                        (
                            allIntI,
                            allPrimitiveInterfaces.get(allIntI)
                        );
                    }


                }

                const bool hasCoeffs =
                (
                    proci == 0
                  ? coeffsOffsets[procIntI+1] > coeffsOffsets[procIntI]
                                    : otherBouCoeffs[proci-1].set(procIntI) != nullptr
                );

                if (!hasCoeffs)
                {
                    continue;
                }

                auto mapCoeffs =
                [&](const UList<scalar>& procBou,
                    const UList<scalar>& procInt)
                {
                    const labelList& map =
                        boundaryFaceMap[proci][procIntI];

                    if (allIntI != -1)
                    {
                        SubField<scalar> allBou
                        (
                            *allInterfaceBouCoeffs,
                            (*allCoeffsOffsets)[allIntI+1]
                              - (*allCoeffsOffsets)[allIntI],
                            (*allCoeffsOffsets)[allIntI]
                        );
                        SubField<scalar> allInt
                        (
                            *allInterfaceIntCoeffs,
                            (*allCoeffsOffsets)[allIntI+1]
                              - (*allCoeffsOffsets)[allIntI],
                            (*allCoeffsOffsets)[allIntI]
                        );

                        forAll(map, i)
                        {
                            label allFacei = map[i];
                            if (allFacei < 0)
                            {
                                FatalErrorInFunction
                                    << "problem." << abort(FatalError);
                            }
                            allBou[allFacei] = procBou[i];
                            allInt[allFacei] = procInt[i];
                        }
                    }
                    else
                    {
                        forAll(map, i)
                        {
                            if (map[i] >= 0)
                            {
                                label allFacei = map[i];

                                if (coarsestMatrix.hasUpper())
                                {
                                    allMatrix.upper()[allFacei] = -procBou[i];
                                }
                                if (coarsestMatrix.hasLower())
                                {
                                    allMatrix.lower()[allFacei] = -procInt[i];
                                }
                            }
                            else
                            {
                                label allFacei = -map[i]-1;

                                if (coarsestMatrix.hasUpper())
                                {
                                    allMatrix.upper()[allFacei] = -procInt[i];
                                }
                                if (coarsestMatrix.hasLower())
                                {
                                    allMatrix.lower()[allFacei] = -procBou[i];
                                }
                            }
                        }
                    }
                };

                if (proci == 0)
                {
                    const label start = coeffsOffsets[procIntI];
                    const label size = coeffsOffsets[procIntI+1] - start;
                    SubField<scalar> procBou(bouCoeffsData, size, start);
                    SubField<scalar> procInt(intCoeffsData, size, start);
                    mapCoeffs(procBou, procInt);
                }
                else
                {
                    mapCoeffs
                    (
                        otherBouCoeffs[proci-1][procIntI],
                        otherIntCoeffs[proci-1][procIntI]
                    );
                }
            }
        }


        //Pout<< "** Assembled allMatrix:" << allMatrix.info() << endl;
        //
        //forAll(allInterfaces, intI)
        //{
        //    if (allInterfaces.set(intI))
        //    {
        //        Pout<< "    patch:" << intI
        //            << " type:" << allInterfaces[intI].type()
        //            << " size:"
        //            << allInterfaces[intI].interface().faceCells().size()
        //            << endl;
        //
        //        const scalarField& bouCoeffs = allInterfaceBouCoeffs[intI];
        //        const scalarField& intCoeffs = allInterfaceIntCoeffs[intI];
        //        forAll(bouCoeffs, facei)
        //        {
        //            Pout<< "        " << facei
        //                << "\tbou:" << bouCoeffs[facei]
        //                << "\tint:" << intCoeffs[facei]
        //                << endl;
        //        }
        //    }
        //}
    }
}


void Foam::GAMGSolver::procAgglomerateMatrix
(
    const labelList& procAgglomMap,
    const List<label>& agglomProcIDs,

    const label levelI
)
{
    autoPtr<lduMatrix> allMatrixPtr;
    autoPtr<scalarField> allInterfaceBouCoeffs;
    autoPtr<scalarField> allInterfaceIntCoeffs;
    autoPtr<labelList> allCoeffsOffsets;
    autoPtr<PtrList<lduInterfaceField>> allPrimitiveInterfaces
    (
        new PtrList<lduInterfaceField>(0)
    );
    autoPtr<lduInterfaceFieldPtrsList> allInterfaces
    (
        new lduInterfaceFieldPtrsList(0)
    );

    procAgglomerateMatrix
    (
        // Agglomeration information
        procAgglomMap,
        agglomProcIDs,

        levelI,

        // Resulting matrix
        allMatrixPtr,
        allInterfaceBouCoeffs,
        allInterfaceIntCoeffs,
        allCoeffsOffsets,
        allPrimitiveInterfaces(),
        allInterfaces()
    );

    matrixLevels_.set(levelI, allMatrixPtr);

    interfaceLevelsBouCoeffs_.set
    (
        levelI,
        allInterfaceBouCoeffs.release()
    );
    interfaceLevelsIntCoeffs_.set
    (
        levelI,
        allInterfaceIntCoeffs.release()
    );
    interfaceLevelsOffsets_.set
    (
        levelI,
        allCoeffsOffsets.release()
    );
    primitiveInterfaceLevels_.set(levelI, allPrimitiveInterfaces);
    interfaceLevels_.set(levelI, allInterfaces);
}


// ************************************************************************* //
