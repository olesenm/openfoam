/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023-2026 OpenCFD Ltd.
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

#include "pairGAMGAgglomeration.H"
#include "lduAddressing.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::pairGAMGAgglomeration::cellCells
(
    const lduAddressing& fineMatrixAddressing,
    const label nCoarseCells,
    const labelList& coarseCellMap,
    labelList& nbrCells,
    labelList& cellOffsets
)
{
    // Construct cell-cell addressing on coarse level

    const labelUList& upperAddr = fineMatrixAddressing.upperAddr();
    const labelUList& lowerAddr = fineMatrixAddressing.lowerAddr();

    // Pass 1: count
    labelList nNbrs(nCoarseCells, 0);

    forAll(upperAddr, facei)
    {
        const label coarseUpper = coarseCellMap[upperAddr[facei]];
        const label coarseLower = coarseCellMap[lowerAddr[facei]];

        if (coarseUpper != coarseLower)
        {
            nNbrs[coarseUpper]++;
            nNbrs[coarseLower]++;
        }
    }

    // Size
    cellOffsets.setSize(nCoarseCells+1);
    cellOffsets[0] = 0;
    forAll(nNbrs, celli)
    {
        cellOffsets[celli+1] = cellOffsets[celli] + nNbrs[celli];
    }
    nbrCells.setSize(cellOffsets.last());

    // Pass2: fill
    nNbrs = 0;

    forAll(upperAddr, facei)
    {
        const label coarseUpper = coarseCellMap[upperAddr[facei]];
        const label coarseLower = coarseCellMap[lowerAddr[facei]];

        if (coarseUpper != coarseLower)
        {
            {
                const label start = cellOffsets[coarseUpper];
                label& nUsed = nNbrs[coarseUpper];
                // SubList<label> used(nbrCells, nUsed, start);

                // if (used.find(coarseLower) == -1)
                {
                    nbrCells[start + nUsed++] = coarseLower;
                }
            }
            {
                const label start = cellOffsets[coarseLower];
                label& nUsed = nNbrs[coarseLower];
                // SubList<label> used(nbrCells, nUsed, start);

                // if (used.find(coarseUpper) == -1)
                {
                    nbrCells[start + nUsed++] = coarseUpper;
                }
            }
        }
    }

    // Filter duplicates
    forAll(nNbrs, celli)
    {
        const label start = cellOffsets[celli];
        const label end = cellOffsets[celli+1];
        const label size = end - start;

        label nUsed = 1;

        // Filter duplicates.
        for (label i = 1; i < size; i++)
        {
            const label nbr = nbrCells[start + i];

            SubList<label> used(nbrCells, i-1, start);
            if (used.find(nbr) == -1)
            {
                nbrCells[start + nUsed++] = nbr;
            }
        }

        // Update offset for next level
        cellOffsets[celli+1] = cellOffsets[celli] + nUsed;
    }
    nbrCells.setSize(cellOffsets.last());
}


void Foam::pairGAMGAgglomeration::agglomerate
(
    const label nCellsInCoarsestLevel,
    const label startLevel,
    const scalarField& startFaceWeights,
    const bool doProcessorAgglomerate
)
{
    if (nCells_.size() < maxLevels_)
    {
        // See compactLevels. Make space if not enough
        nCells_.resize(maxLevels_);
        restrictAddressing_.resize(maxLevels_);
        nFaces_.resize(maxLevels_);
        faceRestrictAddressing_.resize(maxLevels_);
        faceFlipMap_.resize(maxLevels_);
        nPatchFaces_.resize(maxLevels_);
        patchFaceRestrictAddressing_.resize(maxLevels_);
        meshLevels_.resize(maxLevels_);
        // Have procCommunicator_ always, even if not procAgglomerating.
        // Use value -1 to indicate nothing is proc-agglomerated
        procCommunicator_.resize(maxLevels_ + 1, -1);
        if (processorAgglomerate())
        {
            procAgglomMap_.resize(maxLevels_);
            agglomProcIDs_.resize(maxLevels_);
            procCommunicator_.resize(maxLevels_);
            procCellOffsets_.resize(maxLevels_);
            procFaceMap_.resize(maxLevels_);
            procBoundaryMap_.resize(maxLevels_);
            procBoundaryFaceMap_.resize(maxLevels_);
        }
    }


    // Start geometric agglomeration from the given faceWeights
    scalarField faceWeights = startFaceWeights;

    // Agglomerate until the required number of cells in the coarsest level
    // is reached

    label nPairLevels = 0;
    label nCreatedLevels = startLevel;
    bool forward = true;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        if (!hasMeshLevel(nCreatedLevels))
        {
            FatalErrorInFunction<< "No mesh at nCreatedLevels:"
                << nCreatedLevels
                << exit(FatalError);
        }

        const auto& fineMesh = meshLevel(nCreatedLevels);


        label nCoarseCells = -1;

        tmp<labelField> finalAgglomPtr = agglomerate
        (
            nCoarseCells,
            fineMesh.lduAddr(),
            faceWeights,
            renumber_,
            forward
        );

        if
        (
            continueAgglomerating
            (
                nCellsInCoarsestLevel,
                finalAgglomPtr().size(),
                nCoarseCells,
                fineMesh.comm()
            )
        )
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        }
        else
        {
            break;
        }

        // Create coarse mesh
        agglomerateLduAddressing(nCreatedLevels);

        // Agglomerate the faceWeights field for the next level
        {
            scalarField aggFaceWeights
            (
                meshLevels_[nCreatedLevels].upperAddr().size(),
                0.0
            );

            restrictFaceField
            (
                aggFaceWeights,
                faceWeights,
                nCreatedLevels
            );

            faceWeights = std::move(aggFaceWeights);
        }

        if (nPairLevels % mergeLevels_)
        {
            combineLevels(nCreatedLevels);
        }
        else
        {
            nCreatedLevels++;
        }

        nPairLevels++;
    }

    // Shrink the storage of the levels to those created
    compactLevels(nCreatedLevels, doProcessorAgglomerate);


    // Renumber coarse cells at each level (CutHill-McKee) to improve
    // the band profile of the agglomerated matrices
    if (renumber_)
    {
        for (label leveli = 0; leveli < nCreatedLevels; leveli++)
        {
            // Build coarse cell-cell adjacency from the existing coarse mesh
            const lduAddressing& coarseAddr = meshLevels_[leveli].lduAddr();
            const labelUList& upper = coarseAddr.upperAddr();
            const labelUList& lower = coarseAddr.lowerAddr();
            const label nCoarse = nCells_[leveli];

            labelList nNbrs(nCoarse, 0);
            forAll(upper, facei)
            {
                nNbrs[upper[facei]]++;
                nNbrs[lower[facei]]++;
            }

            labelList cellOffsets(nCoarse + 1);
            cellOffsets[0] = 0;
            forAll(nNbrs, celli)
            {
                cellOffsets[celli + 1] = cellOffsets[celli] + nNbrs[celli];
            }

            labelList nbrCells(cellOffsets.last());
            nNbrs = 0;
            forAll(upper, facei)
            {
                nbrCells[cellOffsets[upper[facei]] + nNbrs[upper[facei]]++] =
                    lower[facei];
                nbrCells[cellOffsets[lower[facei]] + nNbrs[lower[facei]]++] =
                    upper[facei];
            }

            // Compute CutHill-McKee ordering
            const labelList newToOld
            (
                meshTools::bandCompression(nbrCells, cellOffsets)
            );
            const labelList oldToNew(invert(newToOld.size(), newToOld));

            // Apply new coarse-cell numbering to the restrict map at this level
            labelField& restrict = restrictAddressing_[leveli];
            restrict = UIndirectList<label>(oldToNew, restrict)();

            // Permute the next level's restrict map so its domain indices
            // match the new coarse-cell ordering of this level
            if (leveli + 1 < nCreatedLevels)
            {
                labelField& nextRestrict = restrictAddressing_[leveli + 1];
                nextRestrict = UIndirectList<label>(nextRestrict, newToOld)();
            }

            // Rebuild the coarse mesh with the new cell ordering
            agglomerateLduAddressing(leveli);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::pairGAMGAgglomeration::agglomerate
(
    label& nCoarseCells,
    const lduAddressing& fineMatrixAddressing,
    const scalarField& faceWeights,
    const bool renumber
)
{
    bool forward = true;

    return agglomerate
    (
        nCoarseCells,
        fineMatrixAddressing,
        faceWeights,
        renumber,
        forward
    );
}


Foam::tmp<Foam::labelField> Foam::pairGAMGAgglomeration::agglomerate
(
    label& nCoarseCells,
    const lduAddressing& fineMatrixAddressing,
    const scalarField& faceWeights,
    const bool renumber,
    bool& forward
)
{
    const label nFineCells = fineMatrixAddressing.size();

    const labelUList& upperAddr = fineMatrixAddressing.upperAddr();
    const labelUList& lowerAddr = fineMatrixAddressing.lowerAddr();

    // For each cell calculate faces
    labelList cellFaces(upperAddr.size() + lowerAddr.size());
    labelList cellFaceOffsets(nFineCells + 1);

    // memory management
    {
        labelList nNbrs(nFineCells, Zero);

        forAll(upperAddr, facei)
        {
            nNbrs[upperAddr[facei]]++;
        }

        forAll(lowerAddr, facei)
        {
            nNbrs[lowerAddr[facei]]++;
        }

        cellFaceOffsets[0] = 0;
        forAll(nNbrs, celli)
        {
            cellFaceOffsets[celli+1] = cellFaceOffsets[celli] + nNbrs[celli];
        }

        // reset the whole list to use as counter
        nNbrs = 0;

        forAll(upperAddr, facei)
        {
            cellFaces
            [
                cellFaceOffsets[upperAddr[facei]] + nNbrs[upperAddr[facei]]
            ] = facei;

            nNbrs[upperAddr[facei]]++;
        }

        forAll(lowerAddr, facei)
        {
            cellFaces
            [
                cellFaceOffsets[lowerAddr[facei]] + nNbrs[lowerAddr[facei]]
            ] = facei;

            nNbrs[lowerAddr[facei]]++;
        }
    }


    // Determine cell visit order. This will lower the effect of having
    // one cell clustering with a neighbouring one whereas it would be better
    // if that neighbour clustered with another one of its neighbours. So
    // instead make sure to start with the cells with highest face weights
    // first. This causes larger clusters to be formed first, but might
    // leave some smaller clusters at the end.

    labelList visitOrder;
    // if (renumber)
    // {
    //     scalarField maxCellWeight(nFineCells, -GREAT);
    //     forAll(upperAddr, facei)
    //     {
    //         scalar& cWeight = maxCellWeight[upperAddr[facei]];
    //         cWeight = max(cWeight, faceWeights[facei]);
    //     }
    //     forAll(lowerAddr, facei)
    //     {
    //         scalar& cWeight = maxCellWeight[lowerAddr[facei]];
    //         cWeight = max(cWeight, faceWeights[facei]);
    //     }

    //     sortedOrder
    //     (
    //         maxCellWeight,
    //         visitOrder,
    //         typename UList<scalar>::greater(maxCellWeight)
    //     );
    // }



    // Go through the faces and create clusters

    auto tcoarseCellMap = tmp<labelField>::New(nFineCells, -1);
    auto& coarseCellMap = tcoarseCellMap.ref();

    // Small tolerance to account for faces potentially having slightly
    // different truncation error in their weights from run to run
    // (e.g. due to offloading). If all the largest faces per cell are
    // within this tolerance use the first one. This guarantees repeatability.
    // Disabled on non-offload situations for now since makes comparison
    // to previous versions harder. TBD.
    #ifdef _OPENMP
    const scalar tol = 1E-10;
    #else
    const scalar tol = 0;
    #endif

    nCoarseCells = 0;
    for (label cellfi=0; cellfi<nFineCells; cellfi++)
    {
        // Change cell ordering depending on direction for this level
        const label celli =
        (
            visitOrder.size()
          ? visitOrder[cellfi]
          : forward ? cellfi : nFineCells - cellfi - 1
        );

        if (coarseCellMap[celli] < 0)
        {
            label matchFaceNo = -1;
            scalar maxFaceWeight = -GREAT;

            // check faces to find ungrouped neighbour with largest face weight
            for
            (
                label faceOs=cellFaceOffsets[celli];
                faceOs<cellFaceOffsets[celli+1];
                faceOs++
            )
            {
                label facei = cellFaces[faceOs];

                // I don't know whether the current cell is owner or neighbour.
                // Therefore I'll check both sides
                if
                (
                    coarseCellMap[upperAddr[facei]] < 0
                 && coarseCellMap[lowerAddr[facei]] < 0
                 && faceWeights[facei] > maxFaceWeight*(1.0 + tol)
                )
                {
                    // Match found. Pick up all the necessary data
                    matchFaceNo = facei;
                    maxFaceWeight = faceWeights[facei];
                }
            }

            if (matchFaceNo >= 0)
            {
                // Make a new group
                coarseCellMap[upperAddr[matchFaceNo]] = nCoarseCells;
                coarseCellMap[lowerAddr[matchFaceNo]] = nCoarseCells;
                nCoarseCells++;
            }
            else
            {
                // No match. Find the best neighbouring cluster and
                // put the cell there
                label clusterMatchFaceNo = -1;
                scalar clusterMaxFaceCoeff = -GREAT;

                for
                (
                    label faceOs=cellFaceOffsets[celli];
                    faceOs<cellFaceOffsets[celli+1];
                    faceOs++
                )
                {
                    label facei = cellFaces[faceOs];

                    if (faceWeights[facei] > clusterMaxFaceCoeff*(1.0 + tol))
                    {
                        clusterMatchFaceNo = facei;
                        clusterMaxFaceCoeff = faceWeights[facei];
                    }
                }

                if (clusterMatchFaceNo >= 0)
                {
                    // Add the cell to the best cluster
                    coarseCellMap[celli] = max
                    (
                        coarseCellMap[upperAddr[clusterMatchFaceNo]],
                        coarseCellMap[lowerAddr[clusterMatchFaceNo]]
                    );
                }
            }
        }
    }

    // Check that all cells are part of clusters,
    // if not create single-cell "clusters" for each
    for (label cellfi=0; cellfi<nFineCells; cellfi++)
    {
        // Change cell ordering depending on direction for this level
        const label celli =
        (
            visitOrder.size()
          ? visitOrder[cellfi]
          : forward ? cellfi : nFineCells - cellfi - 1
        );

        if (coarseCellMap[celli] < 0)
        {
            coarseCellMap[celli] = nCoarseCells;
            nCoarseCells++;
        }
    }

    if (!forward)
    {
        forAll(coarseCellMap, celli)
        {
            coarseCellMap[celli] = nCoarseCells - 1 - coarseCellMap[celli];
        }
    }

    // Reverse the map ordering for the next level
    // to improve the next level of agglomeration
    forward = !forward;

    return tcoarseCellMap;
}


// ************************************************************************* //
