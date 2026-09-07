/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019-2026 OpenCFD Ltd.
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

#include "procFacesGAMGProcAgglomeration.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"
#include "lduMesh.H"
#include "processorLduInterface.H"
#include "processorGAMGInterface.H"
#include "pairGAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(procFacesGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        procFacesGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- copy from GAMGAgglomerateLduAddressing but adapted for lduPrimitiveMesh
//  and to create processorGAMGInterfaces instead of GAMGInterfaces.
Foam::autoPtr<Foam::lduPrimitiveMesh>
Foam::procFacesGAMGProcAgglomeration::agglomerateLduAddressing
(
    const lduPrimitiveMesh& fineMesh,
    const labelUList& restrictMap,
    label& nCoarseFaces,
    labelList& faceRestrictAddr,
    boolList& faceFlipMap
)
{
    const lduAddressing& fineMeshAddr = fineMesh.lduAddr();

    const labelUList& upperAddr = fineMeshAddr.upperAddr();
    const labelUList& lowerAddr = fineMeshAddr.lowerAddr();

    const label nFineFaces = upperAddr.size();

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
    const label nCoarseCells = (restrictMap.size() ? max(restrictMap) + 1 : 0);

    // Storage for coarse cell neighbours and coefficients

    // Guess initial maximum number of neighbours in coarse cell
    label maxNnbrs = 10;

    // Number of faces for each coarse-cell
    labelList cCellnFaces(nCoarseCells, Zero);

    // Setup initial packed storage for coarse-cell faces
    labelList cCellFaces(maxNnbrs*nCoarseCells);

    // Create face-restriction addressing
    faceRestrictAddr.resize_nocopy(nFineFaces);

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineFaces);

    // Counter for coarse faces
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
    faceFlipMap.setSize(nFineFaces, false);

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
    const auto& fineInterfaces = fineMesh.rawInterfaces();

    const label startOfRequests = UPstream::nRequests();

    // Initialise transfer of restrict addressing on the interface
    // The finest mesh uses patchAddr from the original lduAdressing.
    // the coarser levels create their own adressing for faceCells
    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            fineInterfaces[inti].initInternalFieldTransfer
            (
                Pstream::commsTypes::nonBlocking,
                restrictMap,
                fineMeshAddr.patchAddr(inti)
            );
        }
    }

    UPstream::waitRequests(startOfRequests);


    // Add the coarse level
    autoPtr<lduPrimitiveMesh> coarseMeshPtr
    (
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

            restrictMapInternalField =
                fineInterfaces[inti].interfaceInternalField
                (
                    restrictMap,
                    fineMeshAddr.patchAddr(inti)
                );

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
                    fineInterfaces, // fineMesh.rawInterfaces(),
                    fineInterfaces[inti],
                    restrictMapInternalField(),
                    nbrRestrictMapInternalField(),
                    0,  // fineLevelIndex,
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
                    fineMesh.rawInterfaces(),
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
        }
    }

    coarseMeshPtr().addInterfaces
    (
        coarseInterfaces,
        lduPrimitiveMesh::nonBlockingSchedule<processorGAMGInterface>
        (
            coarseInterfaces
        )
    );

    if (debug & 2)
    {
        const auto& coarseAddr = coarseMeshPtr().lduAddr();

        Pout<< "GAMGAgglomeration :"
            // << " agglomerated level " << fineLevelIndex
            << " from nCells:" << fineMeshAddr.size()
            << " nFaces:" << upperAddr.size()
            << " to nCells:" << nCoarseCells
            << " nFaces:" << nCoarseFaces << nl
            << "    lower:" << flatOutput(coarseAddr.lowerAddr()) << nl
            << "    upper:" << flatOutput(coarseAddr.upperAddr()) << nl
            << endl;
    }

    return coarseMeshPtr;
}


Foam::autoPtr<Foam::lduPrimitiveMesh>
Foam::procFacesGAMGProcAgglomeration::singleCellMesh
(
    const label singleCellMeshComm,
    const lduMesh& mesh,
    scalarField& faceWeights
) const
{
    // Return processor-connectivity as an ldu mesh with the weights
    // according to the number of faces on the processor interfaces.
    // Note: valid on master only
    List<Map<label>> procFaces(UPstream::nProcs(mesh.comm()));
    Map<label>& myNeighbours = procFaces[UPstream::myProcNo(mesh.comm())];

    {
        const lduInterfacePtrsList interfaces(mesh.interfaces());
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const processorLduInterface& pp =
                    refCast<const processorLduInterface>
                    (
                        interfaces[intI]
                    );
                label size = interfaces[intI].faceCells().size();
                myNeighbours.insert(pp.neighbProcNo(), size);
            }
        }
    }

    Pstream::allGatherList(procFaces, Pstream::msgType(), mesh.comm());

    autoPtr<lduPrimitiveMesh> singleCellMeshPtr;

    if (Pstream::master(mesh.comm()))
    {
        // I am master
        label nCells = Pstream::nProcs(mesh.comm());

        DynamicList<label> l(3*nCells);
        DynamicList<label> u(3*nCells);
        DynamicList<scalar> weight(3*nCells);

        DynamicList<label> nbrs;
        DynamicList<scalar> weights;

        forAll(procFaces, proci)
        {
            const Map<label>& neighbours = procFaces[proci];

            // Add all the higher processors
            nbrs.clear();
            weights.clear();
            forAllConstIters(neighbours, iter)
            {
                if (iter.key() > proci)
                {
                    nbrs.append(iter.key());
                    weights.append(iter());
                }
                sort(nbrs);
                forAll(nbrs, i)
                {
                    l.append(proci);
                    u.append(nbrs[i]);
                    weight.append(weights[i]);
                }
            }
        }

        faceWeights.transfer(weight);

        PtrList<const lduInterface> primitiveInterfaces(0);
        const lduSchedule ps(0);

        singleCellMeshPtr.reset
        (
            new lduPrimitiveMesh
            (
                nCells,
                l,
                u,
                primitiveInterfaces,
                ps,
                singleCellMeshComm
            )
        );
    }
    return singleCellMeshPtr;
}


Foam::tmp<Foam::labelField>
Foam::procFacesGAMGProcAgglomeration::processorAgglomeration
(
    const lduMesh& mesh,
    const label nCoarsest
) const
{
    UPstream::communicator singleCellMeshComm
    (
        mesh.comm(),
        labelRange(1)  // Processor 0 only
    );

    // Construct ldu mesh from processor connectivity
    scalarField faceWeights;
    autoPtr<lduPrimitiveMesh> singleCellMeshPtr
    (
        singleCellMesh
        (
            singleCellMeshComm.comm(),
            mesh,
            faceWeights
        )
    );

    auto tfineToCoarse = tmp<labelField>::New();
    auto& fineToCoarse = tfineToCoarse.ref();

    if (singleCellMeshPtr)
    {
        // On master call the pair-wise agglomerator
        // Note: no need for renumbering since the cells are processors
        // and there are no algorithms that rely on the ordering of
        // the processors.
        label nCoarseProcs;
        bool forward = true;
        fineToCoarse = pairGAMGAgglomeration::agglomerate
        (
            nCoarseProcs,
            singleCellMeshPtr(),
            faceWeights,
            false,  // no need for renumber
            forward
        );

        // Keep on agglomerating until we've reached the limit. (nCoarsest is
        // labelMax it no iteration needed)
        labelList localFineToCoarse(fineToCoarse);
        while (nCoarseProcs > nCoarsest)
        {
            label nCoarseFaces;
            labelList faceRestrictAddr;
            boolList faceFlipMap;
            singleCellMeshPtr = agglomerateLduAddressing
            (
                singleCellMeshPtr(),
                localFineToCoarse,
                nCoarseFaces,
                faceRestrictAddr,
                faceFlipMap
            );

            // Accumulate face-weights
            scalarField coarseFaceWeights(nCoarseFaces, 0.0);
            forAll(faceWeights, facei)
            {
                label coarseFacei = faceRestrictAddr[facei];
                if (coarseFacei >= 0)
                {
                    coarseFaceWeights[coarseFacei] += faceWeights[facei];
                }
            }

            // Use new face weights to agglomerate cells.
            localFineToCoarse = pairGAMGAgglomeration::agglomerate
            (
                nCoarseProcs,
                singleCellMeshPtr(),
                coarseFaceWeights,
                false,  // no need for renumber
                forward
            );
            faceWeights.transfer(coarseFaceWeights);

            // Update ultimate agglomeration
            forAll(fineToCoarse, celli)
            {
                fineToCoarse[celli] = localFineToCoarse[fineToCoarse[celli]];
            }
        }

        labelList coarseToMaster(nCoarseProcs, labelMax);
        forAll(fineToCoarse, celli)
        {
            label coarseI = fineToCoarse[celli];
            coarseToMaster[coarseI] = min(coarseToMaster[coarseI], celli);
        }

        // Sort according to master and redo restriction
        labelList newToOld(sortedOrder(coarseToMaster));
        labelList oldToNew(invert(newToOld.size(), newToOld));

        fineToCoarse = labelUIndList(oldToNew, fineToCoarse)();
    }

    Pstream::broadcast(fineToCoarse, mesh.comm());

    return tfineToCoarse;
}


bool Foam::procFacesGAMGProcAgglomeration::doProcessorAgglomeration
(
    const lduMesh& mesh
) const
{
    // Check the need for further agglomeration on any processors
    const label nCells = mesh.lduAddr().size();
    bool doAgg =
    (
        nMasters_ < labelMax
      ? UPstream::nProcs(mesh.comm()) > nMasters_
      : nCells < nAgglomeratingCells_
    );
    UPstream::reduceOr(doAgg, mesh.comm());
    return doAgg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::procFacesGAMGProcAgglomeration::procFacesGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    nAgglomeratingCells_(-1),
    nMasters_(labelMax),
    nCellsInMasterLevel_
    (
        controlDict.getOrDefault<label>("nCellsInMasterLevel", -1)
    )
{
    // See masterCoarsest. Here it is easier to work from nMasters
    // instead of nProcessorsPerMaster.
    if
    (
        controlDict.readIfPresent<label>
        (
            "nMasters",
            nMasters_,
            keyType::LITERAL
        )
    )
    {
        if (controlDict.found("nAgglomeratingCells"))
        {
            FatalIOErrorInFunction(controlDict)
                << "Cannot use both 'nMasters' and"
                << " 'nAgglomeratingCells'"
                << exit(FatalIOError);
        }
    }
    else if
    (
        label nProcessorsPerMaster;
        controlDict.readIfPresent<label>
        (
            "nProcessorsPerMaster", nProcessorsPerMaster, keyType::LITERAL
        )
    )
    {
        nMasters_ =
            UPstream::nProcs(agglom.mesh().comm())
          / nProcessorsPerMaster;
    }
    else
    {
        nAgglomeratingCells_ = controlDict.get<label>("nAgglomeratingCells");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::procFacesGAMGProcAgglomeration::agglomerate
(
    const label fineLevelIndex,
    const label nCoarsest
)
{
    if (agglom_.hasMeshLevel(fineLevelIndex))
    {
        // Get the fine mesh
        const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);
        const label levelComm = levelMesh.comm();
        label nProcs = UPstream::nProcs(levelComm);

        if (nProcs > 1 && doProcessorAgglomeration(levelMesh))
        {
            tmp<labelField> tprocAgglomMap
            (
                processorAgglomeration(levelMesh, nCoarsest)
            );
            const labelField& procAgglomMap = tprocAgglomMap();

            // Master processor
            labelList masterProcs;
            // Local processors that agglomerate. agglomProcIDs[0] is in
            // masterProc.
            List<label> agglomProcIDs;
            GAMGAgglomeration::calculateRegionMaster
            (
                levelComm,
                procAgglomMap,
                masterProcs,
                agglomProcIDs
            );

            if (debug)
            {
                if (masterProcs.size())
                {
                    labelListList masterToProcs
                    (
                        invertOneToMany
                        (
                            masterProcs.size(),
                            procAgglomMap
                        )
                    );
                    Info<< typeName << " : agglomerating" << nl
                        << "\tmaster\tnProcs\tprocIDs" << endl;
                    for (const auto& p : masterToProcs)
                    {
                        Info<< '\t' << p[0]
                            << "\t\t" << p.size()
                            << "\t\t"
                            << flatOutput(SubList<label>(p, p.size()-1, 1))
                            << endl;
                    }
                }
            }

            // Communicator for the processor-agglomerated matrix
            comms_.push_back
            (
                UPstream::newCommunicator
                (
                    levelComm,
                    masterProcs
                )
            );

            // Use processor agglomeration maps to do the actual
            // collecting.
            GAMGProcAgglomeration::agglomerate
            (
                fineLevelIndex,
                procAgglomMap,
                masterProcs,
                agglomProcIDs,
                comms_.back()
            );
        }
    }
    return true;
}


bool Foam::procFacesGAMGProcAgglomeration::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        if (nMasters_ < labelMax)
        {
            // Do coarsest level only (see masterCoarsestGAMGProcAgglomeration)
            const label fineLevelIndex = agglom_.size()-1;
            // Agglomerate multiple times until the number of master processors
            // is reached
            agglomerate(fineLevelIndex, nMasters_);
        }
        else
        {
            // Apply processor agglomeration on all levels except the finest
            for
            (
                label fineLevelIndex = 2;
                fineLevelIndex < agglom_.size();
                fineLevelIndex++
            )
            {
                // Agglomerate only once
                agglomerate(fineLevelIndex, labelMax);
            }
        }
    }


    // Optionally restart local agglomeration
    // (from masterCoarsest)
    if (nCellsInMasterLevel_ > 0)
    {
        const label levelI = agglom_.size();
        if (agglom_.hasMeshLevel(levelI))
        {
            const lduMesh& fineMesh = agglom_.meshLevel(levelI);
            const auto& addr = fineMesh.lduAddr();
            const scalarField weights
            (
                addr.lowerAddr().size(),
                1.0
            );
            agglom_.agglomerate
            (
                nCellsInMasterLevel_,
                levelI,
                weights,
                false
            );
        }
    }


    // Print a bit
    if (debug)
    {
        Pout<< nl << "Agglomerated mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    return true;
}


// ************************************************************************* //
