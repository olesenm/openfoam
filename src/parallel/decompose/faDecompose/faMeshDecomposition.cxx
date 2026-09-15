/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2026 OpenCFD Ltd.
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

#include "faMeshDecomposition.H"
#include "faMeshTools.H"
#include "Time.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "Map.H"
#include "globalMeshData.H"
#include "cyclicFaPatch.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::faMeshDecomposition::disallowEdgeEncoding_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMeshDecomposition::distributeFaces()
{
    const word& polyMeshRegionName = faMesh::mesh().name();

    Info<< "\nCalculating distribution of finite-area faces ["
        << polyMesh::regionName(areaName_) << "]" << endl;

    for (label proci = 0; proci < nProcs(); ++proci)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(proci)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );

        IOobject ioFvAddr
        (
            "procAddressing",
            "constant",   // <- facesInstance() ?
            polyMesh::meshSubDir,
            procFvMesh,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        );


        // faceProcAddressing (polyMesh)
        ioFvAddr.resetHeader("faceProcAddressing");
        const labelList fvFaceProcAddressing
        (
            labelIOList::readContents(ioFvAddr)
        );

        labelHashSet fvFaceProcAddrHash;
        fvFaceProcAddrHash.reserve(fvFaceProcAddressing.size());

        // If faMesh's fvPatch is a part of the global face zones, faces of that
        // patch will be present on all processors. Because of that, looping
        // through faceProcAddressing will decompose global faMesh faces to the
        // very last processor regardless of where fvPatch is really decomposed.
        // Since global faces which do not belong to specific processor are
        // located at the end of the faceProcAddressing, cutting it at
        // i = owner.size() will correctly decompose faMesh faces.
        // Vanja Skuric, 2016-04-21
        if (hasGlobalFaceZones_)
        {
            // owner (polyMesh)
            ioFvAddr.resetHeader("owner");

            if
            (
                label ownerSize = labelIOList::readContentsSize(ioFvAddr);
                (ownerSize > 0)
            )
            {
                fvFaceProcAddrHash.insert
                (
                    fvFaceProcAddressing.slice(0, ownerSize)
                );
            }
        }
        else
        {
            fvFaceProcAddrHash.insert(fvFaceProcAddressing);
        }

        forAll(faMesh::faceLabels(), facei)
        {
            // With +1 for lookup in faceMap with flip encoding
            const label index = (faMesh::faceLabels()[facei] + 1);

            if (fvFaceProcAddrHash.contains(index))
            {
                faceToProc_[facei] = proci;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshDecomposition::faMeshDecomposition
(
    const word& areaName,
    const polyMesh& mesh,
    const label nProcessors,
    const dictionary& params
)
:
    faMesh(areaName, mesh),
    areaName_(areaName.empty() ? polyMesh::defaultRegion : areaName),
    nProcs_(nProcessors),
    distributed_(false),
    hasGlobalFaceZones_(false),
    cyclicParallel_(false),
    noEdgeEncoding_(disallowEdgeEncoding_),
    faceToProc_(faMesh::nFaces()),
    procFaceLabels_(nProcs_),
    procMeshEdgesMap_(nProcs_),
    procNInternalEdges_(nProcs_, Foam::zero{}),
    procPatchEdgeLabels_(nProcs_),
    procPatchEdgeLookup_(nProcs_),
    procPointAddressing_(nProcs_),
    procEdgeAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procBoundaryAddressing_(nProcs_),
    procPatchRange_(nProcs_),
    procProcessorPatchRange_(nProcs_),
    procNeighbourProcessors_(nProcs_)
{
    updateParameters(params);
}


Foam::faMeshDecomposition::faMeshDecomposition
(
    const polyMesh& mesh,
    const label nProcessors,
    const dictionary& params
)
:
    faMeshDecomposition
    (
        polyMesh::defaultRegion,
        mesh,
        nProcessors,
        params
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMeshDecomposition::updateParameters
(
    const dictionary& params
)
{
    params.readIfPresent("distributed", distributed_);
    if (params.found("globalFaceZones"))
    {
        hasGlobalFaceZones_ = true;
    }
}


void Foam::faMeshDecomposition::decomposeMesh()
{
    // Decide which face goes to which processor
    distributeFaces();

    const word& polyMeshRegionName = faMesh::mesh().name();

    Info<< "\nDistributing faces to processors ["
        << polyMesh::regionName(areaName_) << "]" << endl;

    labelList nLocalFaces(nProcs_, Foam::zero{});

    // Pass 1: determine local sizes, sanity check

    forAll(faceToProc_, facei)
    {
        const label proci = faceToProc_[facei];

        if (proci < 0 || proci >= nProcs_)
        {
            FatalErrorInFunction
                << "Invalid processor label " << proci
                << " for face " << facei << nl
                << abort(FatalError);
        }
        else
        {
            ++nLocalFaces[proci];
        }
    }

    // Adjust lengths
    forAll(nLocalFaces, proci)
    {
        procFaceAddressing_[proci].resize(nLocalFaces[proci]);
        nLocalFaces[proci] = 0;  // restart list
    }

    // Pass 2: fill in local lists
    forAll(faceToProc_, facei)
    {
        const label proci = faceToProc_[facei];
        const label localFacei = nLocalFaces[proci];
        ++nLocalFaces[proci];

        procFaceAddressing_[proci][localFacei] = facei;
    }


    // Find processor mesh faceLabels and ...

    for (label procI = 0; procI < nProcs(); procI++)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );

        IOobject ioFvAddr
        (
            "procAddressing",
            "constant",
            polyMesh::meshSubDir,
            procFvMesh,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        );


        // pointProcAddressing (polyMesh)
        ioFvAddr.resetHeader("pointProcAddressing");
        const labelList fvPointProcAddressing
        (
            labelIOList::readContents(ioFvAddr)
        );

        Map<label> fvFaceProcAddrHash;

        {
            // faceProcAddressing (polyMesh)
            ioFvAddr.resetHeader("faceProcAddressing");
            const labelList fvFaceProcAddressing
            (
                labelIOList::readContents(ioFvAddr)
            );
            fvFaceProcAddrHash = invertToMap(fvFaceProcAddressing);
        }


        const labelList& curProcFaceAddressing = procFaceAddressing_[procI];
        labelList& curFaceLabels = procFaceLabels_[procI];

        curFaceLabels.resize_fill(curProcFaceAddressing.size(), -1);

        forAll(curProcFaceAddressing, facei)
        {
            const label polyFacei
            (
                faMesh::faceLabels()[curProcFaceAddressing[facei]]
            );
            curFaceLabels[facei] = fvFaceProcAddrHash.find(polyFacei+1).val();
        }


        // Serial mesh. Use the patch() primitive so that we also get
        // the meshPointMap
        const uindirectPrimitivePatch& basePatch = faMesh::patch();
        const Map<label>& meshPointMap = basePatch.meshPointMap();

        // Lookup of edge to edgeId
        EdgeMap<label> edgesHash;
        edgesHash.reserve(basePatch.nEdges());

        const label nIntEdges = basePatch.nInternalEdges();

        for (label edgei = 0; edgei < nIntEdges; ++edgei)
        {
            edgesHash.insert(basePatch.edges()[edgei], edgesHash.size());
        }

        for (const auto& fap : faMesh::boundary())
        {
            // Also include emptyFaPatch etc
            for (label edgei : fap.edgeLabels())
            {
                edgesHash.insert(basePatch.edges()[edgei], edgesHash.size());
            }
        }

        // Processor mesh.
        // Only need the equivalent to faMesh::patch() at this stage
        const uindirectPrimitivePatch procPatch
        (
            UIndirectList<face>(procFvMesh.faces(), curFaceLabels),
            procFvMesh.points()
        );
        const auto& procMeshPoints = procPatch.meshPoints();
        const auto& procEdges = procPatch.edges();

        labelList& curPatchPointAddressing = procPointAddressing_[procI];
        curPatchPointAddressing.resize(procMeshPoints.size(), -1);

        forAll(procMeshPoints, pointi)
        {
            curPatchPointAddressing[pointi] =
                meshPointMap[fvPointProcAddressing[procMeshPoints[pointi]]];
        }

        procNInternalEdges_[procI] = procPatch.nInternalEdges();

        auto& curPatchEdgeLookup = procPatchEdgeLookup_[procI];
        curPatchEdgeLookup.resize(procEdges.size(), -1);

        Map<label>& edgeRemapping = procMeshEdgesMap_[procI];
        edgeRemapping.clear();
        edgeRemapping.reserve(procEdges.size());

        forAll(procEdges, edgei)
        {
            edge curGlobalEdge(curPatchPointAddressing, procEdges[edgei]);

            if (auto iter = edgesHash.cfind(curGlobalEdge); iter.good())
            {
                // The edgeID (not edgeLabel) in serial
                auto serialEdgei = iter.val();

                // (key = serial edgeID; val = proc edgeLabel)
                edgeRemapping.insert(serialEdgei, edgei);

                // For each proc edgeLabel, the serial edgeID.
                // No direction encoding ever!
                curPatchEdgeLookup[edgei] = serialEdgei;
            }
            else
            {
                // Only fails if the polyMesh pointProcAddressing is corrupt
                FatalErrorInFunction
                    << "Failed edge lookup of " << curGlobalEdge << endl
                    << exit(FatalError);
            }
        }
    }
    // After this we now have
    // - procPatchEdgeLookup_ : as a lookup of the local edgeLabel
    //   (primitive patch order) to serial edgeId [no direction encoding!]
    // - procMeshEdgesMap_ : lookup serial edgeId to local edgeLabel


    Info<< "\nDistributing edges to processors" << endl;

    // Loop through all internal edges and decide which processor they
    // belong to. First visit all internal edges.

    // set references to the original mesh
    const faBoundaryMesh& patches = boundary();
    const edgeList& edges = this->edges();
    const labelUList& owner = edgeOwner();
    const labelUList& neighbour = edgeNeighbour();

    // Memory management
    {
        List<DynamicList<label>> procEdgeList(nProcs());

        // Start by transcribing internal edges ids,
        // adding direction information
        forAll(procEdgeList, proci)
        {
            const auto& curEdgeLookup = procPatchEdgeLookup_[proci];
            auto& curProcEdges = procEdgeList[proci];

            // Can directly use the primitive patch edge ids/labels
            // for the internal edges - they are identical

            curProcEdges.reserve(curEdgeLookup.size());
            curProcEdges.push_back
            (
                curEdgeLookup.slice(0, procNInternalEdges_[proci])
            );

            // Imbue with direction encoding
            for (auto& val : curProcEdges)
            {
                val += 1;  // +1 since internal edges are not flipped
            }
        }


        // Detect inter-processor boundaries
        // Track processor boundaries as (neighbour rank, edgeId)
        // lists for each subdomain
        List
        <
            DynamicList<std::pair<label, DynamicList<label>>>
        > interProcBoundaries(nProcs());

        forAll(neighbour, edgeI)
        {
            const label ownProc = faceToProc_[owner[edgeI]];
            const label neiProc = faceToProc_[neighbour[edgeI]];

            if (ownProc != neiProc)
            {
                // inter - processor patch edge found. Go through the list of
                // inside boundaries for the owner processor and try to find
                // this inter-processor patch.

                // The edge index on edge owner vs edge neighbour side
                const label ownEdgeIndex = (edgeI+1);
                const label neiEdgeIndex = -(edgeI+1);

                bool interProcBouFound = false;

                for
                (
                    auto& [ownNbrProc, ownProcEdges]
                  : interProcBoundaries[ownProc]
                )
                {
                    if (ownNbrProc == neiProc)
                    {
                        // Connection from edge owner -> edge neighbour
                        interProcBouFound = true;
                        ownProcEdges.push_back(ownEdgeIndex);

                        bool neighbourFound = false;

                        for
                        (
                            auto& [neiNbrProc, neiProcEdges]
                          : interProcBoundaries[neiProc]
                        )
                        {
                            if (neiNbrProc == ownProc)
                            {
                                // Connection from edge neighbour -> edge owner
                                neighbourFound = true;
                                neiProcEdges.push_back(neiEdgeIndex);
                            }

                            if (neighbourFound) break;
                        }

                        if (interProcBouFound && !neighbourFound)
                        {
                            FatalErrorInFunction
                                << "Inconsistency in inter-processor"
                                << " boundary lists for processors "
                                << ownProc << " and " << neiProc
                                << abort(FatalError);
                        }
                    }

                    if (interProcBouFound) break;
                }

                if (!interProcBouFound)
                {
                    // inter - processor boundaries do not exist
                    // and need to be created

                    // owner -> neighbour
                    auto& [ownNbrProc, ownProcEdges] =
                        interProcBoundaries[ownProc].emplace_back();

                    ownNbrProc = neiProc;
                    ownProcEdges.push_back(ownEdgeIndex);

                    // neighbour -> owner
                    auto& [neiNbrProc, neiProcEdges] =
                        interProcBoundaries[neiProc].emplace_back();

                    neiNbrProc = ownProc;
                    neiProcEdges.push_back(neiEdgeIndex);
                }
            }
        }


        // Loop through patches. For cyclic boundaries detect inter-processor
        // edges; for all other, add edges to the edge list and remember start
        // and size of all patches.

        // Dimension storage for tracking non-processor patches
        for (auto& range : procPatchRange_)
        {
            range.resize(patches.size());
        }

        forAll(patches, patchI)
        {
            const faPatch& fap = patches[patchI];
            const label patchStart = fap.start();

            // Init the start/size ranges for this patch (all procs)
            forAll(procPatchRange_, proci)
            {
                procPatchRange_[proci][patchI]
                    .reset(procEdgeList[proci].size(), 0);

            }

            if (!isA<cyclicFaPatch>(fap))
            {
                // Normal patch. Add edges to processor where the face
                // next to the edge lives

                const labelUList& patchEdgeLabels = fap.edgeLabels();

                forAll(patchEdgeLabels, patchEdgei)
                {
                    const label edgeLabel = patchEdgeLabels[patchEdgei];
                    const label facei = patch().edgeOwner(edgeLabel);
                    const label ownProc = faceToProc_[facei];

                    // Add to the list of edges ids
                    procEdgeList[ownProc].push_back
                    (
                        // +1 : edge remains owner-side on real boundaries
                        (patchStart + patchEdgei + 1)
                    );

                    // Increment the number of edges for this patch
                    procPatchRange_[ownProc][patchI]++;
                }
            }
            else
            {
                // Cyclic patch special treatment

                const faPatch& cPatch = patches[patchI];

                const label cycOffset = cPatch.size()/2;

                // Set reference to faceCells for both patches
                const auto firstEdgeFaces
                (
                    cPatch.edgeFaces().slice(0, cycOffset)
                );

                const auto secondEdgeFaces
                (
                    cPatch.edgeFaces().slice(cycOffset)
                );

                const auto startPart1 = (patchStart);
                const auto startPart2 = (patchStart + cycOffset);

                forAll(firstEdgeFaces, edgeI)
                {
                    const label ownProc = faceToProc_[firstEdgeFaces[edgeI]];
                    const label neiProc = faceToProc_[secondEdgeFaces[edgeI]];

                    const label firstEdgei  = (startPart1 + edgeI);
                    const label secondEdgei = (startPart2 + edgeI);

                    if (ownProc != neiProc)
                    {
                        // This edge becomes an inter-processor boundary edge
                        // inter - processor patch edge found. Go through
                        // the list of inside boundaries for the owner
                        // processor and try to find this inter-processor
                        // patch.

                        cyclicParallel_ = true;

                        bool interProcBouFound = false;

                        // The edge index on edge owner vs edge neighbour side
                        const label ownEdgeIndex = (firstEdgei+1);
                        const label neiEdgeIndex = -(secondEdgei+1);

                        for
                        (
                            auto& [ownNbrProc, ownProcEdges]
                          : interProcBoundaries[ownProc]
                        )
                        {
                            if (ownNbrProc == neiProc)
                            {
                                // Connection from edge owner -> edge neighbour
                                interProcBouFound = true;
                                ownProcEdges.push_back(ownEdgeIndex);

                                bool neighbourFound = false;

                                for
                                (
                                    auto& [neiNbrProc, neiProcEdges]
                                  : interProcBoundaries[neiProc]
                                )
                                {
                                    if (neiNbrProc == ownProc)
                                    {
                                        // From edge neighbour -> edge owner
                                        neighbourFound = true;
                                        neiProcEdges.push_back(neiEdgeIndex);
                                    }

                                    if (neighbourFound) break;
                                }

                                if (interProcBouFound && !neighbourFound)
                                {
                                    FatalErrorInFunction
                                        << "Inconsistency in inter-processor"
                                        << " boundary lists for processors "
                                        << ownProc << " and " << neiProc
                                        << " in cyclic boundary matching"
                                        << abort(FatalError);
                                }
                            }

                            if (interProcBouFound) break;
                        }

                        if (!interProcBouFound)
                        {
                            // inter - processor boundaries do not exist
                            // and need to be created

                            // owner -> neighbour
                            auto& [ownNbrProc, ownProcEdges] =
                                interProcBoundaries[ownProc].emplace_back();

                            ownNbrProc = neiProc;
                            ownProcEdges.push_back(ownEdgeIndex);

                            // neighbour -> owner
                            auto& [neiNbrProc, neiProcEdges] =
                                interProcBoundaries[neiProc].emplace_back();

                            neiNbrProc = ownProc;
                            neiProcEdges.push_back(neiEdgeIndex);
                        }
                    }
                    else
                    {
                        // This cyclic edge remains on the processor

                        // The first edge: +1 direction (stays on processor)
                        procEdgeList[ownProc].push_back(firstEdgei+1);

                        // increment the number of edges for this patch
                        procPatchRange_[ownProc][patchI]++;

                        // Note: I cannot add the other side of the cyclic
                        // boundary here because this would violate the order.
                        // They will be added in a separate loop below
                    }
                }

                // Ordering in cyclic boundaries is important.
                // Add the other half of cyclic edges for cyclic boundaries
                // that remain on the processor
                forAll(secondEdgeFaces, edgeI)
                {
                    const label ownProc = faceToProc_[firstEdgeFaces[edgeI]];
                    const label neiProc = faceToProc_[secondEdgeFaces[edgeI]];

                    //const label firstEdgei  = (startPart1 + edgeI);
                    const label secondEdgei = (startPart2 + edgeI);

                    if (ownProc == neiProc)
                    {
                        // This cyclic edge remains on the processor

                        // The second edge: +1 direction (stays on processor)
                        procEdgeList[ownProc].push_back(secondEdgei+1);

                        // increment the number of edges for this patch
                        procPatchRange_[ownProc][patchI]++;
                    }
                }
            }
        }


        // Sort processor connections according to neighbProcNo.
        // Not needed for functionality, but gives consistently
        // ordered boundaries.

        for (auto& procBoundaries : interProcBoundaries)
        {
            Foam::sort(procBoundaries);
        }

        // Add inter-processor boundaries and remember start indices
        forAll(procEdgeList, procI)
        {
            // Get internal and regular boundary processor edges
            const auto& curProcEdges = procEdgeList[procI];

            // Get reference to processor edge addressing
            labelList& curProcEdgeAddressing = procEdgeAddressing_[procI];

            labelList& curProcNeighbourProcessors =
                procNeighbourProcessors_[procI];

            auto& curProcProcessorPatchRange =
                procProcessorPatchRange_[procI];

            const auto& curInterProcBoundaries = interProcBoundaries[procI];

            const auto nProcPatches = interProcBoundaries[procI].size();

            // Flattened values for processor patches
            curProcNeighbourProcessors.resize_nocopy(nProcPatches);
            curProcProcessorPatchRange.resize_nocopy(nProcPatches);

            // Processor boundaries start after internal and regular boundaries
            label numEdges = curProcEdges.size();

            for (label procPatchi = 0; procPatchi < nProcPatches; ++procPatchi)
            {
                const auto& [bndNbrProc, bndProcEdges] =
                    curInterProcBoundaries[procPatchi];

                const auto count = bndProcEdges.size();

                curProcNeighbourProcessors[procPatchi] = bndNbrProc;
                curProcProcessorPatchRange[procPatchi].reset(numEdges, count);

                numEdges += count;
            }

            // Resize addressing
            curProcEdgeAddressing.resize(numEdges);

            // Fill in the list.
            // The turning index has already been included on all inputs,
            // so just need straight copies.

            // Internal and regular boundary edges
            numEdges = curProcEdges.size();
            curProcEdgeAddressing.slice(0, numEdges) = curProcEdges;

            // Processor boundaries
            for
            (
                const auto& [bndNbrProc, bndProcEdges]
              : curInterProcBoundaries
            )
            {
                const auto count = bndProcEdges.size();

                curProcEdgeAddressing.slice(numEdges, count) = bndProcEdges;
                numEdges += count;

                // Debug and suppress [[maybe_unused]] warning
                if (debug & 2)
                {
                    Info<< "proc" << procI << "to" << bndNbrProc
                        << " edgeLabels "; bndProcEdges.writeList(Info) << endl;
                }
            }
        }
    }

    Info<< "\nCalculating boundary addressing" << endl;
    // For every patch: the original patch index.
    // - identity for non-processor patches (ie, globally identical)
    // - '-1' for processor patches
    forAll(procPatchRange_, proci)
    {
        const auto nNonProcessorPatches = procPatchRange_[proci].size();
        const auto nProcPatches = procProcessorPatchRange_[proci].size();

        auto& curBoundaryAddressing = procBoundaryAddressing_[proci];

        curBoundaryAddressing =
            Foam::identity(nNonProcessorPatches + nProcPatches);

        // Mark processor patches
        curBoundaryAddressing.slice(nNonProcessorPatches) = -1;
    }


    // Gather data about globally shared points

    labelList globallySharedPoints_;

    // Memory management
    {
        labelList pointsUsage(nPoints(), Foam::zero{});

        // Globally shared points are the ones used by more than 2 processors
        // Size the list approximately and gather the points
        labelHashSet gSharedPoints;
        gSharedPoints.reserve(Foam::min(128, nPoints()/1000));

        // Loop through all the processors and mark up points used by
        // processor boundaries.  When a point is used twice, it is a
        // globally shared point

        for (label proci = 0; proci < nProcs(); ++proci)
        {
            // The edgeProcAddressing (edge ids)
            const auto& curEdgeAddr = procEdgeAddressing_[proci];

            // The start/size of processor boundaries
            const auto& curProcessorPatchRange =
                procProcessorPatchRange_[proci];

            // Reset the lookup list
            pointsUsage = 0;

            forAll(curProcessorPatchRange, patchI)
            {
                const auto patchEdgeIds =
                    curEdgeAddr.slice(curProcessorPatchRange[patchI]);

                for (label edgei : patchEdgeIds)
                {
                    // Mark the original points as used
                    // Remember to adjust for turning index
                    const label serialEdgei = (Foam::mag(edgei)-1);

                    const edge& e = edges[serialEdgei];

                    forAll(e, pointI)
                    {
                        if (pointsUsage[e[pointI]] == 0)
                        {
                            // Point not previously used
                            pointsUsage[e[pointI]] = patchI + 1;
                        }
                        else if (pointsUsage[e[pointI]] != patchI + 1)
                        {
                            // Point used by some other patch = global point!
                            gSharedPoints.insert(e[pointI]);
                        }
                    }
                }
            }
        }

        // Grab the result from the hash list
        globallySharedPoints_ = gSharedPoints.sortedToc();
    }


    // Edge label for faPatches

    for (label procI = 0; procI < nProcs(); procI++)
    {
        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );


        // Read volume mesh
        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        // Create processor finite-area mesh
        faMesh procMesh
        (
            areaName_,
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );

        // The serial edgeId values for the local mesh edges
        const labelList& curEdgeAddressing = procEdgeAddressing_[procI];

        const auto& curPatchRange = procPatchRange_[procI];
        const auto& curProcessorPatchRange = procProcessorPatchRange_[procI];

        const label nNonProcessorPatches = curPatchRange.size();
        const label nProcPatches = curProcessorPatchRange.size();

        labelListList& curPatchEdgeLabels = procPatchEdgeLabels_[procI];
        curPatchEdgeLabels.resize_nocopy(nNonProcessorPatches + nProcPatches);

        const Map<label>& edgeRemapping = procMeshEdgesMap_[procI];

        for (label patchi = 0; patchi < nNonProcessorPatches; ++patchi)
        {
            // [output] : edgeLabels for the patch
            auto& curEdgeLabels = curPatchEdgeLabels[patchi];

            // Patch slice of edgeProcAddressing
            const auto edgeProcAddrSlice =
                curEdgeAddressing.slice(curPatchRange[patchi]);

            curEdgeLabels.resize(edgeProcAddrSlice.size(), -1);

            // Transform from encoded serial edgeId to proc-local edgeLabel
            auto iter = curEdgeLabels.begin();
            for (label edgei : edgeProcAddrSlice)
            {
                if (auto fnd = edgeRemapping.cfind(Foam::mag(edgei)-1); fnd.good())
                {
                    *iter = fnd.val();
                }
                ++iter;
            }
        }

        for (label procPatchi = 0; procPatchi < nProcPatches; ++procPatchi)
        {
            // [output] : edgeLabels for the patch
            auto& curEdgeLabels =
                curPatchEdgeLabels[nNonProcessorPatches + procPatchi];

            // Patch slice of edgeProcAddressing
            const auto edgeProcAddrSlice =
                curEdgeAddressing.slice(curProcessorPatchRange[procPatchi]);

            curEdgeLabels.resize(edgeProcAddrSlice.size(), -1);

            // Transform from encoded serial edgeId to proc-local edgeLabel
            auto iter = curEdgeLabels.begin();
            for (label edgei : edgeProcAddrSlice)
            {
                if (auto fnd = edgeRemapping.cfind(Foam::mag(edgei)-1); fnd.good())
                {
                    *iter = fnd.val();
                }
                ++iter;
            }
        }
    }
}


void Foam::faMeshDecomposition::writeProcAddressing
(
    const faMesh& procMesh,
    const labelUList& faceProcAddr,
    const labelUList& edgeProcAddr,
    const labelUList& pointProcAddr,
    const labelUList& boundaryProcAddr,
    bool withoutEdgeEncoding
)
{
    // Processor-local outputs for components
    IOobject ioAddr
    (
        "proc-addressing",
        "constant",
        faMesh::meshSubDir,
        procMesh.thisDb(),
        IOobjectOption::NO_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    // faceProcAddressing
    ioAddr.resetHeader("faceProcAddressing");
    labelIOList::writeContents(ioAddr, faceProcAddr);

    // edgeProcAddressing
    ioAddr.resetHeader("edgeProcAddressing");
    if (withoutEdgeEncoding)
    {
        ioAddr.note() = "no edge encoding";

        // Compat: remove any direction encoding (2512 and earlier)
        const auto& input = edgeProcAddr;
        labelList plainAddressing(input.size());

        for (label i = 0; i < input.size(); ++i)
        {
            plainAddressing[i] = (Foam::mag(input[i])-1);
        }
        labelIOList::writeContents(ioAddr, plainAddressing);
    }
    else
    {
        labelIOList::writeContents(ioAddr, edgeProcAddr);
    }

    // pointProcAddressing
    ioAddr.resetHeader("pointProcAddressing");
    labelIOList::writeContents(ioAddr, pointProcAddr);

    // boundaryProcAddressing
    ioAddr.resetHeader("boundaryProcAddressing");
    labelIOList::writeContents(ioAddr, boundaryProcAddr);
}


bool Foam::faMeshDecomposition::writeDecomposition() const
{
    const word& polyMeshRegionName = faMesh::mesh().name();

    Info<< "\nConstructing processor finite-area meshes" << endl;

    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(invertToMap(globallySharedPoints_));


    label maxProcFaces = 0, totProcFaces = 0;
    label maxProcEdges = 0, totProcEdges = 0;
    label maxProcPatches = 0, totProcPatches = 0;

    // Write out the meshes
    for (label procI = 0; procI < nProcs(); procI++)
    {
        // Create processor mesh without a boundary

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        // Read volume mesh
        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        IOobject ioFvAddr
        (
            "procAddressing",
            "constant",
            polyMesh::meshSubDir,
            procFvMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // boundaryProcAddressing (polyMesh)
        ioFvAddr.resetHeader("boundaryProcAddressing");
        const labelList fvBoundaryProcAddressing
        (
            labelIOList::readContents(ioFvAddr)
        );


        // Create processor finite-area mesh
        faMesh procMesh
        (
            areaName_,
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );

        // Create processor boundary patches
        const labelList& curBoundaryAddressing =
            procBoundaryAddressing_[procI];

        const auto& curPatchRange = procPatchRange_[procI];
        const auto& curProcessorPatchRange = procProcessorPatchRange_[procI];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[procI];

        const label nNonProcessorPatches = curPatchRange.size();
        const label nProcPatches = curProcessorPatchRange.size();

        const labelListList& curPatchEdgeLabels =
            procPatchEdgeLabels_[procI];

        const faPatchList& meshPatches = boundary();

        faPatchList procPatches(nNonProcessorPatches + nProcPatches);

        label nPatches = 0;

        for (label patchi = 0; patchi < nNonProcessorPatches; ++patchi)
        {
            const auto& curEdgeLabels = curPatchEdgeLabels[nPatches];

            const label neiPolyPatchId =
                fvBoundaryProcAddressing.find
                (
                    meshPatches[curBoundaryAddressing[patchi]]
                    .ngbPolyPatchIndex()
                );

            procPatches.set
            (
                nPatches,
                meshPatches[curBoundaryAddressing[patchi]].clone
                (
                    procMesh.boundary(),
                    curEdgeLabels,
                    nPatches,
                    neiPolyPatchId
                )
            );
            ++nPatches;
        }

        for (label procPatchi = 0; procPatchi < nProcPatches; ++procPatchi)
        {
            const auto& curEdgeLabels = curPatchEdgeLabels[nPatches];

            procPatches.set
            (
                nPatches,
                new processorFaPatch
                (
                    curEdgeLabels,
                    nPatches,
                    procMesh.boundary(),
                    -1,     // nbrPolyPatch
                    procI,  // myProcNo
                    curNeighbourProcessors[procPatchi]
                )
            );

            ++nPatches;
        }

        // Add boundary patches
        procMesh.addFaPatches(procPatches);

        // More precision (for points data)
        IOstream::minPrecision(10);

        procMesh.write();

        // Statistics
        Info<< nl << "Processor " << procI;

        if (procMesh.nFaces())
        {
            Info<< nl << "    ";
        }
        else
        {
            Info<< ": ";
        }

        Info<< "Number of faces = " << procMesh.nFaces() << nl;

        if (procMesh.nFaces())
        {
            Info<< "    Number of points = " << procMesh.nPoints() << nl;
        }

        totProcFaces += procMesh.nFaces();
        maxProcFaces = Foam::max(maxProcFaces, procMesh.nFaces());

        label nBoundaryEdges = 0;
        label nProcEdges = 0;

        for (const faPatch& fap : procMesh.boundary())
        {
            if (const auto* ppp = isA<processorFaPatch>(fap); ppp)
            {
                const auto& procPatch = *ppp;

                Info<< "    Number of edges shared with processor "
                    << procPatch.neighbProcNo() << " = "
                    << procPatch.size() << nl;

                nProcEdges += procPatch.size();
            }
            else
            {
                nBoundaryEdges += fap.size();
            }
        }

        if (procMesh.nFaces() && (nBoundaryEdges || nProcEdges))
        {
            Info<< "    Number of processor patches = " << nProcPatches << nl
                << "    Number of processor edges = " << nProcEdges << nl
                << "    Number of boundary edges = " << nBoundaryEdges << nl;
        }

        totProcEdges += nProcEdges;
        totProcPatches += nProcPatches;
        maxProcEdges = Foam::max(maxProcEdges, nProcEdges);
        maxProcPatches = Foam::max(maxProcPatches, nProcPatches);

        // Write the addressing information
        writeProcAddressing
        (
            procMesh,
            procFaceAddressing_[procI],
            procEdgeAddressing_[procI],
            procPointAddressing_[procI],
            procBoundaryAddressing_[procI],
            noEdgeEncoding_
        );
    }


    // Summary stats
    Info<< nl
        << "Number of processor edges = " << (totProcEdges/2) << nl
        << "Max number of faces = " << maxProcFaces;

    if (maxProcFaces != totProcFaces)
    {
        scalar avgValue = scalar(totProcFaces)/nProcs_;

        Info<< " (" << 100.0*(maxProcFaces-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of processor patches = " << maxProcPatches;
    if (totProcPatches)
    {
        scalar avgValue = scalar(totProcPatches)/nProcs_;

        Info<< " (" << 100.0*(maxProcPatches-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of edges between processors = " << maxProcEdges;
    if (totProcEdges)
    {
        scalar avgValue = scalar(totProcEdges)/nProcs_;

        Info<< " (" << 100.0*(maxProcEdges-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl << endl;

    return true;
}


// ************************************************************************* //
