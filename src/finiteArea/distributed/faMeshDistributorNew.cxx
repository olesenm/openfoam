/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2026 OpenCFD Ltd.
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

Description
    The construction of the edge mapping uses the following process:

    Starting from the polyMesh face mapping, can easily establish the
    corresponding faMesh face mapping. The main trouble is the edge
    mapping, since any (internal) edge can be split between processors, and
    processor edges can combine to become internal edges. We have no
    information about the future connectivity on the target mesh, but
    can still use the source mesh to figure it out.

    Use half-edges for all intermediate bookkeeping, which keeps the
    edge information attached to each face so they can be shifted about
    and are not bothered by edge splitting or joining.

    - On the target mesh, get the half-edges per face. Use the face
      mapping to move (per-face) this target half-edge information onto the
      source mesh.
    - Use the existing edge connectivity to create (per face), the edge-wise
      connectivity to the neighbour target faces. This is largely correct
      except at the current processor boundaries. Do an exchange there to
      complete the connectivity information.
    - Now have a complete picture of how the future connectivity will
      look. Since we have global numbering for the target mesh faces,
      can also immediately recognize if the faces will be connected
      as an internal edge or processor edge.

    Edge flipping:
    - current internal edge either remains an internal edge
      (no flipping unless the face order changed), or becomes a processor
      edge. If the half-edge is currently the owner side: no change needed.
      If the half-edge is currently the neighbour side: needs a flip,
      since the processor has an outward pointing normal.

    - current (real) boundary edge: no flipping since it always remains
      an outward normal.

    - current processor edge. This may continue to be a processor
      edge (with the same or different combinations of ranks), in which case
      no flips are needed. Or the processor edge coalesces to become an
      internal edge. Here it needs to flip if the half-edge becomes the
      neighbour side of the internal edge.

    Now walk across the edges of each face on the source mesh
    to combine the half-edge information. The result (per source edge) will
    contain two (face, half-edge) tuples for internal edges, or a single
    (face, half-edge) tuple for boundary edges. For the boundary edges,
    the second tuple is used to encode the patch id.

\*---------------------------------------------------------------------------*/

#include "faMeshDistributor.H"
#include "globalIndex.H"
#include "BitOps.H"
#include "ListOps.H"
#include "mapDistributePolyMesh.H"
#include "processorFaPatch.H"
#include "labelPairHashes.H"

// * * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * * //

namespace
{

// Local bookkeeping tuple (face0 face1 edge0 edge1) to represent the
// joined half-edges, with face0/edge0 being the owner side and face1/edge1
// being the neighbour. For boundary edges, the face1 will be negative
// encoded with the patch id (edge1 has no meaning for this case).
//
// Inheritance from FixedList for IO, equality and contiguous properties
struct edgeTableEntry : Foam::FixedList<Foam::label, 4>
{
    using label = Foam::label;
    using labelPair = Foam::labelPair;

    // Inherit constructors
    using Foam::FixedList<Foam::label, 4>::FixedList;

    // Default construct as 'invalid'
    edgeTableEntry() : Foam::FixedList<Foam::label, 4>(-1) {}

    // Consider empty if both face indices are still invalid
    bool empty() const { return (face0() == -1 && face1() == -1); }

    // Debug: consider 'incomplete' if a -1 face index is present.
    // The current logic should create internal (or processor)
    // connections (both faces +ve) or a boundary edge '-(patchi+2)'
    bool incomplete() const { return (face0() == -1 || face1() == -1); }

    label face0() const noexcept { return this->get<0>(); }
    label face1() const noexcept { return this->get<1>(); }
    label& face0() noexcept { return this->get<0>(); }
    label& face1() noexcept { return this->get<1>(); }

    label edge0() const noexcept { return this->get<2>(); }
    label edge1() const noexcept { return this->get<3>(); }
    label& edge0() noexcept { return this->get<2>(); }
    label& edge1() noexcept { return this->get<3>(); }

    // Get side 0/1 face
    label getFace(int i) const noexcept
    {
        return (i ? this->get<1>() : this->get<0>());
    }

    // Get side 0/1 edge
    label getEdge(int i) const noexcept
    {
        return (i ? this->get<3>() : this->get<2>());
    }

    //- Return (face0 face1) as a pair
    labelPair getFaces() const { return labelPair(face0(), face1()); }

    //- Flip entries to have canonical order (face0 less than face1)
    //- for internal edges. A no-op for boundary edges.
    void canonical()
    {
        if (face1() >= 0)  // An internal edge (non-boundary neighbour)
        {
            if (face0() > face1())
            {
                std::swap(face0(), face1());
                std::swap(edge0(), edge1());
            }
            else if (face0() == face1())
            {
                // A face joined to itself - probably an error...
                if (Foam::mag(edge0()) > Foam::mag(edge1()))
                {
                    std::swap(edge0(), edge1());
                }
            }
        }
    }

    //- Remove edge flip information
    void remove_flips()
    {
        if (face0() >= 0 && edge0() < 0) { edge0() *= -1; }
        if (face1() >= 0 && edge1() < 0) { edge1() *= -1; }
    }

    //- Slot a face-edge pair into a free location in the tuple
    void add(label facei, label edgei)
    {
        if (face0() < 0)         // First slot is free
        {
            face0() = facei;
            edge0() = edgei;
        }
        else if (face1() < 0)    // Second slot is free
        {
            face1() = facei;
            edge1() = edgei;
        }
    }

    //- Partition order (for face1).
    //  In increasing order for positive values, reverse order
    //  (ie, increasing magnitude) for negative values.
    //  Useful for comparing face1 and sorting the patches to the end.
    template<class Int>
    inline static int partition_order(Int a, Int b) noexcept
    {
        if (a == b)
        {
            return 0;
        }

        auto a_neg = (a < 0);
        auto b_neg = (b < 0);

        if (a_neg != b_neg)
        {
            // Positives first
            return a_neg ? +1 : -1;
        }

        // Forward order for positive, reverse order for negative
        return (a_neg ? (a > b) : (a < b)) ? -1 : +1;
    }

    //- Comparison operator (tuple must already be in canonical order).
    bool operator<(const edgeTableEntry& other) const
    {
        // Compare face0 (owners) first
        if (this->face0() < other.face0()) return true;
        if (this->face0() > other.face0()) return false;

        // Compare face1 (neighbours) next, but sort with encoded patches
        // in increasing magnitude, not value.
        auto cmp = partition_order(this->face1(), other.face1());
        if (cmp < 0) return true;
        if (cmp > 0) return false;

        // Both faces are identical (normally does not happen but possible
        // when two edges of the owner faces are boundaries).

        // Compare edges (without flip)
        auto a = Foam::mag(this->edge0());
        auto b = Foam::mag(other.edge0());
        if (a < b) return true;
        if (b > a) return false;

        return (Foam::mag(this->edge1()) < Foam::mag(other.edge1()));
    }
};

//- Combine operation for edgeTableEntry tuples
struct edgeTableCombineOp
{
    void operator()(edgeTableEntry& x, const edgeTableEntry& y) const
    {
        // // Debug
        // if (!x.empty() && !y.empty())
        // {
        //     Foam::Pout<< "Combine " << x << " and " << y << Foam::endl;
        // }
        if (y.empty() || x == y)
        {
            // Nothing to do
        }
        else if (x.empty())
        {
            x = y;
        }
        else  // Both non-empty, but non-identical (should not occur)
        {
            using namespace Foam;
            FatalErrorInFunction << "Unexpected face edge matching: "
                << x << " vs. " << y << '\n'
                << Foam::exit(Foam::FatalError);
        }
    }
};

} // End anonymous namespace


namespace Foam
{

// The edgeTableEntry is contiguous since FixedList<label,..> is.
template<> struct is_contiguous<edgeTableEntry> : std::true_type {};

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::mapDistributePolyMesh
Foam::faMeshDistributor::distribute
(
    const faMesh& oldMesh,
    const mapDistributePolyMesh& distMap,
    const polyMesh& tgtPolyMesh,
    autoPtr<faMesh>& newMeshPtr
)
{
    newMeshPtr.reset(nullptr);

    const uindirectPrimitivePatch& oldAreaPatch = oldMesh.patch();

    // Original (patch) sizes before any changes
    const label nOldPoints = oldAreaPatch.nPoints();
    const label nOldEdges = oldAreaPatch.nEdges();
    const label nOldFaces = oldAreaPatch.nFaces();


    // ------------------------
    // Step 1: The face mapping
    // ------------------------
    // Relatively straightforward.
    // - a subset of face selections on the polyMesh faceMap

    mapDistributeBase faFaceMap(distMap.faceMap());

    {
        // Retain faces needed for the faMesh - preserve original order

        // Note can use compactLocalData without problem since
        // finiteArea is only defined on real boundary faces so there
        // is no danger of sending an internal or processor face.

        labelList oldToNewSub;
        labelList oldToNewConstruct;

        faFaceMap.compactLocalData
        (
            oldMesh.faceLabels(),
            oldToNewSub,
            oldToNewConstruct,
            distMap.nOldFaces(),
            UPstream::msgType()
        );


        // The receiving side:
        // Mapped faces (>= 0) are the polyMesh face labels that define
        // the faMesh. Since the compact order is implicitly sorted in
        // ascending order, it tends to form contiguous ranges on the
        // polyPatches, which serves our purpose nicely.

        labelList newFaceLabels
        (
            ListOps::findIndices(oldToNewConstruct, labelRange::ge0())
        );

        // Set up to read-if-present
        IOobject io(tgtPolyMesh);
        io.readOpt(IOobject::READ_IF_PRESENT);

        newMeshPtr.reset
        (
            new faMesh
            (
                oldMesh.name(),  // preserve the area-region
                tgtPolyMesh,
                std::move(newFaceLabels),
                io
            )
        );
    }

    // The face map is now complete.

    // The new faMesh and the corresponding primitive patch.
    // Cannot use any boundary addressing from the new faMesh at this
    // point, so all accounting is done with the primitive patch.
    auto& newMesh = newMeshPtr();
    const uindirectPrimitivePatch& newAreaPatch = newMesh.patch();


    // ------------------------
    // Step 2: The edge mapping
    // ------------------------
    // Use globally unique addressing to identify both sides of the edges.

    // ------------------------
    // (Edge mapping)
    // ------------------------
    // 1.
    // Create a list of destination edges for each face,
    // appended by a unique face identifier.
    // Using global face numbering from the *target* mesh.

    const auto myProci = UPstream::myProcNo();

    const globalIndex dstFaceGlobalIdx(newAreaPatch.nFaces());

    // The faceEdges() from the destination mesh, with an additional entry
    // marking its unique face id.
    labelListList dstFaceEdges(newAreaPatch.nFaces());

    forAll(newAreaPatch.faceEdges(), facei)
    {
        const auto& fEdges = newAreaPatch.faceEdges()[facei];
        const auto dstFacei = dstFaceGlobalIdx.toGlobal(myProci, facei);

        // The destination (patch) edge indices, with unique destination
        // face identifier:

        auto& dstEdges = dstFaceEdges[facei];
        dstEdges.resize(fEdges.size() + 1);
        dstEdges[fEdges.size()] = dstFacei;

        // Edge (flip) encoding is done later
        dstEdges.slice(0, fEdges.size()) = fEdges;
    }

    // Send back to the original mesh locations
    faFaceMap.reverseDistribute(nOldFaces, dstFaceEdges);

    // Establish information about the destination faces and their
    // neighbours.

    // First grab the destination face id (unique, globally consistent)
    // from the end of the lists that just arrived via reverseDistribute

    labelList dstFaceIds(nOldFaces);
    for (label facei = 0; facei < nOldFaces; ++facei)
    {
        dstFaceIds[facei] = dstFaceEdges[facei].back();
    }

    // Next use the edge connectivity of the current (source) mesh
    // to define the destination neighbour faces for each faceEdge

    labelListList dstFaceNeighbours(nOldFaces);
    for (label facei = 0; facei < nOldFaces; ++facei)
    {
        // Transform the input list of edges (source mesh)
        // into a list of neighbour faces (target mesh)
        const auto& srcEdges = oldAreaPatch.faceEdges()[facei];

        auto& dstNbrs = dstFaceNeighbours[facei];
        dstNbrs.resize(srcEdges.size());

        forAll(srcEdges, i)
        {
            auto srcEdgei = srcEdges[i];
            auto& dstNbr = dstNbrs[i];

            label srcNbr = oldAreaPatch.edgeNeighbour(srcEdgei);
            if (srcNbr == facei)
            {
                srcNbr = oldAreaPatch.edgeOwner(srcEdgei);
            }

            // Neighbour face (target mesh)
            dstNbr = (srcNbr >= 0 ? dstFaceIds[srcNbr] : -1);
        }
    }

    // Note that 'dstFaceNeighbours' will encode the owner/neighbour
    // relationship for the target mesh, but be represented on the source
    // mesh. It acts like full halo information for the faces.
    // This is augmented by 'boundaryEdgeNbrs', which represents the
    // associated destination edge index, also a halo information.

    // The neighbour edge information
    labelList boundaryEdgeNbrs(oldAreaPatch.nBoundaryEdges());
    boundaryEdgeNbrs = -1;

    // Make neighbour face information consistent across ranks
    if (UPstream::parRun())
    {
        const auto startOfRequests = UPstream::nRequests();

        const label oldBndStart = oldAreaPatch.nInternalEdges();

        // Setup sends
        for (const faPatch& fap : oldMesh.boundary())
        {
            if (const auto* fapp = isA<processorFaPatch>(fap); fapp)
            {
                // Send the face/edge information
                List<labelPair> nbrInfo(fap.edgeLabels().size());

                auto* iter = nbrInfo.begin();
                for (auto srcEdgei : fap.edgeLabels())
                {
                    auto srcFacei = oldAreaPatch.edgeOwner(srcEdgei);
                    auto dstFacei = dstFaceIds[srcFacei];
                    auto dstEdgei = srcEdgei;

                    // The index of this edge on the face
                    auto i = oldAreaPatch.faceEdges()[srcFacei].find(srcEdgei);
                    if (i >= 0)
                    {
                        dstEdgei = dstFaceEdges[srcFacei][i];
                    }
                    else
                    {
                        // Failure should be impossible
                        FatalErrorInFunction
                            << "Bad indexing for edge:"
                            << srcEdgei << " on face:" << srcFacei << endl
                            << exit(FatalError);
                    }

                    *iter = labelPair(dstFacei, dstEdgei);
                    ++iter;
                }

                fapp->nonblocking_send(nbrInfo);
            }
        }

        UPstream::waitRequests(startOfRequests);

        // Receive values
        for (const faPatch& fap : oldMesh.boundary())
        {
            if (const auto* fapp = isA<processorFaPatch>(fap); fapp)
            {
                // Receive the face/edge ids (in target mesh ids) and
                // update the face neighbour and boundary edge information
                // accordingly

                List<labelPair> nbrInfo(fap.edgeLabels().size());
                fapp->nonblocking_receive(nbrInfo);

                const auto* iter = nbrInfo.cbegin();
                for (auto srcEdgei : fap.edgeLabels())
                {
                    auto srcFacei = oldAreaPatch.edgeOwner(srcEdgei);
                    auto nbrFacei = (*iter).first();
                    auto nbrEdgei = (*iter).second();
                    ++iter;

                    // The index of this edge on the face
                    auto i = oldAreaPatch.faceEdges()[srcFacei].find(srcEdgei);
                    if (i >= 0)
                    {
                        dstFaceNeighbours[srcFacei][i] = nbrFacei;
                        boundaryEdgeNbrs[srcEdgei-oldBndStart] = nbrEdgei;
                    }
                    else
                    {
                        // Failure should be impossible
                        FatalErrorInFunction
                            << "Bad indexing for edge:"
                            << srcEdgei << " on face:" << srcFacei << endl
                            << exit(FatalError);
                    }
                }
            }
        }
    }

    // Adjust neighbour face information for boundaries.
    // Real (non-processor) patches : imprint with encoded patch index.

    label nNonProcessor(0);
    {
        const faBoundaryMesh& oldBndMesh = oldMesh.boundary();
        const label oldBndStart = oldAreaPatch.nInternalEdges();

        forAll(oldBndMesh, patchi)
        {
            // Encoded as -(patch+2) to distinguish from unconnected (-1)
            const label encodedPatchId = -(patchi+2);

            const faPatch& fap = oldBndMesh[patchi];

            if (isA<processorFaPatch>(fap))
            {
                break;
            }
            ++nNonProcessor;

            for (auto srcEdgei : fap.edgeLabels())
            {
                auto srcFacei = oldAreaPatch.edgeOwner(srcEdgei);

                // The index of this edge on the face
                auto i = oldAreaPatch.faceEdges()[srcFacei].find(srcEdgei);
                if (i >= 0)
                {
                    dstFaceNeighbours[srcFacei][i] = encodedPatchId;
                    boundaryEdgeNbrs[srcEdgei-oldBndStart] = encodedPatchId;
                }
                else
                {
                    // Failure should be impossible
                    FatalErrorInFunction
                        << "Bad indexing for edge:"
                        << srcEdgei << " on face:" << srcFacei << endl
                        << exit(FatalError);
                }
            }
        }
    }
    Pstream::broadcast(nNonProcessor);

    #ifdef FULLDEBUG
    if (auto badIndex = boundaryEdgeNbrs.find(-1); badIndex >= 0)
    {
        FatalErrorInFunction
            << "boundaryEdgeNbrs has unvisited location(s): "
            << badIndex
            << exit(FatalError);
    }
    #endif

    // Join everything together

    // 2.
    // Walk all original faces and their edges to generate a edge lookup
    // table with the destination face/edge information.
    // Eg ((globFacei, localEdgei), (globFacei, localEdgei))

    // NB: currently no provision for polyMesh face flips,
    // which would reverse the face edge order.

    List<edgeTableEntry> edgeTable(nOldEdges, edgeTableEntry());
    DynamicList<bool> edgeFlips;

    for (label facei = 0; facei < nOldFaces; ++facei)
    {
        // Source mesh information
        const auto& srcEdges = oldAreaPatch.faceEdges()[facei];

        // Target mesh information
        const auto& dstEdges = dstFaceEdges[facei];
        const auto& dstNbrs = dstFaceNeighbours[facei];
        const auto dstFacei = dstFaceIds[facei];
        const auto dstProc = dstFaceGlobalIdx.whichProcID(dstFacei);

        const label nFaceEdges = srcEdges.size();

        // Decide on the edge flipping needed.
        // A separate loop just to keep the logic somewhat clearer.

        edgeFlips.clear();

        forAll(srcEdges, i)
        {
            auto srcEdgei = srcEdges[i];
            auto nbrFacei = dstNbrs[i];
            bool flip = false;

            if (srcEdgei < oldAreaPatch.nInternalEdges())
            {
                // Internal edge (entirely on one rank)

                bool srcOwner = (facei == oldAreaPatch.edgeOwner(srcEdgei));
                bool dstOwner = (dstFacei < nbrFacei);
                auto dstProcNbr = dstFaceGlobalIdx.whichProcID(nbrFacei);

                if (dstProc == dstProcNbr)
                {
                    // Remains an internal edge, but owner/neigh order may
                    // change (on the same proc, or moved entirely to a new
                    // proc).

                    flip = (srcOwner != dstOwner);
                }
                else
                {
                    // Was internal edge, now a processor connection.
                    // Needs flip if not originally the owner.

                    flip = (!srcOwner);
                }
            }
            else if (nbrFacei < 0)  // A real boundary edge
            {
                // Nothing to do.
            }
            else
            {
                // A processor boundary edge: so srcOwner==true

                bool dstOwner = (dstFacei < nbrFacei);
                auto dstProcNbr = dstFaceGlobalIdx.whichProcID(nbrFacei);

                if (dstProc == dstProcNbr)
                {
                    // Became an internal connection.
                    // Flip if it becomes the neighbour side.
                    flip = (!dstOwner);
                }
                else
                {
                    // Was a processor boundary and continued to be one.
                    // Nothing to do.
                }
            }

            edgeFlips.push_back(flip);
        }

        const auto oldBndStart = oldAreaPatch.nInternalEdges();

        // Now add into edgeTable information
        for (label edgei = 0; edgei < nFaceEdges; ++edgei)
        {
            label srcEdgei = srcEdges[edgei];
            label dstEdgei =
            (
                edgeFlips[edgei]
              ? -(dstEdges[edgei]+1)
              : (dstEdges[edgei]+1)
            );

            edgeTable[srcEdgei].add(dstFacei, dstEdgei);

            // For boundary edges, add real patch information or processor
            // halo information into the second slot

            if (srcEdgei >= oldAreaPatch.nInternalEdges())
            {
                auto nbrFacei = dstNbrs[edgei];
                auto nbrEdgei = (boundaryEdgeNbrs[srcEdgei-oldBndStart]+1);

                if (nbrFacei < -1)
                {
                    // Real boundary
                    edgeTable[srcEdgei].add(nbrFacei, -1);
                }
                else
                {
                    // Processor boundary (halo information)
                    edgeTable[srcEdgei].add(nbrFacei, nbrEdgei);
                }
            }
        }
    }


    // Globally consistent order
    for (auto& tuple : edgeTable)
    {
        tuple.canonical();
    }


    // Now ready to actually construct the edgeMap

    mapDistributeBase faEdgeMap
    (
        newAreaPatch.nEdges(),  // constructSize
        labelListList(),        // subMap
        labelListList(),        // constructMap
        true,                   // subHasFlip
        true,                   // constructHasFlip
        faFaceMap.comm()
    );

    {
        // Pass 1.
        // Count the number of edges to be sent to each proc

        labelList nProcEdges(UPstream::nProcs(faFaceMap.comm()), Foam::zero{});

        for (const auto& tuple : edgeTable)
        {
            labelPair dstProcs(-1, -1);

            for (int sidei = 0; sidei < 2; ++sidei)
            {
                const label dstFacei = tuple.getFace(sidei);
                //const label dstEdgei = tuple.getEdge(sidei);

                if (dstFacei < 0)
                {
                    // Neighbour is a patch
                    continue;
                }

                label proci = dstFaceGlobalIdx.whichProcID(dstFacei);
                dstProcs[sidei] = proci;

                // Always includes side0, but ignores side1 if it has
                // the same target rank (avoids double counting)
                if (proci >= 0 && dstProcs[0] != dstProcs[1])
                {
                    ++nProcEdges[proci];
                }
            }
        }

        auto& edgeSubMap = faEdgeMap.subMap();
        auto& edgeCnstrMap = faEdgeMap.constructMap();

        edgeSubMap.resize(nProcEdges.size());
        edgeCnstrMap.resize(nProcEdges.size());

        labelListList remoteEdges(nProcEdges.size());

        forAll(nProcEdges, proci)
        {
            edgeSubMap[proci].resize(nProcEdges[proci], -1);
            remoteEdges[proci].resize(nProcEdges[proci], -1);
        }


        // Pass 2.
        // Fill in the maps

        nProcEdges = Zero;  // Reset counter

        forAll(edgeTable, edgei)
        {
            const auto& tuple = edgeTable[edgei];
            const auto srcEdgei = (edgei+1);  // with flip encoding

            labelPair dstProcs(-1, -1);

            for (int sidei = 0; sidei < 2; ++sidei)
            {
                auto dstFacei = tuple.getFace(sidei);
                auto dstEdgei = tuple.getEdge(sidei);

                if (dstFacei < 0)
                {
                    // Neighbour is a patch
                    continue;
                }

                label proci = dstFaceGlobalIdx.whichProcID(dstFacei);
                dstProcs[sidei] = proci;

                // Always includes side0, but ignores side1 if it has
                // the same target rank (avoids double counting)
                if (proci >= 0 && dstProcs[0] != dstProcs[1])
                {
                    edgeSubMap[proci][nProcEdges[proci]] = srcEdgei;
                    remoteEdges[proci][nProcEdges[proci]] = dstEdgei;
                    ++nProcEdges[proci];
                }
            }
        }

        // The remoteEdges are what we know locally about what will be
        // received, but not what is actually received.
        // So need an all-to-all exchange

        Pstream::exchange<labelList, label>
        (
            remoteEdges,
            edgeCnstrMap,
            UPstream::msgType(),
            faEdgeMap.comm()
        );
    }

    // The edge map is now complete [in PrimitivePatch edge order]
    //
    // IMPORTANT:
    // The edge flip information in edgeTable was useful for creating the
    // edge map itself, but now must be removed. The half-edge flipping
    // is already incorporated into edgeMap and if not removed here,
    // the edgeTable will not be globally consistent.
    //
    // Note that the edge numbers are still needed to ensure that both
    // processor sides are consistently ordered (cannot just use the faces).

    for (auto& tuple : edgeTable)
    {
        tuple.remove_flips();
    }

    if (oldMesh.hasInternalEdgeLabels())
    {
        // If there are gaps in the edge numbering or the
        // internal edge labels are out of sequence would
        // have to use compactLocalData etc before sending
        // But just issue an error for now

        FatalErrorInFunction
            << "Originating faMesh has gaps in the edge addressing"
            << " this is currently unsupported"
            << abort(FatalError);
    }


    // ------------------------
    // Patch edge labels
    // ------------------------

    // Distribute the edge lookups.
    // Needs full version for the combine operation

    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        UList<labelPair>::null(),
        faEdgeMap.constructSize(),

        faEdgeMap.subMap(),
        faEdgeMap.subHasFlip(),

        faEdgeMap.constructMap(),
        faEdgeMap.constructHasFlip(),

        edgeTable,
        edgeTableEntry(),       // nullValue
        edgeTableCombineOp{},   // CombineOp
        identityOp{},           // No negation (already globally consistent)

        UPstream::msgType(),
        faEdgeMap.comm()
    );

    faPatchList newFaPatches;
    {
        // Pass 1.
        // Count the number of edges for each patch type

        labelList nEdgeLabels(nNonProcessor, Foam::zero{});
        LabelPairMap<label> nProcEdges;

        forAll(edgeTable, edgei)
        {
            const auto& tuple = edgeTable[edgei];

            labelPair target(tuple.getFaces());

            if (target[1] < -1)
            {
                // Neighbour face was patchId encoded value
                label patchi = (Foam::mag(target[1])-2);
                ++nEdgeLabels[patchi];
            }
            else if (target[0] >= 0 && target[1] >= 0)
            {
                // From global face to proc id
                target[0] = dstFaceGlobalIdx.whichProcID(target[0]);
                target[1] = dstFaceGlobalIdx.whichProcID(target[1]);

                // A processor-processor connection involving my rank...
                if ((target[0] != target[1]) && target.contains(myProci))
                {
                    // Should already be sorted, since the faces were
                    target.sort();
                    ++nProcEdges(target);
                }
            }
        }

        const label nPatches = (nNonProcessor + nProcEdges.size());

        newFaPatches.resize(nPatches);
        nEdgeLabels.resize(nPatches, Foam::zero{});

        labelListList newEdgeLabels(nPatches);
        LabelPairMap<label> procPairToPatchId;

        // Presize edgeLabels arrays,
        // map processor pairs to an index in patches
        {
            for (label patchi = 0; patchi < nNonProcessor; ++patchi)
            {
                label nLabels = nEdgeLabels[patchi];
                newEdgeLabels[patchi].resize(nLabels);
            }

            // Processor patches - sorted order
            label patchi = nNonProcessor;
            for (const labelPair& twoProcs : nProcEdges.sortedToc())
            {
                label nLabels = nProcEdges.lookup(twoProcs, 0);

                nEdgeLabels[patchi] = nLabels;
                newEdgeLabels[patchi].resize(nLabels);

                procPairToPatchId.set(twoProcs, patchi);
                ++patchi;
            }
        }
        nEdgeLabels = Zero;  // Reset counter


        // Populate edgeLabels arrays - walk in canonically sorted
        // order to ensure that both sides of processor edges
        // correspond.

        // Assumes that each edgeTable element is in canonical order
        const labelList order(Foam::sortedOrder(edgeTable));

        for (const label edgei : order)
        {
            const auto& tuple = edgeTable[edgei];

            labelPair target(tuple.getFaces());

            label patchi = -1;

            if (target[1] < -1)
            {
                // Neighbour face was patchId encoded value
                patchi = (Foam::mag(target[1])-2);
            }
            else if (target[0] >= 0 && target[1] >= 0)
            {
                // From global face to proc id
                target[0] = dstFaceGlobalIdx.whichProcID(target[0]);
                target[1] = dstFaceGlobalIdx.whichProcID(target[1]);

                // A processor-processor connection involving my rank...
                if ((target[0] != target[1]) && target.contains(myProci))
                {
                    // Should already be sorted, since the faces were
                    target.sort();
                    patchi = procPairToPatchId.lookup(target, -1);
                }
            }

            if (patchi >= 0)
            {
                auto idx = nEdgeLabels[patchi]++;
                newEdgeLabels[patchi][idx] = edgei;
            }
        }


        // Clone all non-processor patches
        for (label patchi = 0; patchi < nNonProcessor; ++patchi)
        {
            const auto& oldPatch = oldMesh.boundary()[patchi];

            newFaPatches.set
            (
                patchi,
                oldPatch.clone
                (
                    newMesh.boundary(),
                    newEdgeLabels[patchi],  // edgeLabels
                    patchi,
                    oldPatch.ngbPolyPatchIndex()
                )
            );
        }

        // Create any processor patches
        forAllConstIters(procPairToPatchId, iter)
        {
            const auto& twoProcs = iter.key();
            const auto patchi = iter.val();

            auto nbrProcNo =
            (
                twoProcs.first() != myProci
              ? twoProcs.first() : twoProcs.second()
            );

            newFaPatches.set
            (
                patchi,
                new processorFaPatch
                (
                    newEdgeLabels[patchi],  // edgeLabels
                    patchi,
                    newMesh.boundary(),
                    -1,                     // nbrPolyPatchi
                    myProci,
                    nbrProcNo
                )
            );
        }
    }


    newMesh.addFaPatches(newFaPatches);
    newMesh.init(true);


    // At this stage we now have a complete mapping overview in terms of
    // the PrimitivePatch edge ordering. Now need to adjust for the
    // different boundary edge order (internal edges are one-to-one)

    {
        const auto& boundaries = oldMesh.boundary();
        auto& updateMap = faEdgeMap.subMap();
        const bool hasFlip = faEdgeMap.subHasFlip();

        // Boundary edgeLabels -> edgeId
        label start = boundaries.start();
        Map<label> oldToNew;

        for (const faPatch& p : boundaries)
        {
            for (auto edgei : p.edgeLabels())
            {
                oldToNew.insert(edgei, start++);
            }
        }

        mapDistributeBase::renumberMap(updateMap, oldToNew, hasFlip);
    }

    {
        const auto& boundaries = newMesh.boundary();
        auto& updateMap = faEdgeMap.constructMap();
        const bool hasFlip = faEdgeMap.constructHasFlip();

        // Boundary edgeLabels -> edgeId
        label start = boundaries.start();
        Map<label> oldToNew;

        for (const faPatch& p : boundaries)
        {
            for (auto edgei : p.edgeLabels())
            {
                oldToNew.insert(edgei, start++);
            }
        }

        mapDistributeBase::renumberMap(updateMap, oldToNew, hasFlip);
    }


    // ------------------------
    // Patch mapping
    // ------------------------

    mapDistributeBase faPatchMap
    (
        newMesh.boundary().size(),  // constructSize
        labelListList(),        // subMap
        labelListList(),        // constructMap
        false,                  // subHasFlip
        false,                  // constructHasFlip
        faFaceMap.comm()
    );

    // For patch maps, would normally transcribe from patchMapInfo
    // gathered earlier. However, existing practice (APR-2022) for
    // faMesh decomposition is to map all non-processor patches

    {
        // Map all non-processor patches
        const label nProcs = UPstream::nProcs(faPatchMap.comm());

        faPatchMap.subMap().resize(nProcs, identity(nNonProcessor));
        faPatchMap.constructMap().resize(nProcs, identity(nNonProcessor));
    }


    // ------------------------
    // Point mapping
    // ------------------------

    mapDistributeBase faPointMap(distMap.pointMap());

    {
        // Retain meshPoints needed for the faMesh - preserve original order
        // Need both sides (local/remote) for correct compaction maps
        // without dangling points.

        labelList oldToNewSub;
        labelList oldToNewConstruct;

        faPointMap.compactData
        (
            oldAreaPatch.meshPoints(),
            newAreaPatch.meshPoints(),
            oldToNewSub,
            oldToNewConstruct,
            distMap.nOldPoints(),
            UPstream::msgType()
        );
    }


    return mapDistributePolyMesh
    (
        // Mesh before changes
        nOldPoints,
        nOldEdges,          // area: nOldEdges (volume: nOldFaces)
        nOldFaces,          // area: nOldFaces (volume: nOldCells)

        labelList(oldMesh.boundary().patchStarts()),
        labelList(),        // oldPatchNMeshPoints [unused]

        mapDistribute(std::move(faPointMap)),
        mapDistribute(std::move(faEdgeMap)), // area: edgeMap (volume: faceMap)
        mapDistribute(std::move(faFaceMap)), // area: faceMap (volume: cellMap)
        mapDistribute(std::move(faPatchMap))
    );
}


Foam::mapDistributePolyMesh
Foam::faMeshDistributor::distribute
(
    const faMesh& oldMesh,
    const mapDistributePolyMesh& distMap,
    autoPtr<faMesh>& newMeshPtr
)
{
    return faMeshDistributor::distribute
    (
        oldMesh,
        distMap,
        oldMesh.mesh(),  // polyMesh
        newMeshPtr
    );
}


// ************************************************************************* //
