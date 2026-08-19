/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2026 OpenCFD Ltd.
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

#include "foamVtuSizing.H"
#include "foamVtkCore.H"
#include "cellShape.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "manifoldCellsMeshObject.H"

// Only used in this file
#include "foamVtuSizingImpl.cxx"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::vtuSizing::presizeMaps(foamVtkMeshMaps& maps) const
{
    maps.cellMap().resize_nocopy(this->nFieldCells());
    maps.additionalIds().resize_nocopy(this->nAddPoints());
}


void Foam::vtk::vtuSizing::checkSizes
(
    const vtk::vtuSizing& sizing,

    const label cellTypes_size,
    const label vertLabels_size,
    const label vertOffset_size,
    const label faceLabels_size,
    const label faceOffset_size,
    const label polyFaceIds_size,
    const label polyFaceOffset_size,

    const enum contentType output,
    const label cellMap_size,
    const label addPointsIds_size
)
{
    label nErrors = 0;

    #undef  CHECK_SIZING
    #define CHECK_SIZING(what, sizeInput, sizeExpected)        \
    if (sizeInput != sizeExpected)                             \
    {                                                          \
        if (!nErrors++)                                        \
        {                                                      \
            FatalErrorInFunction << "VTK sizing error" << nl;  \
        }                                                      \
        FatalError                                             \
            << "    " << what << " size=" << sizeInput         \
            << " expected " << sizeExpected << nl;             \
    }


    CHECK_SIZING("cellTypes", cellTypes_size, sizing.nFieldCells());
    CHECK_SIZING("cellMap", cellMap_size, sizing.nFieldCells());
    CHECK_SIZING("addPointsIds", addPointsIds_size, sizing.nAddPoints());

    switch (output)
    {
        case contentType::LEGACY :
        {
            // Legacy uses cells for everything
            CHECK_SIZING
            (
                "legacy",
                vertLabels_size,
                sizing.sizeOf<slotType::CELLS>(output)
            );
            break;
        }

        case contentType::XML :
        case contentType::INTERNAL2 :
        {
            // XML, INTERNAL uses connectivity/offset pairs
            CHECK_SIZING
            (
                "connectivity",
                vertLabels_size,
                sizing.sizeOf<slotType::CELLS>(output)
            );
            CHECK_SIZING
            (
                "offsets",
                vertOffset_size,
                sizing.sizeOf<slotType::CELLS_OFFSETS>(output)
            );
            if (sizing.hasPolyCells())
            {
                CHECK_SIZING
                (
                    "faces",
                    faceLabels_size,
                    sizing.sizeOf<slotType::FACES>(output)
                );
                CHECK_SIZING
                (
                    "faceOffsets",
                    faceOffset_size,
                    sizing.sizeOf<slotType::FACES_OFFSETS>(output)
                );
            }
            break;
        }

        case contentType::HDF :
        {
            // VTKHDF connectivity/offset pairs
            CHECK_SIZING
            (
                "Connectivity",
                vertLabels_size,
                sizing.sizeOf<slotType::CELLS>(output)
            );
            CHECK_SIZING
            (
                "Offsets",
                vertOffset_size,
                sizing.sizeOf<slotType::CELLS_OFFSETS>(output)
            );
            if (sizing.hasPolyCells())
            {
                CHECK_SIZING
                (
                    "FaceConnectivity",
                    faceLabels_size,
                    sizing.sizeOf<slotType::FACES>(output)
                );
                CHECK_SIZING
                (
                    "FaceOffsets",
                    faceOffset_size,
                    sizing.sizeOf<slotType::FACES_OFFSETS>(output)
                );
                CHECK_SIZING
                (
                    "PolyhedronOffsets",
                    polyFaceOffset_size,
                    sizing.sizeOf<slotType::POLY_FACEIDS_OFFSETS>(output)
                );

                // Note: polyFaceIds may be zero-sized (if demand-driven)
                // So ignore any size mismatch there.
                if (polyFaceIds_size)
                {
                    CHECK_SIZING
                    (
                        "PolyhedronToFaces",
                        polyFaceIds_size,
                        sizing.sizeOf<slotType::POLY_FACEIDS>(output)
                    );
                }
            }
            break;
        }
    }

    if (nErrors)
    {
        FatalError
            << nl
            << "Total of " << nErrors << " sizing errors encountered!"
            << exit(FatalError);
    }

    #undef CHECK_SIZING
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList
Foam::vtk::vtuSizing::dummyFaceOffsets
(
    const label numCells,
    const enum contentType output,
    label beginOffset
)
{
    labelList offsets;

    switch (output)
    {
        case contentType::LEGACY :
        {
            // Not applicable
            break;
        }

        case contentType::XML :
        case contentType::INTERNAL2 :
        {
            // Primitive cells: -1 placeholder
            offsets.resize(numCells, -1);
            break;
        }

        case contentType::HDF :
        {
            // Primitive cells: zero-sized local sizes
            // NB: this entry point is likely unused
            if (numCells)
            {
                offsets.resize(numCells+1, beginOffset);
            }
            break;
        }
    }
    return offsets;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::vtuSizing::vtuSizing() noexcept
{
    clear();
}


Foam::vtk::vtuSizing::vtuSizing
(
    const polyMesh& mesh,
    const bool decompose
)
{
    clear();
    reset(mesh, decompose);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::vtk::vtuSizing::nFaceLabels(contentType output) const
{
    if (contentType::HDF == output)
    {
        // Just the number of face labels used
        // == sum of [id0,id1, ... ] entries
        return nFaceLabels_;
    }
    else
    {
        // The size for [nFaces, nFace0Pts, id0,id1,..., ] for all cells
        return
        (
            nCellsPoly_   // Sum all [nFaces] entries
          + nFacesPoly_   // Sum all [nFace0Pts, ..., nFace1Pts, ...] entries
          + nFaceLabels_  // Sum all [id0,id1, ... ] entries
        );
    }
}


void Foam::vtk::vtuSizing::clear() noexcept
{
    decompose_   = false;
    manifold_    = false;
    selectionMode_ = selectionModeType::FULL_MESH;
    nCells_      = 0;
    nPoints_     = 0;
    nVertLabels_ = 0;

    nFaceLabels_ = 0;
    nCellsPoly_  = 0;
    nFacesPoly_  = 0;
    nVertPoly_   = 0;

    nAddCells_   = 0;
    nAddPoints_  = 0;
    nAddVerts_   = 0;
}


void Foam::vtk::vtuSizing::reset
(
    const polyMesh& mesh,
    const bool decompose
)
{
    reset(mesh, labelUList::null(), decompose);
}


void Foam::vtk::vtuSizing::reset
(
    const polyMesh& mesh,
    const labelUList& subsetCellsIds,
    const bool decompose
)
{
    // References to cell shape models
    const cellModel& tet   = cellModel::ref(cellModel::TET);
    const cellModel& pyr   = cellModel::ref(cellModel::PYR);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& hex   = cellModel::ref(cellModel::HEX);
    const cellModel& wedge = cellModel::ref(cellModel::WEDGE);
    const cellModel& tetWedge = cellModel::ref(cellModel::TETWEDGE);

    const cellShapeList& shapes = mesh.cellShapes();
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();

    // Special treatment for mesh subsets.
    const bool isSubsetMesh
    (
        notNull(subsetCellsIds)
    );

    if (isSubsetMesh)
    {
        decompose_  = false;  // Disallow decomposition for subset mode
        selectionMode_ = selectionModeType::SUBSET_MESH;
    }
    else
    {
        decompose_  = decompose;  // Permit decomposition (if requested)
        selectionMode_ = selectionModeType::FULL_MESH;
    }

    // Manifold cells detected?
    manifold_ = manifoldCellsMeshObject::New(mesh).manifold();

    const label nInputCells =
    (
        isSubsetMesh
      ? subsetCellsIds.size()
      : shapes.size()
    );

    nCells_    = nInputCells;
    nPoints_   = mesh.nPoints();
    nAddCells_ = 0;
    nAddVerts_ = 0;

    nCellsPoly_  = 0;
    nFacesPoly_  = 0;
    nVertLabels_ = 0;
    nFaceLabels_ = 0;
    nVertPoly_   = 0;

    // The requested polyhedral decomposition was actually needed?
    bool neededDecompose(false);

    // Unique vertex labels per polyhedral
    labelHashSet hashUniqId;
    if (!decompose_) { hashUniqId.reserve(256); }

    // For polyhedral decomposition, ensure that the face tri/quad splits
    // on the neighbour side of a processor patch follow that of the owner side.
    bitSet flipBoundaryFace;
    if (decompose_ && UPstream::parRun())
    {
        flipBoundaryFace.resize(mesh.nBoundaryFaces());

        for (const polyPatch& pp : mesh.boundaryMesh())
        {
            if
            (
                const auto* ppp = isA<processorPolyPatch>(pp);
                (ppp && ppp->neighbour())
            )
            {
                // Neighbour-side of processor patch
                flipBoundaryFace.set(labelRange(pp.offset(), pp.size()));
            }
        }
    }

    for (label inputi = 0; inputi < nInputCells; ++inputi)
    {
        const label celli(isSubsetMesh ? subsetCellsIds[inputi] : inputi);

        const cellShape& shape = shapes[celli];
        const cellModel& model = shape.model();

        if
        (
            model == tet
         || model == pyr
         || model == prism
         || model == hex
        )
        {
            // Normal primitive - not a poly
            nVertLabels_ += shape.size();
        }
        else if (model == tetWedge && decompose_)
        {
            neededDecompose = true;
            nVertLabels_ += 6;  // Treat as squeezed prism (VTK_WEDGE)
        }
        else if (model == wedge && decompose_)
        {
            neededDecompose = true;
            nVertLabels_ += 8;  // Treat as squeezed hex
        }
        else if (decompose_)
        {
            // Polyhedral: Decompose into tets + pyramids.
            neededDecompose = true;
            ++nAddPoints_;

            // Count vertices into first decomposed cell
            bool firstCell = true;

            for (const label facei : meshCells[celli])
            {
                // Face decomposed into triangles and quads
                // Tri -> Tet, Quad -> Pyr
                const face& f = meshFaces[facei];
                label nTria = 0, nQuad = 0;

                // Need to use a flipped face on the boundary?
                if
                (
                    (f.size() > 4)
                 && flipBoundaryFace.test(facei - mesh.nInternalFaces())
                )
                {
                    const face flipped(f.reverseFace());
                    flipped.nTrianglesQuads(mesh.points(), nTria, nQuad);
                }
                else
                {
                    f.nTrianglesQuads(mesh.points(), nTria, nQuad);
                }

                nAddCells_ += nTria + nQuad;
                nAddVerts_ += (nTria * 4) + (nQuad * 5);

                if (firstCell)
                {
                    firstCell = false;
                    --nAddCells_;

                    const label nvrt = (nQuad ? 5 : 4);
                    nAddVerts_   -= nvrt;
                    nVertLabels_ += nvrt;
                }
            }
        }
        else
        {
            // Polyhedral: Not decomposed

            // The face stream is often like this:
            // number of faces, size of each face, vertices per face
            // [nFaces, nFace0Pts, id0,id1,..., nFace1Pts, id0,...]
            //
            // but with VTKHDF the sizing is handled separately,
            // so keep bookkeeping separate at this stage.

            const labelList& cFaces = meshCells[celli];

            ++nCellsPoly_;
            nFacesPoly_ += cFaces.size();

            // Unique node ids used
            hashUniqId.clear();

            for (const label facei : cFaces)
            {
                const face& f = meshFaces[facei];
                nFaceLabels_ += f.size();

                hashUniqId.insert(f);
            }

            // Legacy format only uses the face-stream.
            // - track what *NOT* to use for legacy
            nVertLabels_ += hashUniqId.size();
            nVertPoly_   += hashUniqId.size();
        }
    }

    // Requested and actually required
    decompose_ = (decompose_ && neededDecompose);
}


// Synchronize changes here with the following:
// - vtuSizing::resetShapes
// - vtuSizing::populateArrays
//
void Foam::vtk::vtuSizing::resetShapes
(
    const UList<cellShape>& shapes
)
{
    const cellModel& tet   = cellModel::ref(cellModel::TET);
    const cellModel& pyr   = cellModel::ref(cellModel::PYR);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& hex   = cellModel::ref(cellModel::HEX);

    decompose_ = false;  // Disallow decomposition
    manifold_ = false;   // Assume no manifold cells possible

    selectionMode_ = selectionModeType::SHAPE_MESH;

    const label nInputCells = shapes.size();

    nCells_    = nInputCells;
    nPoints_   = 0;
    nAddCells_ = 0;
    nAddVerts_ = 0;

    nCellsPoly_  = 0;
    nFacesPoly_  = 0;
    nVertLabels_ = 0;
    nFaceLabels_ = 0;
    nVertPoly_   = 0;

    label nIgnored = 0;

    for (label inputi = 0; inputi < nInputCells; ++inputi)
    {
        const cellShape& shape = shapes[inputi];
        const cellModel& model = shape.model();

        if
        (
            model == tet
         || model == pyr
         || model == prism
         || model == hex
        )
        {
            nVertLabels_ += shape.size();

            // Guess for number of addressed points
            nPoints_ = Foam::max(nPoints_, shape.max());
        }
        else
        {
            --nCells_;
            ++nIgnored;
        }
    }

    if (nIgnored)
    {
        FatalErrorInFunction
            << "Encountered " << nIgnored << " unsupported cell shapes"
            << " ... this is likely not good" << nl
            << exit(FatalError);
    }

    if (nCells_)
    {
        ++nPoints_;
    }
}


Foam::label Foam::vtk::vtuSizing::sizeOf
(
    const enum contentType output,
    const enum slotType slot
) const noexcept
{
    switch (output)
    {
        case contentType::LEGACY:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // legacy uses connectivity for primitives, but directly
                    // stores face streams into connectivity as well.
                    // size-prefix per cell
                    return
                    (
                        nVertLabels() + nAddVerts() - nVertPoly() // primitives
                      + nFaceLabels(output)  // face-stream (poly)
                      + nFieldCells()     // nFieldCells (size prefix)
                    );
                    break;

                // Legacy format does everything via cell connectivity!
                default:
                    break;
            }
            break;
        }

        case contentType::XML:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // Cell connectivity and extra cell centres
                    return (nVertLabels() + nAddVerts());
                    break;

                case slotType::CELLS_OFFSETS:
                    // End offset per cell connectivity
                    return nFieldCells();
                    break;

                case slotType::FACES:
                    // Face stream with embedded sizes ...
                    return nFaceLabels(output);
                    break;

                case slotType::FACES_OFFSETS:
                    // End offset per face connectivity
                    return hasPolyCells() ? nFieldCells() : 0;
                    break;

                // HDF only
                case slotType::POLY_FACEIDS:
                case slotType::POLY_FACEIDS_OFFSETS:
                    break;
            }
            break;
        }

        case contentType::INTERNAL2:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // Cell connectivity and extra cell centres
                    return (nVertLabels() + nAddVerts());
                    break;

                case slotType::CELLS_OFFSETS:
                    // The begin/end offsets for cell connectivity
                    return (nFieldCells() + 1);
                    break;

                case slotType::FACES:
                    // Face stream with various prefixing
                    return nFaceLabels(output);
                    break;

                case slotType::FACES_OFFSETS:
                    // The per-cell begin location of each face stream
                    return hasPolyCells() ? nFieldCells() : 0;
                    break;

                // HDF only
                case slotType::POLY_FACEIDS:
                case slotType::POLY_FACEIDS_OFFSETS:
                    break;
            }
            break;
        }

        case contentType::HDF :
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // Cell connectivity and extra cell centres
                    return (nVertLabels() + nAddVerts());
                    break;

                case slotType::CELLS_OFFSETS:
                    // The begin/end offsets for cell connectivity
                    return (nFieldCells() ? (nFieldCells() + 1) : 0);
                    break;

                case slotType::FACES:
                    // Face vertices (no prefixing!)
                    return nFaceLabels();
                    break;

                case slotType::FACES_OFFSETS:
                    // The begin/end offsets for face connectivity
                    return hasPolyCells() ? (nFacesPoly() + 1) : 0;
                    break;

                case slotType::POLY_FACEIDS:
                    // Polyhedral face lookup (mapping)
                    return hasPolyCells() ? nFacesPoly() : 0;
                    break;

                case slotType::POLY_FACEIDS_OFFSETS:
                    // The per-cell begin/end offsets for polyhedral face lookup
                    // - an empty range for primitive types
                    return hasPolyCells() ? (nFieldCells() + 1) : 0;
                    break;
            }
            break;
        }
    }

    return 0;
}


// * * * * * * * * * * * * * *  Populate Lists * * * * * * * * * * * * * * * //

void Foam::vtk::vtuSizing::populateLegacy
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    labelUList& vertLabels,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        vertLabels,
        unused,  // offsets
        unused,  // faces
        unused,  // facesOffsets
        unused,  // polyFaceIds
        unused,  // polyFaceOffsets
        contentType::LEGACY,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateShapesLegacy
(
    const UList<cellShape>& shapes,
    UList<uint8_t>& cellTypes,
    labelUList& vertLabels,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        shapes,
        *this,
        cellTypes,
        vertLabels,
        unused,  // offsets
        unused,  // faces
        unused,  // facesOffsets
        unused,  // polyFaceIds
        unused,  // polyFaceOffsets
        contentType::LEGACY,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateXml
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    labelUList& connectivity,
    labelUList& offsets,
    labelUList& faces,
    labelUList& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        unused,  // polyFaceIds
        unused,  // polyFaceOffsets
        contentType::XML,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateShapesXml
(
    const UList<cellShape>& shapes,
    UList<uint8_t>& cellTypes,
    labelUList& connectivity,
    labelUList& offsets,
    labelUList& faces,
    labelUList& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        shapes,
        *this,
        cellTypes,
        connectivity,
        offsets,
        unused,  // faces
        unused,  // facesOffsets
        unused,  // polyFaceIds
        unused,  // polyFaceOffsets
        contentType::XML,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateHdf
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    labelUList& connectivity,
    labelUList& offsets,
    labelUList& faces,
    labelUList& facesOffsets,
    labelUList& polyFaceIds,
    labelUList& polyFaceOffsets,
    foamVtkMeshMaps& maps
) const
{
    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        polyFaceIds,
        polyFaceOffsets,
        contentType::HDF,
        maps.cellMap(),
        maps.additionalIds()
    );
}


#undef  definePopulateInternalMethod
#define definePopulateInternalMethod(Type)                                   \
                                                                             \
    void Foam::vtk::vtuSizing::populateInternal                              \
    (                                                                        \
        const polyMesh& mesh,                                                \
        UList<uint8_t>& cellTypes,                                           \
        UList<Type>& connectivity,                                           \
        UList<Type>& offsets,                                                \
        UList<Type>& faces,                                                  \
        UList<Type>& facesOffsets,                                           \
        foamVtkMeshMaps& maps,                                               \
        const enum contentType output                                        \
    ) const                                                                  \
    {                                                                        \
        List<Type> unused;                                                   \
        presizeMaps(maps);                                                   \
                                                                             \
        populateArrays                                                       \
        (                                                                    \
            mesh,                                                            \
            *this,                                                           \
            cellTypes,                                                       \
            connectivity,                                                    \
            offsets,                                                         \
            faces,                                                           \
            facesOffsets,                                                    \
            unused,  /* polyFaceIds */                                       \
            unused,  /* polyFaceOffsets */                                   \
            output,                                                          \
            maps.cellMap(),                                                  \
            maps.additionalIds()                                             \
        );                                                                   \
    }                                                                        \
                                                                             \
    void Foam::vtk::vtuSizing::populateInternal                              \
    (                                                                        \
        const polyMesh& mesh,                                                \
        UList<uint8_t>& cellTypes,                                           \
        UList<Type>& connectivity,                                           \
        UList<Type>& offsets,                                                \
        UList<Type>& faces,                                                  \
        UList<Type>& facesOffsets,                                           \
        labelUList& cellMap,                                                 \
        labelUList& addPointsIds,                                            \
        const enum contentType output                                        \
    ) const                                                                  \
    {                                                                        \
        List<Type> unused;                                                   \
                                                                             \
        populateArrays                                                       \
        (                                                                    \
            mesh,                                                            \
            *this,                                                           \
            cellTypes,                                                       \
            connectivity,                                                    \
            offsets,                                                         \
            faces,                                                           \
            facesOffsets,                                                    \
            unused,  /* polyFaceIds */                                       \
            unused,  /* polyFaceOffsets */                                   \
            output,                                                          \
            cellMap,                                                         \
            addPointsIds                                                     \
        );                                                                   \
    }


definePopulateInternalMethod(int);
definePopulateInternalMethod(long);
definePopulateInternalMethod(long long);


#undef definePopulateInternalMethod

// * * * * * * * * * * * * * * * Renumbering * * * * * * * * * * * * * * * * //

namespace
{

template<class IntListType, class OffsetIntType>
void renumberVerts_legacy
(
    IntListType& vertLabels,
    const OffsetIntType pointOffset
)
{
    if (vertLabels.empty() || !pointOffset)
    {
        return;
    }

    // LEGACY vertLabels = "cells" contains
    // - connectivity
    // [nLabels, vertex labels...]
    // - face-stream
    // [nLabels nFaces, nFace0Pts, id0,id1,..., nFace1Pts, id0,...]

    // Note the simplest volume cell is a tet (4 points, 4 faces)
    // As a poly-face stream this would have
    // 2 for nLabels, nFaces
    // 4 labels (size + ids) per face * 4 == 16 labels
    //
    // Therefore anything with 18 labels or more must be a poly

    auto iter = vertLabels.begin();
    const auto last = vertLabels.end();

    while (iter < last)
    {
        auto nLabels = *iter;  // nLabels (for this cell)
        ++iter;

        if (nLabels < 18)
        {
            // Normal primitive type

            while (nLabels--)
            {
                *iter += pointOffset;
                ++iter;
            }
        }
        else
        {
            // Polyhedral face-stream

            auto nFaces = *iter;
            ++iter;

            while (nFaces--)
            {
                nLabels = *iter;  // nLabels (for this face)
                ++iter;

                while (nLabels--)
                {
                    *iter += pointOffset;
                    ++iter;
                }
            }
        }
    }
}


template<class IntListType, class IntListType2>
void renumberVertLabels_legacy
(
    IntListType& vertLabels,
    const IntListType2& oldToNew
)
{
    if (vertLabels.empty() || oldToNew.empty())
    {
        return;
    }

    // LEGACY vertLabels = "cells" contains
    // - connectivity
    // [nLabels, vertex labels...]
    // - face-stream
    // [nLabels nFaces, nFace0Pts, id0,id1,..., nFace1Pts, id0,...]

    // Note the simplest volume cell is a tet (4 points, 4 faces)
    // As a poly-face stream this would have
    // 2 for nLabels, nFaces
    // 4 labels (size + ids) per face * 4 == 16 labels
    //
    // Therefore anything with 18 labels or more must be a poly

    auto iter = vertLabels.begin();
    const auto last = vertLabels.end();

    while (iter < last)
    {
        auto nLabels = *iter;  // nLabels (for this cell)
        ++iter;

        if (nLabels < 18)
        {
            // Normal primitive type

            while (nLabels--)
            {
                *iter = oldToNew[*iter];
                ++iter;
            }
        }
        else
        {
            // Polyhedral face-stream

            auto nFaces = *iter;
            ++iter;

            while (nFaces--)
            {
                nLabels = *iter;  // nLabels (for this face)
                ++iter;

                while (nLabels--)
                {
                    *iter = oldToNew[*iter];
                    ++iter;
                }
            }
        }
    }
}

} // End anonymous namespace


void Foam::vtk::vtuSizing::renumberVertLabels
(
    labelUList& vertLabels,
    const label pointOffset,
    const contentType output
)
{
    if (vertLabels.empty() || pointOffset <= 0)
    {
        return;
    }

    switch (output)
    {
        case contentType::LEGACY :
        {
            renumberVerts_legacy(vertLabels, pointOffset);
            break;
        }

        default :
        {
            // XML, INTERNAL2, HDF etc
            // vertLabels = "connectivity" contains
            // [cell1-verts, cell2-verts, ...]

            for (auto& id : vertLabels)
            {
                id += pointOffset;
            }
            break;
        }
    }
}


void Foam::vtk::vtuSizing::renumberVertLabels
(
    labelUList& vertLabels,
    const labelUList& oldToNew,
    const contentType output
)
{
    if (vertLabels.empty() || oldToNew.empty())
    {
        return;
    }

    switch (output)
    {
        case contentType::LEGACY :
        {
            renumberVertLabels_legacy(vertLabels, oldToNew);
            break;
        }

        default :
        {
            // XML, INTERNAL2, HDF etc
            // vertLabels = "connectivity" contains
            // [cell1-verts, cell2-verts, ...]

            for (auto& id : vertLabels)
            {
                id = oldToNew[id];
            }
            break;
        }
    }
}


void Foam::vtk::vtuSizing::renumberFaceLabels
(
    labelUList& faceLabels,
    const label pointOffset,
    const contentType output
)
{
    if (faceLabels.empty() || !pointOffset)
    {
        return;
    }

    switch (output)
    {
        case contentType::LEGACY :
        {
            // Not applicable
            break;
        }

        case contentType::HDF :
        {
            // HDF faces
            // [id1,id2,..., id1,id2,...]

            for (auto& id : faceLabels)
            {
                id += pointOffset;
            }
            break;
        }

        default :
        {
            // XML, INTERNAL face-stream
            // [nFaces, nFace0Pts, id1,id2,..., nFace1Pts, id1,id2,...]

            auto iter = faceLabels.begin();
            const auto last = faceLabels.end();

            while (iter < last)
            {
                auto nFaces = *iter;
                ++iter;

                while (nFaces--)
                {
                    auto nLabels = *iter;
                    ++iter;

                    while (nLabels--)
                    {
                        *iter += pointOffset;
                        ++iter;
                    }
                }
            }
            break;
        }
    }
}


void Foam::vtk::vtuSizing::renumberFaceLabels
(
    labelUList& faceLabels,
    const labelUList& oldToNew,
    const contentType output
)
{
    if (faceLabels.empty() || oldToNew.empty())
    {
        return;
    }

    switch (output)
    {
        case contentType::LEGACY :
        {
            // Not applicable
            break;
        }

        case contentType::HDF :
        {
            // HDF faces
            // [id1,id2,..., id1,id2,...]

            for (auto& id : faceLabels)
            {
                id = oldToNew[id];
            }
            break;
        }

        default :
        {
            // XML, INTERNAL face-stream
            // [nFaces, nFace0Pts, id1,id2,..., nFace1Pts, id1,id2,...]

            auto iter = faceLabels.begin();
            const auto last = faceLabels.end();

            while (iter < last)
            {
                auto nFaces = *iter;
                ++iter;

                while (nFaces--)
                {
                    auto nLabels = *iter;
                    ++iter;

                    while (nLabels--)
                    {
                        *iter = oldToNew[*iter];
                        ++iter;
                    }
                }
            }
            break;
        }
    }
}


void Foam::vtk::vtuSizing::renumberFaceOffsets
(
    labelUList& faceOffsets,
    const label beginOffset,
    const contentType output
)
{
    if (faceOffsets.empty() || !beginOffset)
    {
        return;
    }

    switch (output)
    {
        case contentType::LEGACY :
        {
            // Not applicable
            break;
        }

        default :
        {
            // XML, INTERNAL2, etc offsets
            // [-1, off1, off2, ... -1, ..]

            // HDF offsets
            // [off1, off2, ...]

            for (auto& off : faceOffsets)
            {
                // Leave -1 placeholders untouched
                if (off >= 0)
                {
                    off += beginOffset;
                }
            }
            break;
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

void Foam::vtk::vtuSizing::info(Ostream& os) const
{
    os  << "nFieldCells:" << nFieldCells();
    if (nAddCells_)
    {
        os  << " (" << nCells_ << "+" << nAddCells_ << ")";
    }
    else
    {
        os  << " (poly:" << nCellsPoly_ << ")";
    }

    os  << " nFieldPoints:" << nFieldPoints();
    if (nAddPoints_)
    {
        os  << " (" << nPoints_ << "+" << nAddPoints_ << ")";
    }

    os  << " nVertLabels:" << (nVertLabels_ + nAddVerts_);
    if (nAddVerts_)
    {
        os  << " (" << nVertLabels_ << "+" << nAddVerts_ << ")";
    }
    else if (nVertPoly_)
    {
        os  << " (poly:" << nVertPoly_ << ")";
    }

    os << " nFaceLabels:" << nFaceLabels_;
    os << " legacy-count:" << sizeOf<slotType::CELLS>(contentType::LEGACY);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::vtk::vtuSizing::operator==(const vtuSizing& rhs) const
{
    return
    (
        decompose()   == rhs.decompose()
        // required?  && pointOffset() == rhs.pointOffset()
     && nCells()      == rhs.nCells()
     && nPoints()     == rhs.nPoints()
     && nVertLabels() == rhs.nVertLabels()
     && nFaceLabels() == rhs.nFaceLabels()
     && nCellsPoly()  == rhs.nCellsPoly()
     && nFacesPoly()  == rhs.nFacesPoly()
     && nVertPoly()   == rhs.nVertPoly()
     && nAddCells()   == rhs.nAddCells()
     && nAddPoints()  == rhs.nAddPoints()
     && nAddVerts()   == rhs.nAddVerts()
    );
}


bool Foam::vtk::vtuSizing::operator!=(const vtuSizing& rhs) const
{
    return !operator==(rhs);
}


// ************************************************************************* //
