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

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Adjust \p vertOffset for all cells.
// On input it contains the cell sizes.
//
// On output, the cell connectivity offsets :
// - XML format = end-offsets
// - INTERNAL2 = begin/end-offsets
// - HDF = begin/end-offsets
// .
template<class IntType>
void adjustCellOffsets
(
    Foam::UList<IntType>& vertOffset,
    const Foam::vtk::vtuSizing::contentType output
)
{
    using namespace Foam;
    using namespace Foam::vtk;

    switch (output)
    {
        case vtuSizing::contentType::LEGACY :  // Nothing to do
            break;

        case vtuSizing::contentType::XML :
        {
            // Transform cell sizes (vertOffset) into end offsets

            // vertOffset[0] already contains its size, leave untouched
            for (label i = 1; i < vertOffset.size(); ++i)
            {
                vertOffset[i] += vertOffset[i-1];
            }
            break;
        }

        case vtuSizing::contentType::INTERNAL2 :
        case vtuSizing::contentType::HDF :
        {
            // Transform cell sizes (vertOffset) into begin/end offsets
            // input    [n1, n2, n3, ..., 0]
            // becomes  [0, n1, n1+n2, n1+n2+n3, ..., nTotal]

            if (!vertOffset.empty())
            {
                vertOffset.back() = 0;  // safety (if uninitialized)
                IntType total(0);

                for (IntType& off : vertOffset)
                {
                    const IntType sz(off);
                    off = total;
                    total += sz;
                }
            }
            break;
        }
    }
}


// Adjust \p faceOffset for all polyhedral faces
// On input it contains the face or face-stream sizes,
// possibly with primitive cells marked as -1 placeholders
//
// On output, the face connectivity offsets :
// - XML format = end-offsets, with -1 placeholders
// - INTERNAL2 = begin/end-offsets, with -1 placeholders
// - HDF = begin/end-offsets
// .
template<class IntType>
void adjustFaceOffsets
(
    Foam::UList<IntType>& faceOffset,
    const Foam::vtk::vtuSizing::contentType output
)
{
    using namespace Foam;
    using namespace Foam::vtk;

    switch (output)
    {
        case vtuSizing::contentType::LEGACY : // Nothing to do
            break;

        case vtuSizing::contentType::XML :
        {
            // Transform face sizes (faceOffset) into end offsets,
            // leaving -1 placeholders untouched
            if (!faceOffset.empty())
            {
                IntType total(0);

                for (IntType& off : faceOffset)
                {
                    const IntType sz(off);
                    if (sz >= 0)
                    {
                        total += sz;
                        off = total;
                    }
                }
            }
            break;
        }

        case vtuSizing::contentType::INTERNAL2 :
        {
            // Transform face sizes (faceOffset) into begin locations,
            // leaving -1 placeholders untouched
            if (!faceOffset.empty())
            {
                IntType beg(0);

                for (IntType& off : faceOffset)
                {
                    const IntType sz(off);
                    if (sz >= 0)
                    {
                        off = beg;
                        beg += sz;
                    }
                }
            }
            break;
        }

        case vtuSizing::contentType::HDF :
        {
            // Transform face sizes (faceOffset) into begin/end offsets,
            // treat any -1 placeholders like size = 0

            if (!faceOffset.empty())
            {
                faceOffset.back() = 0;  // safety (if uninitialized)
                IntType total(0);

                for (IntType& off : faceOffset)
                {
                    const IntType sz(off);
                    off = total;
                    if (sz >= 0)
                    {
                        total += sz;
                    }
                }
            }
            break;
        }
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class LabelType>
void Foam::vtk::vtuSizing::populateArrays
(
    const polyMesh& mesh,
    const vtk::vtuSizing& sizing,

    UList<uint8_t>& cellTypes,
    UList<LabelType>& vertLabels,
    UList<LabelType>& vertOffset,
    UList<LabelType>& faceLabels,
    UList<LabelType>& faceOffset,

    UList<LabelType>& polyFaceIds,      // HDF-only
    UList<LabelType>& polyFaceOffsets,  // HDF-only

    const enum contentType output,
    labelUList& cellMap,
    labelUList& addPointsIds
)
{
    if (sizing.selectionMode() == selectionModeType::SHAPE_MESH)
    {
        FatalErrorInFunction
            << "Programming error ... attempting to populate a VTU mesh"
            << " but it was originally sized using independent cell shapes"
            << exit(FatalError);
    }

    // Verify storage sizes
    checkSizes
    (
        sizing,

        cellTypes.size(),
        vertLabels.size(), vertOffset.size(),
        faceLabels.size(), faceOffset.size(),
        polyFaceIds.size(), polyFaceOffsets.size(),

        output,

        cellMap.size(),
        addPointsIds.size()
    );

    // Characteristics

    // Are vertLabels prefixed with the size?
    // Also use as the size of the prefixed information
    const int prefix = (output == contentType::LEGACY) ? 1 : 0;


    // Initialization
    // ~~~~~~~~~~~~~~

    // The "HDF" mode : uses face vertices/offsets, not a face-stream
    const bool isHDFmode = (!polyFaceOffsets.empty());

    if (isHDFmode)
    {
        faceOffset = 0;
    }
    else
    {
        // For XML + INTERNAL formats, tag with -1 for primitives
        faceOffset = -1;
    }

    // Some formats have begin/end offsets (eg, nFieldCells+1),
    // which means that either end may be unvisited. Set as zero now.

    if (vertOffset.size())
    {
        vertOffset.front() = 0;
        vertOffset.back() = 0;
    }

    // Looks strange, but we only generate HDF faces by walking
    // each polyhedron. So the corresponding face ids are always
    // just an identity map.
    // Thus polyFaceIds can also be treated as demand-driven content
    // elsewhere.

    std::iota(polyFaceIds.begin(), polyFaceIds.end(), 0);
    polyFaceOffsets = 0;
    label nPolyFaces = 0;  // Counter into faceOffset (HDF-mode)

    const cellModel& tet      = cellModel::ref(cellModel::TET);
    const cellModel& pyr      = cellModel::ref(cellModel::PYR);
    const cellModel& prism    = cellModel::ref(cellModel::PRISM);
    const cellModel& hex      = cellModel::ref(cellModel::HEX);
    const cellModel& wedge    = cellModel::ref(cellModel::WEDGE);
    const cellModel& tetWedge = cellModel::ref(cellModel::TETWEDGE);

    const cellShapeList& shapes = mesh.cellShapes();
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();

    // The face owner is needed to determine the face orientation
    const labelList& owner = mesh.faceOwner();


    // Index into vertLabels, faceLabels for normal cells
    label nVertLabels = 0;
    label nFaceLabels = 0;

    // Index into vertLabels for decomposed polys
    label nVertDecomp = sizing.nVertLabels() + prefix*sizing.nCells();

    // Placement of additional decomposed cells
    label nCellDecomp = mesh.nCells();

    // Placement of additional point labels (apex of decomposed cells)
    label apexVertLabel = mesh.nPoints();

    // Count of additional point labels
    label nAddPoints = 0;

    // Non-decomposed polyhedral are represented as a face-stream.
    // For legacy format, this stream replaces the normal connectivity
    // information. Use references to alias where the face output should land.

    UList<LabelType>& faceOutput =
    (
        output == contentType::LEGACY
      ? vertLabels
      : faceLabels
    );

    label& faceIndexer =
    (
        output == contentType::LEGACY
      ? nVertLabels
      : nFaceLabels
    );

    // ===========================================
    // STAGE 2: Rewrite in VTK form
    // During this stage, the vertOffset contains the *size* associated with
    // the per-cell vertLabels entries, and the faceOffset contains the *size*
    // associated with the per-cell faceLabels.
    // Similarly, the polyFaceOffsets will contain the per-cell
    // *number* of polyhedral faces at this stage.


    // Special treatment for mesh subsets
    // Here the cellMap is the list of input cells!

    const bool isSubsetMesh
    (
        sizing.selectionMode() == selectionModeType::SUBSET_MESH
    );

    const label nInputCells =
    (
        isSubsetMesh
      ? cellMap.size()
      : shapes.size()
    );

    // Unique vertex labels per polyhedral
    labelHashSet hashUniqId;
    if (!sizing.decompose()) { hashUniqId.reserve(256); }

    // For polyhedral decomposition, ensure that the face tri/quad splits
    // on the neighbour side of a processor patch follow that of the owner side.
    bitSet flipBoundaryFace;
    if (sizing.decompose() && UPstream::parRun())
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

    for
    (
        label inputi = 0, cellIndex = 0; // cellIndex: the ouput location
        inputi < nInputCells;
        ++inputi, ++cellIndex
    )
    {
        const label celli(isSubsetMesh ? cellMap[inputi] : inputi);

        const cellShape& shape = shapes[celli];
        const cellModel& model = shape.model();

        if (!isSubsetMesh)
        {
            cellMap[cellIndex] = celli;
        }

        if (model == tet)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_TETRA;
            constexpr label nShapePoints = 4; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == pyr)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_PYRAMID;
            constexpr label nShapePoints = 5; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == hex)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_HEXAHEDRON;
            constexpr label nShapePoints = 8; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == prism)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_WEDGE;
            constexpr label nShapePoints = 6; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            // VTK_WEDGE triangles point outwards (swap 1<->2, 4<->5)
            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[4];
        }
        else if (model == tetWedge && sizing.decompose())
        {
            // Treat as squeezed prism
            cellTypes[cellIndex] = vtk::cellType::VTK_WEDGE;
            constexpr label nShapePoints = 6;

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[4];
            vertLabels[nVertLabels++] = shape[3];
        }
        else if (model == wedge && sizing.decompose())
        {
            // Treat as squeezed hex
            cellTypes[cellIndex] = vtk::cellType::VTK_HEXAHEDRON;
            constexpr label nShapePoints = 8;

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[4];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[6];
        }
        else if (sizing.decompose())
        {
            // Polyhedral cell - decompose into tet/pyr.

            // Ensure we have the correct orientation for the base of the
            // primitive cell shape.
            // If the cell is face owner, the orientation needs to be flipped
            // to avoid defining negative cells.
            // VTK may not care, but we'll do it anyhow for safety.

            // The cell-centre corresponding to the apexVertLabel
            addPointsIds[nAddPoints] = celli;

            // Whether to insert cell in place of original or not.
            bool firstCell = true;

            // Count triangles/quads in decomposition
            label nTria = 0, nQuad = 0;
            DynamicList<face> faces3, faces4;

            for (const label facei : meshCells[celli])
            {
                // Face decomposed into triangles and quads
                // Tri -> Tet, Quad -> Pyr
                const face& f = meshFaces[facei];
                bool isOwner = (owner[facei] == celli);

                // Need to use a flipped face on the boundary?
                if
                (
                    (f.size() > 4)
                 && flipBoundaryFace.test(facei - mesh.nInternalFaces())
                )
                {
                    isOwner = !isOwner;

                    const face flipped(f.reverseFace());

                    // Count triangles/quads in decomposition
                    nTria = nQuad = 0;
                    flipped.nTrianglesQuads(mesh.points(), nTria, nQuad);

                    // Do actual decomposition
                    faces3.resize_nocopy(nTria);
                    faces4.resize_nocopy(nQuad);
                    nTria = nQuad = 0;

                    flipped.trianglesQuads
                    (
                        mesh.points(), nTria, nQuad, faces3, faces4
                    );
                }
                else
                {
                    // Count triangles/quads in decomposition
                    nTria = nQuad = 0;
                    f.nTrianglesQuads(mesh.points(), nTria, nQuad);

                    // Do actual decomposition
                    faces3.resize_nocopy(nTria);
                    faces4.resize_nocopy(nQuad);
                    nTria = nQuad = 0;

                    f.trianglesQuads
                    (
                        mesh.points(), nTria, nQuad, faces3, faces4
                    );
                }

                for (const face& quad : faces4)
                {
                    // Quad becomes a pyramid

                    constexpr label nShapePoints = 5;  // pyr (5 vertices)

                    label celLoc, vrtLoc;
                    if (firstCell)
                    {
                        firstCell = false;
                        celLoc = cellIndex;
                        vrtLoc = nVertLabels;
                        nVertLabels += prefix + nShapePoints;
                    }
                    else
                    {
                        celLoc = nCellDecomp++;
                        vrtLoc = nVertDecomp;
                        nVertDecomp += prefix + nShapePoints;
                    }
                    cellMap[celLoc] = celli;

                    cellTypes[celLoc] = vtk::cellType::VTK_PYRAMID;
                    if (vertOffset.size())
                    {
                        vertOffset[celLoc] = nShapePoints;
                    }
                    if (prefix)
                    {
                        vertLabels[vrtLoc++] = nShapePoints;
                    }

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels[vrtLoc++] = quad[0];
                        vertLabels[vrtLoc++] = quad[3];
                        vertLabels[vrtLoc++] = quad[2];
                        vertLabels[vrtLoc++] = quad[1];
                    }
                    else
                    {
                        vertLabels[vrtLoc++] = quad[0];
                        vertLabels[vrtLoc++] = quad[1];
                        vertLabels[vrtLoc++] = quad[2];
                        vertLabels[vrtLoc++] = quad[3];
                    }

                    // The apex
                    vertLabels[vrtLoc++] = apexVertLabel;
                }

                for (const face& tria : faces3)
                {
                    // Triangle becomes a tetrahedral

                    constexpr label nShapePoints = 4;  // tet (4 vertices)

                    label celLoc, vrtLoc;
                    if (firstCell)
                    {
                        firstCell = false;
                        celLoc = cellIndex;
                        vrtLoc = nVertLabels;
                        nVertLabels += prefix + nShapePoints;
                    }
                    else
                    {
                        celLoc = nCellDecomp++;
                        vrtLoc = nVertDecomp;
                        nVertDecomp += prefix + nShapePoints;
                    }
                    cellMap[celLoc] = celli;

                    cellTypes[celLoc] = vtk::cellType::VTK_TETRA;
                    if (vertOffset.size())
                    {
                        vertOffset[celLoc] = nShapePoints;
                    }
                    if (prefix)
                    {
                        vertLabels[vrtLoc++] = nShapePoints;
                    }

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels[vrtLoc++] = tria[0];
                        vertLabels[vrtLoc++] = tria[2];
                        vertLabels[vrtLoc++] = tria[1];
                    }
                    else
                    {
                        vertLabels[vrtLoc++] = tria[0];
                        vertLabels[vrtLoc++] = tria[1];
                        vertLabels[vrtLoc++] = tria[2];
                    }

                    // The apex
                    vertLabels[vrtLoc++] = apexVertLabel;
                }
            }

            // Increment values
            ++apexVertLabel;
            ++nAddPoints;
        }
        else
        {
            // Polyhedral cell - not decomposed
            hashUniqId.clear();  // unique node ids used (XML, INTERNAL)

            // face-stream
            //   [nFaces, nFace0Pts, id0,id1,..., nFace1Pts, id0,...]
            // but HDF only stores the vertices!
            //   [id0,id12, ..., id0,...]

            cellTypes[cellIndex] = vtk::cellType::VTK_POLYHEDRON;

            const labelList& cFaces = meshCells[celli];

            const label startLabel = faceIndexer;

            if (output == contentType::LEGACY)
            {
                faceOutput[startLabel] = 0; // placeholder for total size
                ++faceIndexer;
            }

            if (isHDFmode)
            {
                // With separate handling of poly faces
                // - eg, VTKHDF format
                polyFaceOffsets[cellIndex] = cFaces.size();
            }
            else
            {
                // With embedded handling of number of poly faces
                // - currently all formats other than VTKHDF
                faceOutput[faceIndexer++] = cFaces.size();
            }

            for (const label facei : cFaces)
            {
                const face& f = mesh.faces()[facei];
                const bool isOwner = (owner[facei] == celli);
                const label nFacePoints = f.size();

                hashUniqId.insert(f);

                // The number of labels for this face
                if (isHDFmode)
                {
                    // With separate handling of poly faces
                    // - eg, for VTKHDF format
                    faceOffset[nPolyFaces] = nFacePoints;
                }
                else
                {
                    // With embedded handling of poly faces
                    faceOutput[faceIndexer++] = nFacePoints;
                }
                ++nPolyFaces;

                // The actual face vertices:
                faceOutput[faceIndexer++] = f[0];
                if (isOwner)
                {
                    for (label fp = 1; fp < nFacePoints; ++fp)
                    {
                        faceOutput[faceIndexer++] = f[fp];
                    }
                }
                else
                {
                    for (label fp = nFacePoints - 1; fp > 0; --fp)
                    {
                        faceOutput[faceIndexer++] = f[fp];
                    }
                }
            }

            if (output == contentType::LEGACY)
            {
                // Update size for legacy face stream
                // (subtract 1 to avoid counting the storage location)
                faceOutput[startLabel] = (faceIndexer - 1 - startLabel);
            }
            else
            {
                if (!isHDFmode)
                {
                    // Size of the face stream
                    faceOffset[cellIndex] = (faceIndexer - startLabel);
                }

                vertOffset[cellIndex] = hashUniqId.size();
                if (prefix)
                {
                    vertLabels[nVertLabels++] = hashUniqId.size();
                }

                for (const label pointi : hashUniqId.sortedToc())
                {
                    vertLabels[nVertLabels++] = pointi;
                }
            }
        }
    }

    // ===========================================
    // STAGE 3: Adjust sizes into offsets

    adjustCellOffsets<LabelType>(vertOffset, output);

    if (sizing.hasPolyCells())
    {
        adjustFaceOffsets<LabelType>(faceOffset, output);
        adjustFaceOffsets<LabelType>(polyFaceOffsets, output);
    }
}


// Synchronize changes here with the following:
// - vtuSizing::resetShapes
// - vtuSizing::populateArrays

template<class LabelType>
void Foam::vtk::vtuSizing::populateArrays
(
    const UList<cellShape>& shapes,
    const vtk::vtuSizing& sizing,

    UList<uint8_t>& cellTypes,
    UList<LabelType>& vertLabels,
    UList<LabelType>& vertOffset,
    UList<LabelType>& faceLabels,
    UList<LabelType>& faceOffset,
    UList<LabelType>& polyFaceIds,
    UList<LabelType>& polyFaceOffsets,
    const enum contentType output,
    labelUList& cellMap,
    labelUList& addPointsIds
)
{
    if (sizing.selectionMode() != selectionModeType::SHAPE_MESH)
    {
        FatalErrorInFunction
            << "Programming error ... attempting to populate a VTU mesh"
            << " from cell shapes, but sizing originated from a different"
            << " representation" << nl
            << exit(FatalError);
    }

    // Verify storage sizes
    checkSizes
    (
        sizing,

        cellTypes.size(),
        vertLabels.size(), vertOffset.size(),
        faceLabels.size(), faceOffset.size(),
        polyFaceIds.size(), polyFaceOffsets.size(),

        output,

        cellMap.size(),
        addPointsIds.size()
    );

    // Characteristics

    // Are vertLabels prefixed with the size?
    // Also use as the size of the prefixed information
    const int prefix = (output == contentType::LEGACY) ? 1 : 0;


    // Initialization
    // ~~~~~~~~~~~~~~

    // For XML + INTERNAL formats, tag with -1 for primitives
    faceOffset = -1;

    // Some formats have begin/end offsets (eg, nFieldCells+1),
    // which means that either end may be unvisited. Set as zero now.

    if (vertOffset.size())
    {
        vertOffset.front() = 0;
        vertOffset.back() = 0;
    }

    // Don't actually support polyhedrals, so just tag as "bad"
    polyFaceIds = -1;
    polyFaceOffsets = 0;

    const cellModel& tet      = cellModel::ref(cellModel::TET);
    const cellModel& pyr      = cellModel::ref(cellModel::PYR);
    const cellModel& prism    = cellModel::ref(cellModel::PRISM);
    const cellModel& hex      = cellModel::ref(cellModel::HEX);

    // Index into vertLabels for normal cells
    label nVertLabels = 0;

    // ===========================================
    // STAGE 2: Rewrite in VTK form
    // During this stage, the vertOffset contains the *size* associated with
    // the per-cell vertLabels entries, and the faceOffset contains the *size*
    // associated with the per-cell faceLabels.

    const label nInputCells = shapes.size();

    // label nIgnored = 0;

    for
    (
        label inputi = 0, cellIndex = 0; // cellIndex: the ouput location
        inputi < nInputCells;
        ++inputi, ++cellIndex
    )
    {
        const cellShape& shape = shapes[inputi];
        const cellModel& model = shape.model();

        if (model == tet)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_TETRA;
            constexpr label nShapePoints = 4; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == pyr)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_PYRAMID;
            constexpr label nShapePoints = 5; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == hex)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_HEXAHEDRON;
            constexpr label nShapePoints = 8; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            for (const label cpi : shape)
            {
                vertLabels[nVertLabels++] = cpi;
            }
        }
        else if (model == prism)
        {
            cellTypes[cellIndex] = vtk::cellType::VTK_WEDGE;
            constexpr label nShapePoints = 6; // OR shape.size();

            if (vertOffset.size())
            {
                vertOffset[cellIndex] = nShapePoints;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = nShapePoints;
            }

            // VTK_WEDGE triangles point outwards (swap 1<->2, 4<->5)
            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[4];
        }
        else
        {
            // Silent here.
            // - already complained (and skipped) during initial sizing
            --cellIndex;
            // ++nIgnored;
        }
    }

    // May have been done by caller,
    // but for additional safety set an identity mapping
    std::iota(cellMap.begin(), cellMap.end(), 0);

    // ===========================================
    // Adjust sizes into offsets

    adjustCellOffsets<LabelType>(vertOffset, output);

    // This should be a no-op for shape-based conversions
    if (sizing.hasPolyCells())
    {
        adjustFaceOffsets<LabelType>(faceOffset, output);
        adjustFaceOffsets<LabelType>(polyFaceOffsets, output);
    }
}


//unused template<class LabelType, class LabelType2>
//unused void Foam::vtk::vtuSizing::renumberVertLabelsInternalImpl
//unused (
//unused     UList<uint8_t>& cellTypes,
//unused     UList<LabelType>& vertLabels,
//unused     const LabelType2 globalPointOffset
//unused )
//unused {
//unused     // INTERNAL vertLabels = "connectivity" contain
//unused     // [nLabels, vertex labels...]
//unused
//unused     auto iter = vertLabels.begin();
//unused     const auto last = vertLabels.end();
//unused
//unused     while (iter < last)
//unused     {
//unused         LabelType nLabels = *iter;
//unused         ++iter;
//unused
//unused         while (nLabels--)
//unused         {
//unused             *iter += globalPointOffset;
//unused             ++iter;
//unused         }
//unused     }
//unused }


// ************************************************************************* //
