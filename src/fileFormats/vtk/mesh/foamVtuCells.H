/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2025 OpenCFD Ltd.
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

Class
    Foam::vtk::vtuCells

Description
    A deep-copy description of an OpenFOAM volume mesh in data structures
    suitable for VTK UnstructuredGrid, including the possibility of
    decomposing polyhedral cells into primitive cell types.

    Knowledge of the vtkUnstructuredGrid and the corresponding \c .vtu
    xml file-format aids in understanding this class.
    The class can be used for the VTK xml format, legacy format, as well as a
    VTK internal representation. The internal representation is somewhat
    related to the xml format, but not entirely.

SeeAlso
    Foam::vtk::vtuSizing

SourceFiles
    foamVtuCells.cxx
    foamVtuCells.txx
    foamVtuCellsI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtuCells_H
#define Foam_vtk_vtuCells_H

#include "foamVtkCore.H"
#include "foamVtkMeshMaps.H"
#include "foamVtuSizing.H"
#include "labelList.H"
#include "pointField.H"
#include "refPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

namespace vtk
{

// Forward Declarations
class outputOptions;

/*---------------------------------------------------------------------------*\
                     Class Foam::vtk::vtuCells Declaration
\*---------------------------------------------------------------------------*/

class vtuCells
:
    public vtk::vtuSizing
{
    // Private Data

    // Requested output types

        //- Output content type (enum)
        vtk::vtuSizing::contentType output_;

        //- Bookkeeping for polyhedral cell decomposition
        bool decomposeRequest_;

        //- Bookkeeping for point merge status
        vtk::vtuSizing::mergeType merged_;


    // Storage of output

        //- Cell types (including added cells) in vtk numbering
        //  Range is 1-255
        List<uint8_t> cellTypes_;

        //- Vertices per cell (including added cells) in vtk ordering
        List<label> vertLabels_;

        //- Connectivity (vertices) offsets for each cell.
        //- Can be begin, or end or begin/end (depending on the format)
        List<label> vertOffset_;

        //- Face lists per polyhedral cell
        List<label> faceLabels_;

        //- Face label offsets.
        //- Can be begin, or end or begin/end (depending on the format)
        List<label> faceOffset_;

        //- Face ids per polyhedral cell (HDF only).
        //  Possibly demand-driven content
        mutable List<label> polyFaceIds_;

        //- Offsets into the face ids (HDF only)
        List<label> polyFaceOffset_;

        //- Bookkeeping for polyhedral cell decomposition,
        //- also used by the manual SUBSET_MESH selection mode.
        foamVtkMeshMaps maps_;


    // Private Member Functions

        //- Resize storage
        void resize_all();

        //- Create the geometry entries
        void populateOutput(const polyMesh& mesh);

        //- Create the geometry entries from shapes
        void populateOutput(const UList<cellShape>& shapes);

public:

    // Generated Methods

        //- No copy construct
        vtuCells(const vtuCells&) = delete;

        //- No copy assignment
        void operator=(const vtuCells&) = delete;


    // Constructors

        //- Default construct (XML format, no polyhedral decomposition)
        explicit vtuCells
        (
            const enum contentType output = contentType::XML,
            const bool decompose = false
        );

        //- Construct with output options,
        //- optionally with polyhedral decomposition.
        explicit vtuCells
        (
            const vtk::outputOptions opts,
            const bool decompose = false
        );

        //- Construct from components, create output information immediately
        explicit vtuCells
        (
            const polyMesh& mesh,
            const enum contentType output = contentType::XML,
            const bool decompose = false
        );

        //- Construct from components, create output information immediately
        vtuCells
        (
            const polyMesh& mesh,
            const vtk::outputOptions opts,
            const bool decompose = false
        );

        //- [ADVANCED USAGE]
        //- Construct from primitive shapes,
        //- create output information immediately.
        //  No polyhedral decomposition.
        explicit vtuCells
        (
            const UList<cellShape>& shapes,
            const enum contentType output = contentType::XML
        );

        //- [ADVANCED USAGE]
        //- Construct from primitive shapes,
        //- create output information immediately.
        //  No polyhedral decomposition.
        explicit vtuCells
        (
            const UList<cellShape>& shapes,
            const vtk::outputOptions opts
        );


    // Member Functions

    // Access

        //- The output content type
        inline enum contentType content() const noexcept;

        //- Query the polyhedral decompose requested flag.
        inline bool decomposeRequested() const noexcept;

        //- Query the point merge status.
        inline bool merged() const noexcept;

        //- True if no cellTypes are populated.
        inline bool empty() const noexcept;

        //- The size of populated cellTypes (the number of cells)
        inline label size() const noexcept;

        //- The cellMap is non-identity for a SUBSET_MESH or
        //- when there are decomposed cells.
        inline bool useCellMap() const noexcept;

        //- The pointMap is available and non-identity (eg, after point merge)
        inline bool usePointMap() const noexcept;

        //- Test for hdf content type
        inline bool is_hdf() const noexcept;


    // Edit

        //- Reset all sizes to zero.
        void clear();

        //- Create the geometry using the previously requested output and
        //- decomposition types.
        void reset(const polyMesh& mesh);

        //- Create the geometry for a mesh subset,
        //- using previously requested output and decomposition types.
        void reset
        (
            const polyMesh& mesh,
            const labelUList& subsetCellsIds
        );

        //- Respecify requested output and decomposition type prior to
        //- creating the geometry
        void reset
        (
            const polyMesh& mesh,
            const enum contentType output,
            const bool decompose
        );

        //- Reset sizing using primitive shapes only (ADVANCED USAGE)
        //  Effectively removes any polyhedrals!
        void resetShapes(const UList<cellShape>& shapes);

        //- Renumber cell ids to account for subset meshes
        void renumberCells(const labelUList& mapping);

        //- Renumber point ids to account for subset meshes
        void renumberPoints(const labelUList& mapping);

        //- Define which additional cell-centres are to be used (ADVANCED!)
        void addPointCellLabels(const labelUList& cellIds);

        //- Determine topological point merging, renumbers contents accordingly
        void mergePoints(const polyMesh& mesh);


    // Storage Access (legacy and non-legacy)

        //- The cell \c "types" (legacy: \c "CELL_TYPES")
        inline const List<uint8_t>& cellTypes() const noexcept;

        //- The cell \c "connectivity" (legacy: \c "CELLS")
        inline const labelList& vertLabels() const noexcept;


    // Storage Access (non-legacy only)

        //- The \c "offsets" into the connectivity list vertLabels()
        //  The content depends on the format:
        //  - XML : end offset
        //  - INTERNAL2 : begin/end offsets
        //  - HDF : begin/end offsets
        inline const labelList& vertOffsets() const noexcept;

        //- The polyhedral \c "faces" or \c "FaceConnectivity".
        //  The content depends on the format:
        //  - XML / INTERNAL : a face stream with embedded sizing
        //  - HDF : the vertex list of polyhedral faces
        inline const labelList& faceLabels() const noexcept;

        //- The \c "facesoffset" or the \c "FaceOffsets".
        //  The content depends on the format:
        //  - XML : the per-cell \c end offset into faceLabels(),
        //  - INTERNAL : the per-cell \c begin offset into faceLabels(),
        //    with -1 for primitive types
        //  - HDF : the per-cell begin/end offset into faceLabels(),
        //    with an empty range for primitive types
        inline const labelList& faceOffsets() const noexcept;

        //- Face ids per polyhedral cell (HDF only)
        //  Possibly demand-driven content
        const labelList& polyFaceIds() const;

        //- Offsets into the polyhedral face ids (HDF only)
        inline const labelList& polyFaceOffsets() const noexcept;

        //- Additional point addressing (from added point to original cell)
        inline const labelList& addPointCellLabels() const noexcept;

        //- The mesh cell ids for all cells (regular and decomposed).
        //  For a volume field these would be the cell ids to extract
        //  from the field.
        inline const labelList& cellMap() const noexcept;

        //- The mesh point ids [FUTURE USE]
        inline const labelUList& pointMap() const noexcept;


    // Adjusted versions

        //- The "connectivity" (legacy: "CELLS")
        //- adjusted by the specified (global) point offset
        refPtr<labelList> vertLabels(label pointOffset) const;

        //- Connectivity (vertices) offsets for each cell
        //- adjusted by the specified (global) cell offset
        refPtr<labelList> vertOffsets
        (
            //! For formats with begin/end, can specify -1 for
            //! an automatic calculation (in parallel)
            label beginOffset,
            //! For begin/end offsets for parallel concatenation
            bool syncPar = false
        ) const;

        //- The "faces" stream/connectivity (non-legacy only)
        //- adjusted by the specified (global) point offset
        refPtr<labelList> faceLabels(label pointOffset) const;

        //- The offsets into "faces" (non-legacy only)
        //- adjusted by the specified (global) face offset.
        //  Returns a dummy set of offsets if there are no "faces"
        refPtr<labelList> faceOffsets
        (
            //! For formats with begin/end, can specify -1 for
            //! an automatic calculation (in parallel)
            label beginOffset,
            //! For begin/end offsets for parallel concatenation
            bool syncPar = false
        ) const;

        //- The "PolyhedronToFaces" ids (HDF only)
        //- adjusted by the specified (global) face offset.
        refPtr<labelList> polyFaceIds(label beginOffset) const;

        //- The "PolyhedronOffsets" (HDF only), which are the
        //- offsets into the polyhedral face ids,
        //- adjusted by the specified (global) cell offset.
        refPtr<labelList> polyFaceOffsets
        (
            //! For formats with begin/end, can specify -1 for
            //! an automatic calculation (in parallel)
            label beginOffset,
            //! For begin/end offsets for parallel concatenation
            bool syncPar = false
        ) const;


    // Mapping/Generation Methods

        //- Return the mesh points,
        //- possibly with cell centres for addPointCellLabels()
        //  \note not valid for a SHAPE_MESH
        template<class PointType = Foam::point>
        tmp<Field<PointType>> points(const polyMesh& mesh) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtuCellsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtuCells.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
