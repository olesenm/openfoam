/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
    Foam::vtk::vtuWriter

Description
    Write an OpenFOAM volume (internal) geometry and internal fields
    as a vtu file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

See Also
    Foam::vtk::internalMeshWriter
    Foam::vtk::internalWriter

SourceFiles
    foamVtuWriter.cxx
    foamVtuWriter.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtuWriter_H
#define Foam_vtk_vtuWriter_H

#include "foamVtkFileWriter.H"
#include "foamVtuCells.H"
#include "pointField.H"
#include "refPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                      Class vtk::vtuWriter Declaration
\*---------------------------------------------------------------------------*/

class vtuWriter
:
    public vtk::fileWriter
{
protected:

    // Protected Data

        //- Reference to the points (eg, for cell shapes)
        refPtr<pointField> points_;

        //- The volume cells (internalMesh)
        const vtk::vtuCells& vtuCells_;

        //- Slab addressing for field points of the current Piece
        OffsetRange<label> pointSlab_;

        //- Slab addressing for field cells of the current Piece
        OffsetRange<label> cellSlab_;


    // Protected Member Functions

        //- Total number of field points for the current Piece
        label nTotalPoints() const noexcept { return pointSlab_.total(); }

        //- Total number of field cells for the current Piece
        label nTotalCells() const noexcept { return cellSlab_.total(); }

        //- Write mesh topology (implementation).
        bool writeGeometry_impl
        (
            const UList<point>& points,
            const UList<point>& cellCentres  //< for decomposed poly
        );


    // Constructors

        //- Construct without point reference (default format INLINE_BASE64),
        explicit vtuWriter
        (
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

private:

    // Private Member Functions

        //- Construction helper
        void setPointReference(const pointField& points);

        //- Begin piece
        void beginPiece();

        //- Write mesh points
        void writePoints
        (
            const UList<point>& points,
            const UList<point>& cellCentres  //< for decomposed poly
        );

        //- Write cells (connectivity and type), legacy format
        void writeCells_legacy();

        //- Write cells connectivity
        void writeCellsConnectivity();

        //- Write cells face streams
        void writeCellsFaces();

public:

    //- Debug information
    static int debug;

//- File extension for given output type
//  constexpr static vtk::fileTag content = vtk::fileTag::UNSTRUCTURED_GRID;

    // Generated Methods

        //- No copy construct
        vtuWriter(const vtuWriter&) = delete;

        //- No copy assignment
        void operator=(const vtuWriter&) = delete;


    // Constructors

        //- Construct from components (default format INLINE_BASE64),
        vtuWriter
        (
            const pointField& points,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        vtuWriter
        (
            const pointField& points,
            const vtk::vtuCells& cells,
            const fileName& file,
            bool parallel = UPstream::parRun()
        );

        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        vtuWriter
        (
            const pointField& points,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = UPstream::parRun()
        );


    //- Destructor
    virtual ~vtuWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type
        inline static word ext(vtk::outputOptions opts)
        {
            return opts.ext(vtk::fileTag::UNSTRUCTURED_GRID);
        }

        //- Write mesh topology.
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry();

        //- Begin CellData output section for specified number of fields.
        //  Must be called prior to writing any cell data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | POINT_DATA).
        //
        //  \return True if the state changed
        virtual bool beginCellData(label nFields = 0);

        //- Begin PointData for specified number of fields.
        //  Must be called prior to writing any point data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | CELL_DATA).
        //
        //  \return True if the state changed
        virtual bool beginPointData(label nFields = 0);


        //- Write cell ids as CellData.
        //  Must be called within the CELL_DATA state.
        void writeCellIDs();

        //- Write processor ids as CellData. This is no-op in serial.
        //  Must be called within the CELL_DATA state.
        bool writeProcIDs();

        //- Write point ids as PointData.
        //  Must be called within the POINT_DATA state.
        void writePointIDs();


    // Write Fields

        //- Write a uniform field of Cell or Point values
        template<class Type>
        void writeUniform(const word& fieldName, const Type& val);

        //- Write primitive field of CellData.
        //- Uses cellMap (for decomposed poly cells).
        template<class Type>
        void writeCellData(const word& fieldName, const UList<Type>& field);

        //- Write primitive field of PointData.
        //- Uses pointMap (for merged data).
        template<class Type>
        void writePointData(const word& fieldName, const UList<Type>& field);

        //- Write primitive field of PointData.
        //- Uses pointMap (for merged data).
        template<class Type>
        void writePointData
        (
            const word& fieldName,
            const UList<Type>& field,
            //! Extra values corresponding to the decomposed cell centres.
            //! Ignored if UList<Type>::null()
            const UList<Type>& extra
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtuWriter.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
