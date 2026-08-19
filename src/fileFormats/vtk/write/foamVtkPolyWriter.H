/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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
    Foam::vtk::polyWriter

Description
    Write faces/points (optionally with fields)
    as a vtp file or a legacy vtk file.
    The class is also slightly expanded to allow reuse as a line-writer,
    but cannot be used to write both faces and lines at once.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkPolyWriter.cxx
    foamVtkPolyWriter.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_polyWriter_H
#define Foam_vtk_polyWriter_H

#include "foamVtkFileWriter.H"
#include "pointField.H"
#include "edgeList.H"
#include "faceList.H"
#include "labelledTri.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                       Class vtk::polyWriter Declaration
\*---------------------------------------------------------------------------*/

class polyWriter
:
    public vtk::fileWriter
{
protected:

    // Protected Data

        //- Slab addressing for field points of the current Piece
        OffsetRange<label> pointSlab_;

        //- Slab addressing for field cells (edges or faces)
        //- of the current Piece
        OffsetRange<label> cellSlab_;


    // Protected Member Functions

        //- Total number of field points for the current Piece
        label nTotalPoints() const noexcept { return pointSlab_.total(); }

        //- Total number of field cells (edges or faces) for the current Piece
        label nTotalCells() const noexcept { return cellSlab_.total(); }

        //- Write verts (legacy or non-legacy format).
        void writeVerts(const label nTotalVerts);

private:

    // Private Member Functions

        //- Determine sizes and begin piece (CellData == lines)
        void beginPiece(const pointField& points, const UList<edge>& edges);

        //- Determine sizes and begin piece (CellData == polys)
        template<class FaceType>
        void beginPiece_impl(const pointField&, const UList<FaceType>&);

        //- Determine sizes and begin piece (CellData == verts)
        void beginPiece(const pointField& points, const bool useVerts);

        //- Write points
        void writePoints(const pointField& points);

        //- Write lines, legacy format
        void writeLines_legacy(const UList<edge>& edges);

        //- Write lines
        void writeLines(const UList<edge>& edges);

        //- Write faces, legacy format
        template<class FaceType>
        void writePolys_legacy_impl(const UList<FaceType>& faces);

        //- Write faces
        template<class FaceType>
        void writePolys_impl(const UList<FaceType>& faces);

        //- Low-level write face/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        template<class FaceType>
        bool writePolyGeometry_impl
        (
            const pointField& points,
            const UList<FaceType>& faces
        );


public:

    // Generated Methods

        //- No copy construct
        polyWriter(const polyWriter&) = delete;

        //- No copy assignment
        void operator=(const polyWriter&) = delete;


    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        explicit polyWriter
        (
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        explicit polyWriter
        (
            const fileName& file,
            bool parallel = UPstream::parRun()
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        polyWriter
        (
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = UPstream::parRun()
        );


    //- Destructor
    virtual ~polyWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type
        static word ext(vtk::outputOptions opts)
        {
            return opts.ext(vtk::fileTag::POLY_DATA);
        }

        //- Dummy write mesh topology method - Fatal if called
        //  Provided for special-purpose writers that do not hold
        //  their own geometry, but use writePolyGeometry() directly
        virtual bool writeGeometry();

        //- Low-level write edge/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        bool writeLineGeometry
        (
            const pointField& points,
            const UList<edge>& edges
        );

        //- Low-level write face/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        bool writePolyGeometry
        (
            const pointField& points,
            const UList<face>& faces
        );

        //- Low-level write face/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        bool writePolyGeometry
        (
            const pointField& points,
            const UList<triFace>& faces
        );

        //- Low-level write face/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        bool writePolyGeometry
        (
            const pointField& points,
            const UList<labelledTri>& faces
        );

        //- Low-level write vertex/point topology.
        //- Normally used by writeGeometry() in a derived class
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        bool writeVertGeometry(const pointField& points);

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


    // Write

        //- In parallel, write "procID" field as CellData or PointData
        //- (depending on the current context). A no-op in serial.
        bool writeProcIDs();

        //- Write rank-local identity map as CellData or PointData
        //- (depending on the current context).
        void writeLocalIDs(const word& fieldName);

        //- Write global identity map as CellData or PointData
        //- (depending on the current context).
        void writeGlobalIDs(const word& fieldName);

        //- Write a uniform field of CellData or PointData
        //- (depending on the current context).
        template<class Type>
        void writeUniform(const word& fieldName, const Type& val);

        //- Write primitive field of CellData or PointData values
        //- (size depending on the current context).
        template<class Type>
        void write(const word& fieldName, const UList<Type>& field);

        //- Write primitive field of CellData
        template<class Type>
        void writeCellData(const word& fieldName, const UList<Type>& field);

        //- Write primitive field of PointData
        template<class Type>
        void writePointData(const word& fieldName, const UList<Type>& field);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkPolyWriter.txx"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
