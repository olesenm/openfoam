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
    Foam::vtk::surfaceWriter

Description
    Write faces/points (optionally with fields)
    as a vtp file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

    This writer uses a fat interface to support face/triFace/labelledTri
    types but these must be used consistently in parallel.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkSurfaceWriter.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_surfaceWriter_H
#define Foam_vtk_surfaceWriter_H

#include "foamVtkPolyWriter.H"
#include "instant.H"
#include "refPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceWriter
:
    public vtk::polyWriter
{
    // Private Member Data

        //- Reference to the points
        refPtr<pointField> points_;

        //- Reference to the (regular) faces
        refPtr<List<face>> faces_;

        //- Reference to the tri-faces
        refPtr<List<triFace>> triFaces_;

        //- Reference to the labelledTri-faces
        refPtr<List<labelledTri>> labelledTris_;

        //- Time name/value
        instant instant_;

        //- Output as 1D vertex/point elements instead of faces
        bool vertexOutput_;


    // Private Member Functions

        //- Construct from components. Internal use only
        template<class FaceType>
        surfaceWriter
        (
            const vtk::outputOptions,
            const pointField& pts,
            const List<FaceType>& faces
        );

        //- Set poly, triFace or labelledTri face reference
        template<class FaceType>
        void setFaceReference(const List<FaceType>& faces);

public:

    // Generated Methods

        //- No copy construct
        surfaceWriter(const surfaceWriter&) = delete;

        //- No copy assignment
        void operator=(const surfaceWriter&) = delete;


    // Constructors

#undef  declare_constructors
#define declare_constructors(FaceType)                                        \
        /*! \brief Construct from components */                               \
        surfaceWriter                                                         \
        (                                                                     \
            const pointField& points,                                         \
            const List<FaceType>& faces,                                      \
            /*! Output format (default: INLINE_BASE64) */                     \
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64    \
        );                                                                    \
                                                                              \
        /*! \brief Construct (INLINE_BASE64) and open file for writing. */    \
        surfaceWriter                                                         \
        (                                                                     \
            const pointField& points,                                         \
            const List<FaceType>& faces,                                      \
            /*! The file name (with or without an extension). */              \
            const fileName& file,                                             \
            bool parallel = UPstream::parRun()                                \
        )                                                                     \
        :                                                                     \
            surfaceWriter(points, faces)                                      \
        {                                                                     \
            open(file, parallel);                                             \
        }                                                                     \
                                                                              \
        /*! \brief Construct with given format and open file for writing. */  \
        surfaceWriter                                                         \
        (                                                                     \
            const pointField& points,                                         \
            const List<FaceType>& faces,                                      \
            /*! Output format */                                              \
            const vtk::outputOptions opts,                                    \
            /*! The file name (with or without an extension). */              \
            const fileName& file,                                             \
            bool parallel = UPstream::parRun()                                \
        )                                                                     \
        :                                                                     \
            surfaceWriter(points, faces, opts)                                \
        {                                                                     \
            open(file, parallel);                                             \
        }

        declare_constructors(face);
        declare_constructors(triFace);
        declare_constructors(labelledTri);

        #undef declare_constructors

    //- Destructor
    virtual ~surfaceWriter() = default;


    // Member Functions

        //- Output as 1D vertex/point elements instead of faces
        bool vertexOutput() const noexcept { return vertexOutput_; }

        //- Output as 1D vertex/point elements instead of faces
        bool vertexOutput(bool on) noexcept
        {
            bool old(vertexOutput_);
            vertexOutput_ = on;
            return old;
        }

        //- Define a time name/value for the output
        virtual void setTime(const instant& inst);

        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED).
        virtual bool beginFile(std::string title = "");

        //- Write faces topology
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry();

        //- Write "TimeValue" FieldData (name as per Catalyst output)
        //  Must be called within the FIELD_DATA state.
        //  \note As a convenience this can also be called from
        //      (OPENED | DECLARED) states, in which case it invokes
        //      beginFieldData(1) internally.
        using vtk::fileWriter::writeTimeValue;

        //- Write the currently set time as "TimeValue" FieldData
        void writeTimeValue();

#undef  declare_methods
#define declare_methods(FaceType)                                             \
        /*! \brief Reset point/face references to begin a new piece */        \
        void piece(const pointField& points, const List<FaceType>& faces);

        declare_methods(face);
        declare_methods(triFace);
        declare_methods(labelledTri);

        #undef declare_methods
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
