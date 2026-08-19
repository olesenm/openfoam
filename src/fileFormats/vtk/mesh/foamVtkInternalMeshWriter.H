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
    Foam::vtk::internalMeshWriter

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
    Foam::vtk::internalWriter
    Foam::vtk::vtuWriter

SourceFiles
    foamVtkInternalMeshWriter.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_internalMeshWriter_H
#define Foam_vtk_internalMeshWriter_H

#include "foamVtuWriter.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                   Class vtk::internalMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class internalMeshWriter
:
    public vtk::vtuWriter
{
protected:

    // Protected Data

        //- Reference to the OpenFOAM mesh (or subset)
        const polyMesh& mesh_;

public:

    //- Debug information
    static int debug;

    // Generated Methods

        //- No copy construct
        internalMeshWriter(const internalMeshWriter&) = delete;

        //- No copy assignment
        void operator=(const internalMeshWriter&) = delete;


    // Constructors

        //- Construct from components (default format INLINE_BASE64),
        internalMeshWriter
        (
            const polyMesh& mesh,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        internalMeshWriter
        (
            const polyMesh& mesh,
            const vtk::vtuCells& cells,
            const fileName& file,
            bool parallel = UPstream::parRun()
        );

        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        internalMeshWriter
        (
            const polyMesh& mesh,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = UPstream::parRun()
        );


    //- Destructor
    virtual ~internalMeshWriter() = default;


    // Member Functions

        //- File extension for current format type or given output
        using vtk::vtuWriter::ext;

        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED)
        virtual bool beginFile(std::string title = "");

        //- Write mesh topology.
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
