/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
    Foam::vtk::patchWriter

Description
    Write OpenFOAM patches and patch fields in VTP or legacy vtk format.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

See Also
    Foam::vtk::patchMeshWriter

SourceFiles
    foamVtkPatchWriter.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_patchWriter_H
#define Foam_vtk_patchWriter_H

#include "foamVtkPatchMeshWriter.H"
#include "volFields.H"
#include "pointFields.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                      Class vtk::patchWriter Declaration
\*---------------------------------------------------------------------------*/

class patchWriter
:
    public vtk::patchMeshWriter
{
    // Private Member Data

        //- Use patch internal field values instead of patch value
        bool useNearCellValue_ = false;

public:

    // Generated Methods

        //- No copy construct
        patchWriter(const patchWriter&) = delete;

        //- No copy assignment
        void operator=(const patchWriter&) = delete;


    // Constructors

        //- Inherit all vtk::patchMeshWriter constructors
        using vtk::patchMeshWriter::patchMeshWriter;


    //- Destructor
    virtual ~patchWriter() = default;


    // Member Functions

    // Access / Edit

        //- Use patch internal field values?
        bool useCellValue() const noexcept
        {
            return useNearCellValue_;
        }

        //- Select use of patch internal field values
        //  \returns the previous state
        bool useCellValue(bool on) noexcept
        {
            bool old(useNearCellValue_);
            useNearCellValue_ = on;
            return old;
        }


    // Write

        //- Write point field
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, pointMesh>& field
        );

        //- Write volume field
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, volMesh>& field
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PrimitivePatchInterpolation<primitivePatch>& pInterp
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PrimitivePatchInterpolation<primitivePatch>* pInterp
        );


    // Housekeeping

        //- Construct with \p useNearCellValue setting
        //  \deprecated
        FOAM_DEPRECATED_FOR(2025-11, "use useCellValue() after construction")
        patchWriter
        (
            const polyMesh& mesh,
            const labelUList& patchIDs,
            const vtk::outputOptions opts,
            //! Use patch internal instead of patch face values
            const bool useNearCellValue,
            const fileName& file,
            bool parallel = UPstream::parRun()
        )
        :
            patchWriter(mesh, patchIDs, opts, file, parallel)
        {
            useCellValue(useNearCellValue);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkPatchWriter.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
