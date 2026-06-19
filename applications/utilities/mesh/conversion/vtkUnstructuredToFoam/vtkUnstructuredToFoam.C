/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

Application
    vtkUnstructuredToFoam

Group
    grpMeshConversionUtilities

Description
    Convert legacy VTK file (ascii) containing an unstructured grid
    to an OpenFOAM mesh without boundary information.

Usage
    \b vtkUnstructuredToFoam \<XXX.vtk\>

    Options:
      - \par -no-fields
        Do not attempt to recreate volFields

Note
  - The .vtk format does not contain any boundary information. Boundary fields
    are written with \c zeroGradient type since VTK legacy format cannot
    preserve patch names or boundary condition types.
  - Known field names (T, U, k) are assigned physical dimensions; all other
    imported fields are dimensionless.

    Not extensively tested.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "vtkUnstructuredReader.H"

#include "columnFvMesh.H"
#include "dimensionSets.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "vectorIOField.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const dimensionSet& fieldDimensions(const word& name)
{
    static const HashTable<dimensionSet> dims
    ({
        { "T", dimTemperature },
        { "U", dimVelocity },
        { "k", sqr(dimVelocity) },
    });

    return dims.found(name) ? dims[name] : dimless;
}


template<class Type>
void writeImportedVolField
(
    fvMesh& mesh,
    const word& name,
    const Field<Type>& internal
)
{
    GeometricField<Type, fvPatchField, volMesh> fld
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        fieldDimensions(name),
        internal,
        fvPatchFieldBase::zeroGradientType()
    );

    fld.correctBoundaryConditions();
    fld.write();
}


template<class Type>
void constructVolFieldsFromIO
(
    fvMesh& mesh,
    const vtkUnstructuredReader& reader
)
{
    const labelUList& cellMap = reader.cellMap();

    for (const auto& field : reader.cellData().csorted<IOField<Type>>())
    {
        Info<< "Constructing volField " << field.name() << endl;

        writeImportedVolField
        (
            mesh,
            field.name(),
            Field<Type>(UIndirectList<Type>(field, cellMap))
        );
    }
}


template<class Type>
void constructVolFields(fvMesh& mesh, const vtkUnstructuredReader& reader)
{
    constructVolFieldsFromIO<Type>(mesh, reader);
}


template<>
void constructVolFields<scalar>
(
    fvMesh& mesh,
    const vtkUnstructuredReader& reader
)
{
    constructVolFieldsFromIO<scalar>(mesh, reader);

    const labelUList& cellMap = reader.cellMap();

    for (const auto& field : reader.cellData().csorted<labelIOField>())
    {
        Info<< "Constructing volScalarField " << field.name()
            << " (from label VTK data)" << endl;

        scalarField internal(cellMap.size());

        forAll(cellMap, celli)
        {
            internal[celli] = scalar(field[cellMap[celli]]);
        }

        writeImportedVolField(mesh, field.name(), internal);
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert legacy VTK file (ascii) containing an unstructured grid"
        " to an OpenFOAM mesh without boundary information"
    );

    argList::noParallel();
    argList::addOptionCompat("no-fields", {"noFields", 2106});
    argList::addArgument("vtk-file", "The input legacy ascii vtk file");

    #include "setRootCase.H"
    #include "createTime.H"

    const bool doFields = !args.found("no-fields");

    IFstream mshStream(args.get<fileName>(1));

    vtkUnstructuredReader reader(runTime, mshStream);

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        std::move(reader.points()),
        reader.cells(),
        faceListList(),
        wordList(),
        wordList(),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );

    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< "Writing mesh ..." << endl;

    mesh.removeFiles();
    mesh.write();


    if (doFields)
    {
        // Re-read mesh as fvMesh so we can have fields
        Info<< "Re-reading mesh ..." << endl;
        #include "createMesh.H"

        constructVolFields<scalar>(mesh, reader);
        constructVolFields<vector>(mesh, reader);
        constructVolFields<sphericalTensor>(mesh, reader);
        constructVolFields<symmTensor>(mesh, reader);
        constructVolFields<tensor>(mesh, reader);
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
