/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::functionObjects::processorField

Group
    grpFieldFunctionObjects

Description
    Writes a scalar field whose value is the local processor ID.  The output
    field name is \c processorID.

    Operands:
    \table
      Operand        | Type           | Location
      input          | volScalarField | <time>/inputField
      output file    | -              | -
      output field   | volScalarField | <time>/outputField
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    processorFieldFO
    {
        // Mandatory entries
        type        processorField;
        libs        (fieldFunctionObjects);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: processorField          | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink

SourceFiles
    processorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_processorField_H
#define Foam_functionObjects_processorField_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class processorField Declaration
\*---------------------------------------------------------------------------*/

class processorField
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Processor ID field (not stored on mesh registry to avoid
        //  conflicts with fvMesh::mapFields and fvMeshDistribute)
        autoPtr<volScalarField> procFieldPtr_;


public:

    //- Runtime type information
    TypeName("processorField");


    // Constructors

        //- Construct from name, Time and dictionary
        processorField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        processorField(const processorField&) = delete;

        //- No copy assignment
        void operator=(const processorField&) = delete;


    //- Destructor
    virtual ~processorField() = default;


    // Member Functions

        //- Create the processorID field if not present
        void createProcessorIDField();

        //- Read the function-object dictionary
        virtual bool read(const dictionary& dict);

        //- Execute the function-object operations
        virtual bool execute();

        //- Write the function-object results
        virtual bool write();

        //- Update for changes of mesh
        //  The base implementation is a no-op.
        virtual void updateMesh(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
