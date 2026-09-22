/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
    Copyright (C) 2026 Patrick Wilhelm, Data Lab Hell
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
    Foam::wallHeatFluxModels::wall

Description
    This model computes the wall-heat flux at selected wall patches.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wallHeatFlux1
    {
        // Mandatory entries
        type        wallHeatFlux;
        libs        (fieldFunctionObjects);

        model       wall;

        // Optional entries
        patches     (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");
        qr          <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Reqd | Deflt
      type         | Type name: wallHeatFlux               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patches    | Names of operand patches   | wordList | no | all wall patches
      qr         | Name of radiative heat flux field       | word | no | qr
    \endtable

    The inherited entries are elaborated in:
     - \link fvMeshFunctionObject.H \endlink
     - \link writeFile.H \endlink

Note
  - The model \c wall corresponds to \c wallHeatFlux in OpenFOAM v2506 and
    earlier versions.

SourceFiles
    wallHeatFlux_wall.cxx

\*---------------------------------------------------------------------------*/

#ifndef wallHeatFluxModels_wall_H
#define wallHeatFluxModels_wall_H

#include "wallHeatFluxModel.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatFluxModels
{

/*---------------------------------------------------------------------------*\
                   Class wall Declaration
\*---------------------------------------------------------------------------*/

class wall
:
    public wallHeatFluxModel
{
    // Private Data

        //- Name of the radiative heat-flux field
        word qrName_;

        //- Wall patches to process (optionally filtered by name)
        labelList patchIDs_;

        //- Set once the anisotropic-solid warning has been given
        bool warnedAnisotropic_;


    // Private Member Functions

        //- Write file-header information
        void writeFileHeader(Ostream& os);

        //- Calculate the wall heat-flux
        void calcHeatFlux
        (
            const volScalarField& gamma,
            const volScalarField& vf,
            volScalarField& wallHeatFlux
        );


public:

    //- Runtime type information
    TypeName("wall");


    // Constructors

        //- Construct from components
        wall
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& name,
            const word objName,
            functionObjects::stateFunctionObject& state
        );


    //- Destructor
    virtual ~wall() = default;


    // Member Functions

    // Access

        //- Return const reference to name of radiative heat-flux field
        const word& qrName() const noexcept { return qrName_; }


    // Evaluation

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Calculate the heat-flux data
        virtual bool execute();

        //- Write the heat-flux data
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallHeatFluxModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
