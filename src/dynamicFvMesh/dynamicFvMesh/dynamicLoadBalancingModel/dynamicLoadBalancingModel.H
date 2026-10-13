/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::dynamicLoadBalancingModel

Description
    Runtime-selectable dynamic-load-balancing model executed from
    dynamicFvMesh::controlledUpdate().

    The model is selected from \c constant/dynamicLoadBalanceDict once per
    mesh if the dictionary is present. If absent, dynamic load balancing is
    disabled. Later edits to the dictionary do not reconfigure the model.

SourceFiles
    dynamicLoadBalancingModel.C
    dynamicLoadBalancingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dynamicLoadBalancingModel_H
#define Foam_dynamicLoadBalancingModel_H

#include "autoPtr.H"
#include "className.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class dynamicLoadBalancingModel Declaration
\*---------------------------------------------------------------------------*/

class dynamicLoadBalancingModel
{
    // Private Data

        //- Finite-volume mesh
        fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("dynamicLoadBalancingModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicLoadBalancingModel,
            dictionary,
            (
                fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Return the selected dynamic-load-balancing model, or nullptr
        //  if constant/dynamicLoadBalanceDict is absent or inactive.
        static autoPtr<dynamicLoadBalancingModel> New(fvMesh& mesh);


    // Constructors

        //- Construct from mesh and dictionary
        dynamicLoadBalancingModel
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        dynamicLoadBalancingModel(const dynamicLoadBalancingModel&) = delete;

        //- No copy assignment
        void operator=(const dynamicLoadBalancingModel&) = delete;


    //- Destructor
    virtual ~dynamicLoadBalancingModel();


    // Member Functions

    // Access

        //- Return finite-volume mesh
        fvMesh& mesh() const noexcept { return mesh_; }


    // Evaluation

        //- Execute the model and return mesh-changed status
        virtual bool execute() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
