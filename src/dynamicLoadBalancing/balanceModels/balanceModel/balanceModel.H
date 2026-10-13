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

Namespace
    Foam::balanceModels

Description
    A namespace for various balance model implementations for the dynamic
    load balancing.

Class
    Foam::balanceModel

Description
    A base class for the balance models.

SourceFiles
    balanceModel.C
    balanceModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_balanceModel_H
#define Foam_balanceModel_H

#include "volFields.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class balanceModel Declaration
\*---------------------------------------------------------------------------*/

class balanceModel
{
    // Private Data

        //- Reference to the mesh
        fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("balanceModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            balanceModel,
            dictionary,
            (
                fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected performance model
        static autoPtr<balanceModel> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from mesh and dictionary
        balanceModel
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        balanceModel(const balanceModel&) = delete;

        //- No copy assignment
        void operator=(const balanceModel&) = delete;


    //- Destructor
    virtual ~balanceModel();


    // Member Functions

    // Access

        //- Return reference to the mesh
        fvMesh& mesh() const noexcept { return mesh_; }


    // Evaluation

        //- Balance the load by redistributing the mesh and fields
        virtual bool balance(const scalarField& weights) = 0;


    // I-O

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict) { return true; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
