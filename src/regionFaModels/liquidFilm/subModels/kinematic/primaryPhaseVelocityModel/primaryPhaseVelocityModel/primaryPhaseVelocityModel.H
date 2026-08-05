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
    Foam::regionModels::areaSurfaceFilmModels::primaryPhaseVelocityModel

Description
    Abstract base class for primary-phase velocity models used by
    \c liquidFilmBase::Up().

    The model is selected via the \c primaryPhaseVelocityModel keyword in the
    \c velocityFilmShell boundary-condition dictionary with the default
    value \c extrapolated for backward compatibility.

SourceFiles
    primaryPhaseVelocityModel.C
    primaryPhaseVelocityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_primaryPhaseVelocityModel_H
#define regionFaModels_primaryPhaseVelocityModel_H

#include "liquidFilmBase.H"
#include "areaFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                  Class primaryPhaseVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class primaryPhaseVelocityModel
{
protected:

    // Protected Data

        //- Const reference to the liquid film model
        const liquidFilmBase& film_;

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;


    // Protected Member Functions

        //- Return model coefficients sub-dictionary
        static const dictionary& coeffsDict
        (
            const word& modelType,
            const dictionary& dict
        );

        //- Return a new primary-phase velocity area field
        tmp<areaVectorField> newUpField() const;

        //- Remove the normal component of a surface velocity field
        void projectTangent(vectorField& Up) const;


public:

    //- Runtime type information
    TypeName("primaryPhaseVelocityModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            primaryPhaseVelocityModel,
            dictionary,
            (
                const liquidFilmBase& film,
                const dictionary& dict
            ),
            (film, dict)
        );


    // Generated Methods

        //- No copy construct
        primaryPhaseVelocityModel(const primaryPhaseVelocityModel&) = delete;

        //- No copy assignment
        void operator=(const primaryPhaseVelocityModel&) = delete;


    // Constructors

        //- Construct from model type, film, and dictionary
        primaryPhaseVelocityModel
        (
            const word& modelType,
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected primary-phase velocity model
        static autoPtr<primaryPhaseVelocityModel> New
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~primaryPhaseVelocityModel() = default;


    // Member Functions

    // Access

        //- Return const access to the film model
        const liquidFilmBase& film() const noexcept { return film_; }

        //- Return const access to the finite-area mesh
        const faMesh& mesh() const noexcept { return mesh_; }


    // Evaluation

        //- Return the tangential primary-phase velocity at the film surface
        virtual tmp<areaVectorField> Up() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
