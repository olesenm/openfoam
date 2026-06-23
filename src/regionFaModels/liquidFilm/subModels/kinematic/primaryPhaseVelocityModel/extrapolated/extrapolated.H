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
    Foam::regionModels::areaSurfaceFilmModels::extrapolated

Description
    Primary-phase velocity from linear near-wall extrapolation:

    \f[
        \vec{U}_p = h \, \nabla_n \vec{U}
    \f]

    where
    \vartable
        \vec{U}   | Primary-phase velocity
        h         | Film height
        \vec{U}_p | Tangential primary-phase velocity at the film surface
    \endvartable

Usage
    Minimal example:
    \verbatim
    {
        // Optional entries
        primaryPhaseVelocityModel  extrapolated;

        // Inherited entries
        type                       velocityFilmShell;
        ...
    }
    \endverbatim

SourceFiles
    extrapolated.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_extrapolated_H
#define regionFaModels_extrapolated_H

#include "primaryPhaseVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class extrapolated Declaration
\*---------------------------------------------------------------------------*/

class extrapolated
:
    public primaryPhaseVelocityModel
{
public:

    //- Runtime type information
    TypeName("extrapolated");


    // Constructors

        //- Construct from film and dictionary
        extrapolated
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~extrapolated() = default;


    // Member Functions

        //- Return extrapolated tangential primary-phase velocity
        virtual tmp<areaVectorField> Up() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
