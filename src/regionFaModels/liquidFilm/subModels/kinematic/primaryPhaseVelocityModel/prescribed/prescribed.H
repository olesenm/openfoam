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
    Foam::regionModels::areaSurfaceFilmModels::prescribed

Description
    Explicitly prescribed primary-phase velocity, either as a uniform vector
    or mapped from a named \c volVectorField.

Usage
    Minimal example:
    \verbatim
    {
        // Optional entries
        primaryPhaseVelocityModel   prescribed;

        prescribedCoeffs
        {
            value  <vector>;  // (30 0 0);
            // OR
            field   <word>;   // U.gas;
        }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                 | Type   | Reqd   | Deflt
      value       | Uniform prescribed velocity | vector | cndtnl | -
      field       | Name of the volume field    | word   | cndtnl | -
    \endtable

SourceFiles
    prescribed.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_prescribed_H
#define regionFaModels_prescribed_H

#include "primaryPhaseVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class prescribed Declaration
\*---------------------------------------------------------------------------*/

class prescribed
:
    public primaryPhaseVelocityModel
{
    // Private Data

        //- Uniform prescribed velocity
        vector value_;

        //- Name of the volume field to map from
        word fieldName_;


public:

    //- Runtime type information
    TypeName("prescribed");


    // Constructors

        //- Construct from film and dictionary
        prescribed
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~prescribed() = default;


    // Member Functions

        //- Return prescribed tangential primary-phase velocity
        virtual tmp<areaVectorField> Up() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
