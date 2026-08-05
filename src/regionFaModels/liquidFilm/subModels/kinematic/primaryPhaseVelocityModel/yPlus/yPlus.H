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
    Foam::regionModels::areaSurfaceFilmModels::yPlus

Description
    Primary-phase velocity from a rough-wall logarithmic law evaluated at a
    target \f$y^+\f$. The local friction velocity is obtained from the
    primary-region wall-function yPlus calculation, with direction taken
    from the patch-internal tangential primary-phase velocity.

Usage
    Minimal example:
    \verbatim
    {
        // Optional entries
        primaryPhaseVelocityModel   yPlus;

        yPlusCoeffs
        {
            yPlus   <scalar>;
            kappa   <scalar>;
            E       <scalar>;
            Ks      <scalar>;
            Cs      <scalar>;
        }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type   | Reqd | Deflt
      yPlus     | Target yPlus value             | scalar | yes  | -
      kappa     | von Karman constant            | scalar | no   | 0.41
      E         | Wall roughness parameter       | scalar | no   | 9.8
      Ks        | Roughness height [m]           | scalar | no   | 0
      Cs        | Roughness constant             | scalar | no   | 0
    \endtable

Note
  - The model requires a turbulence model on the primary region.
  - The selected primary-region wall function supplies the local friction
    velocity; \c Ks and \c Cs apply to the target log-law.

SourceFiles
    yPlus.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_yPlus_H
#define regionFaModels_yPlus_H

#include "primaryPhaseVelocityModel.H"
#include "wallFunctionCoefficients.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                           Class yPlus Declaration
\*---------------------------------------------------------------------------*/

class yPlus
:
    public primaryPhaseVelocityModel
{
    // Private Data

        //- Wall-function coefficients
        wallFunctionCoefficients wallCoeffs_;

        //- Target yPlus value
        scalar yPlusTarget_;

        //- Roughness height
        scalar Ks_;

        //- Roughness model constant
        scalar Cs_;


public:

    //- Runtime type information
    TypeName("yPlus");


    // Constructors

        //- Construct from film and dictionary
        yPlus
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~yPlus() = default;


    // Member Functions

        //- Return log-law tangential primary-phase velocity at target yPlus
        virtual tmp<areaVectorField> Up() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
