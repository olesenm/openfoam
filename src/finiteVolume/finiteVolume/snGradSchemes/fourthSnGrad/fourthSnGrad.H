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
    Foam::fv::fourthSnGrad

Group
    grpFvSnGradSchemes

Description
    Fourth-order snGrad scheme with non-orthogonal correction.

    Combines a standard orthogonal difference with a fourth-order correction
    derived from the difference between the corrected surface-normal gradient
    and the face-normal component of the Gauss cell gradient, scaled by a
    factor of 1/6, plus the full non-orthogonal correction.

Usage
    Minimal example by using \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)       fourth;
    }
    \endverbatim

    or as part of a laplacianScheme:
    \verbatim
    laplacianSchemes
    {
        default         Gauss linear fourth;
    }
    \endverbatim

SourceFiles
    fourthSnGrad.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_fourthSnGrad_H
#define Foam_fourthSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class fourthSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fourthSnGrad
:
    public snGradScheme<Type>
{
public:

    //- Runtime type information
    TypeName("fourth");


    // Constructors

        //- Construct from mesh
        fourthSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}

        //- Construct from mesh and data stream
        fourthSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~fourthSnGrad() = default;


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const noexcept { return true; }

        //- Return the explicit correction to the fourthSnGrad
        //- for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fourthSnGrad.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
