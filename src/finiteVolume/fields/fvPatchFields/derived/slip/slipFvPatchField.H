/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::slipFvPatchField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a slip constraint.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            slip;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef Foam_slipFvPatchField_H
#define Foam_slipFvPatchField_H

#include "basicSymmetryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class slipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class slipFvPatchField
:
    public basicSymmetryFvPatchField<Type>
{
    typedef slipFvPatchField<Type> this_bctype;
    typedef basicSymmetryFvPatchField<Type> parent_bctype;

public:

    //- Runtime type information
    TypeName("slip");


    // Constructors

        //- Construct from patch and internal field
        slipFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF
        )
        :
            parent_bctype(p, iF)
        {}

        //- Construct from patch, internal field and dictionary
        slipFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const dictionary& dict
        )
        :
            parent_bctype(p, iF, dict)
        {}

        //- Construct by mapping onto a new patch
        slipFvPatchField
        (
            const this_bctype& ptf,
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const fvPatchFieldMapper& mapper
        )
        :
            parent_bctype(ptf, p, iF, mapper)
        {}

        //- Construct as copy setting internal field reference
        slipFvPatchField
        (
            const this_bctype& ptf,
            const DimensionedField<Type, volMesh>& iF
        )
        :
            parent_bctype(ptf, iF)
        {}

        //- No copy without an internal field
        [[deprecated("2026-03: soon = delete")]]
        slipFvPatchField(const this_bctype& pfld)
        #ifdef Foam_fvPatchField_copyConstruct
        : this_bctype(pfld, pfld.internalField()) {}
        #else
        = delete;
        #endif

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        // Attributes

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }



    // Member Operators

        // Disable assignment operators
        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
