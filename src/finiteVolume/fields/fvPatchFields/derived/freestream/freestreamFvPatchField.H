/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::freestreamFvPatchField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a free-stream condition.  It is a 'mixed'
    condition derived from the \c inletOutlet condition, whereby the mode of
    operation switches between fixed (free stream) value and zero gradient
    based on the sign of the flux.

Usage
    \table
        Property        | Description           | Required    | Default value
        freestreamValue | freestream velocity   | no          |
        freestreamBC    | patchField providing the inlet field | no  |
        phi             | flux field name       | no          | phi
    \endtable

    Example of the boundary condition specification using a fixed value
    inlet condition:
    \verbatim
    <patchName>
    {
        type            freestream;
        freestreamValue uniform (300 0 0);
    }
    \endverbatim

    Example of the boundary condition specification using a separate condition:
    \verbatim
    <patchName>
    {
        type            freestream;
        freestreamBC
        {
            type        timeVaryingMappedFixedValue;
            ...
        }
    }
    \endverbatim

See also
    Foam::mixedFvPatchField
    Foam::inletOutletFvPatchField

SourceFiles
    freestreamFvPatchField.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_freestreamFvPatchField_H
#define Foam_freestreamFvPatchField_H

#include "inletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class freestreamFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class freestreamFvPatchField
:
    public inletOutletFvPatchField<Type>
{
    typedef freestreamFvPatchField<Type> this_bctype;
    typedef inletOutletFvPatchField<Type> parent_bctype;

    // Private Data

        //- BC to supply the freestream value
        tmp<fvPatchField<Type>> freestreamBCPtr_;


public:

    //- Runtime type information
    TypeName("freestream");


    // Constructors

        //- Construct from patch and internal field
        freestreamFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freestreamFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        freestreamFvPatchField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        freestreamFvPatchField
        (
            const this_bctype&,
            const DimensionedField<Type, volMesh>&
        );

        //- No copy without an internal field
        [[deprecated("2026-03: soon = delete")]]
        freestreamFvPatchField(const this_bctype& pfld)
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

        //- Same as refValue()
        const Field<Type>& freestreamValue() const { return this->refValue(); }

        //- Same as refValue()
        Field<Type>& freestreamValue() { return this->refValue(); }


    // Mapping Functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper& mapper);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchField<Type>& ptf,
            const labelList& addr
        );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "freestreamFvPatchField.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
