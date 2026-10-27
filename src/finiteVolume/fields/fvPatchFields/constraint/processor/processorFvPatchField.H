/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    Foam::processorFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition enables processor communication across patches.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            processor;
    }
    \endverbatim

SourceFiles
    processorFvPatchField.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorFvPatchField_H
#define Foam_processorFvPatchField_H

#include "coupledFvPatchField.H"
#include "processorLduInterfaceField.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class processorFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorFvPatchField
:
    public processorLduInterfaceField,
    public coupledFvPatchField<Type>
{
    typedef processorFvPatchField<Type> this_bctype;
    typedef coupledFvPatchField<Type> parent_bctype;

    // Private Data

        //- Local reference cast into the processor patch
        const processorFvPatch& procPatch_;


        // Sending and receiving

            //- Current (non-blocking) send request
            mutable int sendRequest_;

            //- Current (non-blocking) recv request
            mutable int recvRequest_;

            //- Send buffer.
            mutable Field<Type> sendBuf_;

            //- Receive buffer.
            mutable Field<Type> recvBuf_;

            //- Scalar send buffer
            mutable solveScalarField scalarSendBuf_;

            //- Scalar recv buffer
            mutable solveScalarField scalarRecvBuf_;


    // Private Member Functions

        //- Receive and send requests have both completed
        virtual bool all_ready() const;


public:

    //- Runtime type information
    TypeName(processorFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        processorFvPatchField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        processorFvPatchField
        (
            const this_bctype&,
            const DimensionedField<Type, volMesh>&
        );

        //- No copy without an internal field
        [[deprecated("2026-03: soon = delete")]]
        processorFvPatchField(const this_bctype& pfld)
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


    //- Destructor
    ~processorFvPatchField() = default;


    // Member Functions

    // Coupling

        //- The patch field is coupled if running in parallel
        virtual bool coupled() const { return UPstream::parRun(); }

        //- Are all (receive) data available?
        virtual bool ready() const;

        //- Return neighbour field
        virtual tmp<Field<Type>> patchNeighbourField() const;

        //- Retrieve neighbour field
        virtual void patchNeighbourField(UList<Type>& pnf) const;


    // Evaluation

        //- Initialise the evaluation of the patch field
        virtual void initEvaluate(const Pstream::commsTypes commsType);

        //- Evaluate the patch field
        virtual void evaluate(const Pstream::commsTypes commsType);

        //- Initialise the evaluation of the patch field after a local
        //  operation. Dummy since operating on a copy
        virtual void initEvaluateLocal
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::buffered
        )
        {}

        //- Evaluate the patch field after a local operation (e.g. *=).
        //  Dummy since operating on a copy
        virtual void evaluateLocal
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::buffered
        )
        {}

        //- Return patch-normal gradient
        virtual tmp<Field<Type>> snGrad
        (
            const scalarField& deltaCoeffs
        ) const;

        //- Retrieve patch-normal gradient
        virtual void snGrad
        (
            const scalarField& deltaCoeffs,
            UList<Type>& result
        ) const;


        // Coupled interface functionality

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                Field<Type>& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;


    // Processor coupled interface functions

        //- Return communicator used for communication
        virtual label comm() const
        {
            return procPatch_.comm();
        }

        //- Return processor number
        virtual int myProcNo() const
        {
            return procPatch_.myProcNo();
        }

        //- Return neighbour processor number
        virtual int neighbProcNo() const
        {
            return procPatch_.neighbProcNo();
        }

        //- Does the patch field perform the transformation
        virtual bool doTransform() const
        {
            return
            (
                is_rotational_vectorspace_v<Type>
             && !procPatch_.parallel()
            );
        }

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const
        {
            return procPatch_.forwardT();
        }

        //- Return rank of component for transform
        virtual int rank() const
        {
            return pTraits<Type>::rank;
        }


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "processorFvPatchField.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
