/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "cyclicAMIPolyPatch.H"
#include "mapDistributeBase.H"
#include "AMIInterpolation.H"
#include "fvMatrix.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    parent_bctype(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    patchNeighbourFieldPtr_(nullptr)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMILduInterfaceField(),
    parent_bctype(p, iF, dict, IOobjectOption::NO_READ),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p, dict)),
    patchNeighbourFieldPtr_(nullptr)
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (cacheNeighbourField())
    {
        // Handle neighbour value first, before any evaluate()
        const auto* hasNeighbValue =
            dict.findEntry("neighbourValue", keyType::LITERAL);

        if (hasNeighbValue)
        {
            patchNeighbourFieldPtr_.reset
            (
                new Field<Type>(*hasNeighbValue, p.size())
            );
        }
    }

    // Use 'value' supplied, or evaluate (if coupled) or set to internal field
    if (!this->readValueEntry(dict))
    {
        if (this->coupled())
        {
            // Tricky: avoid call to evaluate without call to initEvaluate.
            // For now just disable the localConsistency to make it use the
            // old logic (ultimately calls the fully self contained
            // patchNeighbourField)

            const auto oldConsistency = FieldBase::localBoundaryConsistency(0);

            this->evaluate(UPstream::commsTypes::nonBlocking);

            FieldBase::localBoundaryConsistency(oldConsistency);
        }
        else
        {
            this->extrapolateInternal();  // Zero-gradient patch values
        }
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMILduInterfaceField(),
    parent_bctype(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    patchNeighbourFieldPtr_(nullptr)
{
    //if (ptf.patchNeighbourFieldPtr_ && cacheNeighbourField())
    //{
    //    patchNeighbourFieldPtr_.reset
    //    (
    //        new Field<Type>(ptf.patchNeighbourFieldPtr_(), mapper)
    //    );
    //}

    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << cyclicAMIPatch_.name()
            << abort(FatalError);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const this_bctype& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    parent_bctype(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    patchNeighbourFieldPtr_(nullptr)
{
    if (debug && !ptf.all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on patch " << cyclicAMIPatch_.name()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::all_ready() const
{
    int done = 0;

    if
    (
        UPstream::finishedRequests
        (
            recvRequests_.start(),
            recvRequests_.size()
        )
    &&  UPstream::finishedRequests
        (
            recvRequests1_.start(),
            recvRequests1_.size()
        )
    )
    {
        recvRequests_.clear();
        recvRequests1_.clear();
        ++done;
    }

    if
    (
        UPstream::finishedRequests
        (
            sendRequests_.start(),
            sendRequests_.size()
        )
    &&  UPstream::finishedRequests
        (
            sendRequests1_.start(),
            sendRequests1_.size()
        )
    )
    {
        sendRequests_.clear();
        sendRequests1_.clear();
        ++done;
    }

    return (done == 2);
}


template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::ready() const
{
    if
    (
        UPstream::finishedRequests
        (
            recvRequests_.start(),
            recvRequests_.size()
        )
     && UPstream::finishedRequests
        (
            recvRequests1_.start(),
            recvRequests1_.size()
        )
    )
    {
        recvRequests_.clear();
        recvRequests1_.clear();

        if
        (
            UPstream::finishedRequests
            (
                sendRequests_.start(),
                sendRequests_.size()
            )
         && UPstream::finishedRequests
            (
                sendRequests1_.start(),
                sendRequests1_.size()
            )
        )
        {
            sendRequests_.clear();
            sendRequests1_.clear();
        }

        return true;
    }

    return false;
}



template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    this->parent_bctype::autoMap(mapper);
    patchNeighbourFieldPtr_.reset(nullptr);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    this->parent_bctype::rmap(ptf, addr);
    patchNeighbourFieldPtr_.reset(nullptr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::getNeighbourField
(
    const UList<Type>& internalData
) const
{
    // By pass polyPatch to get nbrId. Instead use cyclicAMIFvPatch virtual
    // neighbPatch()
    const auto& neighbPatch = cyclicAMIPatch_.neighbPatch();
    const labelUList& nbrFaceCells = neighbPatch.faceCells();

    Field<Type> pnf(internalData, nbrFaceCells);
    Field<Type> defaultValues;

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        defaultValues = Field<Type>(internalData, cyclicAMIPatch_.faceCells());
    }

    tmp<Field<Type>> tpnf = cyclicAMIPatch_.interpolate(pnf, defaultValues);

    if (doTransform())
    {
        transform(tpnf.ref(), forwardT(), tpnf());
    }

    return tpnf;
}


template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::cacheNeighbourField() const
{
    // const auto& AMI = this->ownerAMI();

    // if (AMI.cacheActive())
    // {
    //     return false;
    // }
    // else
    {
        return (FieldBase::localBoundaryConsistency() != 0);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::getPatchNeighbourField
(
    bool checkCommunicator
) const
{
    const auto& AMI = this->ownerAMI();

    if
    (
        AMI.distributed() && cacheNeighbourField()
     && (!checkCommunicator || AMI.comm() != -1)
    )
    {
        if (!this->ready())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicAMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        const auto& fvp = this->patch();
        const auto& mesh = fvp.boundaryMesh().mesh();

        if
        (
            patchNeighbourFieldPtr_
        && !mesh.upToDatePoints(this->internalField())
        )
        {
            //DebugPout
            //    << "cyclicAMIFvPatchField::patchNeighbourField() :"
            //    << " field:" << this->internalField().name()
            //    << " patch:" << fvp.name()
            //    << " CLEARING patchNeighbourField"
            //    << endl;
            patchNeighbourFieldPtr_.reset(nullptr);
        }

        // Initialise if not done in construct-from-dictionary
        if (!patchNeighbourFieldPtr_)
        {
            //DebugPout
            //    << "cyclicAMIFvPatchField::patchNeighbourField() :"
            //    << " field:" << this->internalField().name()
            //    << " patch:" << fvp.name()
            //    << " caching patchNeighbourField"
            //    << endl;

            // Do interpolation and store result
            patchNeighbourFieldPtr_.reset
            (
                getNeighbourField(this->primitiveField()).ptr()
            );
        }
        else
        {
            // Have cached value. Check
            //if (debug)
            //{
            //    tmp<Field<Type>> tpnf
            //    (
            //        getNeighbourField(this->primitiveField())
            //    );
            //    if (tpnf() != patchNeighbourFieldPtr_())
            //    {
            //        FatalErrorInFunction
            //            << "On field " << this->internalField().name()
            //            << " patch " << fvp.name() << endl
            //            << "Cached patchNeighbourField    :"
            //            << flatOutput(patchNeighbourFieldPtr_()) << endl
            //            << "Calculated patchNeighbourField:"
            //            << flatOutput(tpnf()) << exit(FatalError);
            //    }
            //}
        }

        return patchNeighbourFieldPtr_();
    }
    else
    {
        // Do interpolation
        return getNeighbourField(this->primitiveField());
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    // checkCommunicator = true
    return this->getPatchNeighbourField(true);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField
(
    UList<Type>& pnf
) const
{
    // checkCommunicator = false
    auto tpnf = this->getPatchNeighbourField(false);
    pnf.deepCopy(tpnf());
}


template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::neighbourPatchField() const
{
    const auto& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const this_bctype>
    (
        fld.boundaryField()[cyclicAMIPatch_.neighbPatchID()]
    );
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const auto& AMI = this->ownerAMI();

    if (AMI.distributed() && cacheNeighbourField() && AMI.comm() != -1)
    {
        //DebugPout
        //    << "*** cyclicAMIFvPatchField::initEvaluate() :"
        //    << " field:" << this->internalField().name()
        //    << " patch:" << this->patch().name()
        //    << " sending patchNeighbourField"
        //    << endl;

        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            // Invalidate old field - or flag as fatal?
            patchNeighbourFieldPtr_.reset(nullptr);
            return;
        }

        // Start sending

        // Bypass polyPatch to get nbrId.
        // - use cyclicACMIFvPatch::neighbPatch() virtual instead
        const cyclicAMIFvPatch& neighbPatch = cyclicAMIPatch_.neighbPatch();
        const labelUList& nbrFaceCells = neighbPatch.faceCells();
        const Field<Type> pnf(this->primitiveField(), nbrFaceCells);

        const cyclicAMIPolyPatch& cpp = cyclicAMIPatch_.cyclicAMIPatch();

        // Assert that all receives are known to have finished
        if (!recvRequests_.empty() || !recvRequests1_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicAMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        // Assume that sends are also OK
        sendRequests_.clear();
        sendRequests1_.clear();

        cpp.initInterpolate
        (
            pnf,
            sendRequests_,
            recvRequests_,
            sendBufs_,
            recvBufs_,

            sendRequests1_,
            recvRequests1_,
            sendBufs1_,
            recvBufs1_
        );
    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const auto& AMI = this->ownerAMI();

    if (AMI.distributed() && cacheNeighbourField() && AMI.comm() != -1)
    {
        // Calculate patchNeighbourField
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        patchNeighbourFieldPtr_.reset(nullptr);

        const cyclicAMIPolyPatch& cpp = cyclicAMIPatch_.cyclicAMIPatch();

        Field<Type> defaultValues;
        if (AMI.applyLowWeightCorrection())
        {
            defaultValues = this->patchInternalField();
        }

        //DebugPout
        //    << "*** cyclicAMIFvPatchField::evaluate() :"
        //    << " field:" << this->internalField().name()
        //    << " patch:" << this->patch().name()
        //    << " receiving&caching patchNeighbourField"
        //    << endl;

        patchNeighbourFieldPtr_.reset
        (
            cpp.interpolate
            (
                Field<Type>::null(),    // Not used for distributed
                recvRequests_,
                recvBufs_,
                recvRequests1_,
                recvBufs1_,
                defaultValues
            ).ptr()
        );

        // Receive requests all handled by last function call
        recvRequests_.clear();
        recvRequests1_.clear();

        if (doTransform())
        {
            // In-place transform
            auto& pnf = *patchNeighbourFieldPtr_;
            transform(pnf, forwardT(), pnf);
        }
    }

    // Use patchNeighbourField() and patchInternalField() to obtain face value
    coupledFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const auto& AMI = this->ownerAMI();

    if (AMI.distributed() && AMI.comm() != -1)
    {
        // Start sending
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicAMIPatch_.neighbPatchID());

        solveScalarField pnf(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf, cmpt);

        const cyclicAMIPolyPatch& cpp = cyclicAMIPatch_.cyclicAMIPatch();

        // Assert that all receives are known to have finished
        if (!recvRequests_.empty() || !recvRequests1_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicAMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        // Assume that sends are also OK
        sendRequests_.clear();
        sendRequests1_.clear();

        cpp.initInterpolate
        (
            pnf,
            sendRequests_,
            recvRequests_,
            scalarSendBufs_,
            scalarRecvBufs_,

            sendRequests1_,
            recvRequests1_,
            scalarSendBufs1_,
            scalarRecvBufs1_
        );
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    //DebugPout<< "cyclicAMIFvPatchField::updateInterfaceMatrix() :"
    //    << " field:" << this->internalField().name()
    //    << " patch:" << this->patch().name()
    //    << endl;

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    const auto& AMI = this->ownerAMI();

    solveScalarField pnf;

    if (AMI.distributed() && AMI.comm() != -1)
    {
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        solveScalarField defaultValues;
        if (AMI.applyLowWeightCorrection())
        {
            defaultValues = solveScalarField(psiInternal, faceCells);
        }

        const cyclicAMIPolyPatch& cpp = cyclicAMIPatch_.cyclicAMIPatch();

        pnf =
            cpp.interpolate
            (
                solveScalarField::null(),   // Not used for distributed
                recvRequests_,
                scalarRecvBufs_,
                recvRequests1_,
                scalarRecvBufs1_,
                defaultValues
            );

        // Receive requests all handled by last function call
        recvRequests_.clear();
        recvRequests1_.clear();
    }
    else
    {
        solveScalarField defaultValues;
        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            defaultValues = solveScalarField(psiInternal, faceCells);
        }

        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicAMIPatch_.neighbPatchID());

        pnf = solveScalarField(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf, cmpt);

        pnf = cyclicAMIPatch_.interpolate(pnf, defaultValues);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);

    this->updatedMatrix(true);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    const auto& AMI = this->ownerAMI();

    if (AMI.distributed() && AMI.comm() != -1)
    {
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicAMIPatch_.neighbPatchID());

        Field<Type> pnf(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf);

        const cyclicAMIPolyPatch& cpp = cyclicAMIPatch_.cyclicAMIPatch();

        // Assert that all receives are known to have finished
        if (!recvRequests_.empty() || !recvRequests1_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicAMIPatch_.name()
                << " field " << this->internalField().name()
                << abort(FatalError);
        }

        // Assume that sends are also OK
        sendRequests_.clear();
        sendRequests1_.clear();

        cpp.initInterpolate
        (
            pnf,
            sendRequests_,
            recvRequests_,
            sendBufs_,
            recvBufs_,

            sendRequests1_,
            recvRequests1_,
            sendBufs1_,
            recvBufs1_
        );
    }

    this->updatedMatrix(false);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    //DebugPout<< "cyclicAMIFvPatchField::updateInterfaceMatrix() :"
    //    << " field:" << this->internalField().name()
    //    << " patch:" << this->patch().name()
    //    << endl;

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    const auto& AMI = this->ownerAMI();

    Field<Type> pnf;

    if (AMI.distributed() && AMI.comm() != -1)
    {
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        const cyclicAMIPolyPatch& cpp = cyclicAMIPatch_.cyclicAMIPatch();

        Field<Type> defaultValues;
        if (AMI.applyLowWeightCorrection())
        {
            defaultValues = Field<Type>(psiInternal, faceCells);
        }

        pnf =
            cpp.interpolate
            (
                Field<Type>::null(),  // Not used for distributed
                recvRequests_,
                recvBufs_,
                recvRequests1_,
                recvBufs1_,
                defaultValues
            );

        // Receive requests all handled by last function call
        recvRequests_.clear();
        recvRequests1_.clear();
    }
    else
    {
        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicAMIPatch_.neighbPatchID());

        pnf = Field<Type>(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf);

        Field<Type> defaultValues;
        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            defaultValues = Field<Type>(psiInternal, faceCells);
        }

        pnf = cyclicAMIPatch_.interpolate(pnf, defaultValues);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);

    this->updatedMatrix(true);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const label mat,
    const direction cmpt
)
{
    if (this->cyclicAMIPatch().owner())
    {
        const label index = this->patch().index();

        // Save old matrix boundary coefficients
        const Field<Type> oldIntCoeffs(matrix.internalCoeffs());
        const Field<Type> oldBndCoeffs(matrix.boundaryCoeffs());
        const labelList oldCoeffsOffset(matrix.coeffsOffset());

        const label globalPatchID =
            matrix.lduMeshAssembly().patchLocalToGlobalMap()[mat][index];

        const Field<scalar> intCoeffsCmpt
        (
            matrix.internalCoeffs(globalPatchID).component(cmpt)
        );

        const Field<scalar> boundCoeffsCmpt
        (
            matrix.boundaryCoeffs(globalPatchID).component(cmpt)
        );

        tmp<Field<scalar>> tintCoeffs(coeffs(matrix, intCoeffsCmpt, mat));
        tmp<Field<scalar>> tbndCoeffs(coeffs(matrix, boundCoeffsCmpt, mat));
        const Field<scalar>& intCoeffs = tintCoeffs.ref();
        const Field<scalar>& bndCoeffs = tbndCoeffs.ref();

        const labelUList& u = matrix.lduAddr().upperAddr();
        const labelUList& l = matrix.lduAddr().lowerAddr();

        label subFaceI = 0;

        const labelList& faceMap =
            matrix.lduMeshAssembly().faceBoundMap()[mat][index];

        forAll (faceMap, j)
        {
            label globalFaceI = faceMap[j];

            const scalar boundCorr = -bndCoeffs[subFaceI];
            const scalar intCorr = -intCoeffs[subFaceI];

            matrix.upper()[globalFaceI] += boundCorr;
            matrix.diag()[u[globalFaceI]] -= intCorr;
            matrix.diag()[l[globalFaceI]] -= boundCorr;

            if (matrix.asymmetric())
            {
                matrix.lower()[globalFaceI] += intCorr;
            }
            subFaceI++;
        }

        // Set internalCoeffs and boundaryCoeffs in the assembly matrix
        // on cyclicAMI patches to be used in the individual matrix by
        // matrix.flux()
        if (matrix.psi(mat).mesh().fluxRequired(this->internalField().name()))
        {
            // Copy the internal and boundary coefficients to the assembly
            // matrix. Note that the sizes of intCoeffs and bndCoeffs have
            // changed!

            const label nbrPathID =
                cyclicAMIPatch_.cyclicAMIPatch().neighbPatchID();

            const label nbrGlobalPatchID =
                matrix.lduMeshAssembly().patchLocalToGlobalMap()
                    [mat][nbrPathID];


            // Work offsets into sizes, correct, and update
            labelList newOffset(matrix.coeffsOffset().size());
            {
                labelList newSizes(matrix.nPatches());
                forAll(newSizes, patchi)
                {
                    newSizes[patchi] =
                       matrix.coeffsOffset()[patchi+1]
                     - matrix.coeffsOffset()[patchi];
                }
                newSizes[globalPatchID] = intCoeffs.size();
                newSizes[nbrGlobalPatchID] = intCoeffs.size();
                newOffset[0] = 0;
                for (label i = 1; i < newOffset.size(); ++i)
                {
                    newOffset[i] = newOffset[i-1] + newSizes[i-1];
                }
            }

            // New offsets and resize internal and boundary coefficients
            matrix.coeffsOffset() = newOffset;
            matrix.internalCoeffs().resize_nocopy(matrix.coeffsOffset().last());
            matrix.boundaryCoeffs().resize_nocopy(matrix.coeffsOffset().last());

            // Insert new coefficients
            for (label patchi = 0; patchi < matrix.nPatches(); ++patchi)
            {
                if (patchi == globalPatchID)
                {
                    matrix.internalCoeffs(globalPatchID) =
                        intCoeffs*pTraits<Type>::one;
                    matrix.boundaryCoeffs(globalPatchID) =
                        bndCoeffs*pTraits<Type>::one;
                }
                else if (patchi == nbrGlobalPatchID)
                {
                    matrix.internalCoeffs(nbrGlobalPatchID) =
                        intCoeffs*pTraits<Type>::one;
                    matrix.boundaryCoeffs(nbrGlobalPatchID) =
                        bndCoeffs*pTraits<Type>::one;
                }
                else
                {
                    // Existing patch from old coefficients
                    const label oldOffset = oldCoeffsOffset[patchi];
                    matrix.internalCoeffs(patchi) =
                        SubList<Type>
                        (
                            oldIntCoeffs,
                            oldCoeffsOffset[patchi+1] - oldOffset,
                            oldOffset
                        );
                    matrix.boundaryCoeffs(patchi) =
                        SubList<Type>
                        (
                            oldBndCoeffs,
                            oldCoeffsOffset[patchi+1] - oldOffset,
                            oldOffset
                        );
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::cyclicAMIFvPatchField<Type>::coeffs
(
    fvMatrix<Type>& matrix,
    const Field<scalar>& coeffs,
    const label mat
) const
{
    const auto& matrixAssembly = matrix.lduMeshAssembly();
    const label index(this->patch().index());

    const label nSubFaces(matrixAssembly.cellBoundMap()[mat][index].size());

    auto tmapCoeffs = tmp<Field<scalar>>::New(nSubFaces, Zero);
    auto& mapCoeffs = tmapCoeffs.ref();

    const scalarListList& srcWeight =
        cyclicAMIPatch_.cyclicAMIPatch().AMI().srcWeights();

    label subFaceI = 0;
    forAll(*this, faceI)
    {
        const scalarList& w = srcWeight[faceI];
        for(label i=0; i<w.size(); i++)
        {
            const label localFaceId =
                matrixAssembly.facePatchFaceMap()[mat][index][subFaceI];
            mapCoeffs[subFaceI] = w[i]*coeffs[localFaceId];
            subFaceI++;
        }
    }

    return tmapCoeffs;
}


template<class Type>
template<class Type2>
void Foam::cyclicAMIFvPatchField<Type>::collectStencilData
(
    const refPtr<mapDistribute>& mapPtr,
    const labelListList& stencil,
    const Type2& data,
    List<Type2>& expandedData
)
{
    expandedData.resize_nocopy(stencil.size());
    if (mapPtr)
    {
        Type2 work(data);
        mapPtr().distribute(work);

        forAll(stencil, facei)
        {
            const auto& slots = stencil[facei];
            expandedData[facei].push_back
            (
                UIndirectList<typename Type2::value_type>(work, slots)
            );
        }
    }
    else
    {
        forAll(stencil, facei)
        {
            const auto& slots = stencil[facei];
            expandedData[facei].push_back
            (
                UIndirectList<typename Type2::value_type>(data, slots)
            );
        }
    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    fvPatchField<Type>::writeValueEntry(os);

    if (patchNeighbourFieldPtr_)
    {
        patchNeighbourFieldPtr_->writeEntry("neighbourValue", os);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=(ptf);

    //Pout<< "cyclicAMIFvPatchField::operator= :"
    //    << " field:" << this->internalField().name()
    //    << " patch:" << this->patch().name()
    //    << " copying from field:" << ptf.internalField().name()
    //    << endl;

    const auto* cycPtr = isA<this_bctype>(ptf);

    if
    (
        cycPtr
     && cycPtr->patchNeighbourFieldPtr_
     && cycPtr->patchNeighbourFieldPtr_->size() == this->size()
    )
    {
        if (!patchNeighbourFieldPtr_)
        {
            patchNeighbourFieldPtr_ = autoPtr<Field<Type>>::New();
        }

        // Copy values
        *patchNeighbourFieldPtr_ = *(cycPtr->patchNeighbourFieldPtr_);
    }
    else
    {
        patchNeighbourFieldPtr_.reset(nullptr);
    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator==(ptf);

    //Pout<< "cyclicAMIFvPatchField::operator== :"
    //    << " field:" << this->internalField().name()
    //    << " patch:" << this->patch().name()
    //    << " copying from field:" << ptf.internalField().name()
    //    << endl;

    const auto* cycPtr = isA<this_bctype>(ptf);

    if
    (
        cycPtr
     && cycPtr->patchNeighbourFieldPtr_
     && cycPtr->patchNeighbourFieldPtr_->size() == this->size()
    )
    {
        if (!patchNeighbourFieldPtr_)
        {
            patchNeighbourFieldPtr_ = autoPtr<Field<Type>>::New();
        }

        // Copy values
        *patchNeighbourFieldPtr_ = *(cycPtr->patchNeighbourFieldPtr_);
    }
    else
    {
        patchNeighbourFieldPtr_.reset(nullptr);
    }
}


// ************************************************************************* //
