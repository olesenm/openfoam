/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2026 UKRI-STFC.
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

    Some parts here added have been inspired by the work done at the following
    repository https://github.com/SimFlowCFD/RapidCFD-dev
\*---------------------------------------------------------------------------*/

#include "GAMGInterface.H"
#ifdef FOAM_OFFLOAD
#include "GAMGAgglomeration.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGInterface, 0);
    defineRunTimeSelectionTable(GAMGInterface, lduInterface);
    defineRunTimeSelectionTable(GAMGInterface, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGInterface::GAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    index_(index),
    coarseInterfaces_(coarseInterfaces),
    faceCells_(is),
    faceRestrictAddressing_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef FOAM_OFFLOAD
void Foam::GAMGInterface::updateAddressing()
{
    labelList restrictTmp(faceRestrictAddressing_);

    GAMGAgglomeration::createSort
    (
        restrictTmp,
        faceRestrictSortAddressing_
    );

    GAMGAgglomeration::createTarget
    (
        restrictTmp,
        faceRestrictSortAddressing_,
        faceRestrictTargetAddressing_,
        faceRestrictTargetStartAddressing_
    );

    labelList sort(faceCells_.size());
    GAMGAgglomeration::createSort
    (
        faceCells_,
        cellFaces_
    );

    GAMGAgglomeration::createTarget
    (
        faceCells_,
        cellFaces_,
        sortCells_,
        cellFacesStart_
    );
}
#endif

void Foam::GAMGInterface::combine(const GAMGInterface& coarseGi)
{
    const labelList& coarseFra = coarseGi.faceRestrictAddressing_;

    forAll(faceRestrictAddressing_, ffi)
    {
        faceRestrictAddressing_[ffi] = coarseFra[faceRestrictAddressing_[ffi]];
    }

    faceCells_ = coarseGi.faceCells_;

    #ifdef FOAM_OFFLOAD
    updateAddressing();
    #endif
}


Foam::tmp<Foam::labelField> Foam::GAMGInterface::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return interfaceInternalField<label>(internalData);
}


Foam::tmp<Foam::labelField> Foam::GAMGInterface::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    return interfaceInternalField<label>(internalData, faceCells);
}


Foam::tmp<Foam::scalarField> Foam::GAMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    auto tcoarseCoeffs = tmp<scalarField>::New(size(), Zero);
    auto& coarseCoeffs = tcoarseCoeffs.ref();

    if (fineCoeffs.size() != faceRestrictAddressing_.size())
    {
        FatalErrorInFunction
            << "Size of coefficients " << fineCoeffs.size()
            << " does not correspond to the size of the restriction "
            << faceRestrictAddressing_.size()
            << abort(FatalError);
    }
    if (debug && max(faceRestrictAddressing_) > size())
    {
        FatalErrorInFunction
            << "Face restrict addressing addresses outside of coarse interface"
            << " size. Max addressing:" << max(faceRestrictAddressing_)
            << " coarse size:" << size()
            << abort(FatalError);
    }

    #ifdef FOAM_OFFLOAD
    labelRange coarseFaceRange(faceRestrictTargetStartAddressing_.size() - 1);

    apply_policy(coarseFaceRange.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            coarseFaceRange.begin(),
            coarseFaceRange.end(),
            [&](const label coarseFacei)
            {
                const label start =
                    faceRestrictTargetStartAddressing_[coarseFacei];
                const label end =
                    faceRestrictTargetStartAddressing_[coarseFacei + 1];

                auto sum = 0.0;

                for (label i = start; i < end; ++i)
                {
                    const label fineFacei = faceRestrictSortAddressing_[i];

                    sum += fineCoeffs[fineFacei];
                }

                const label targetFacei =
                    faceRestrictTargetAddressing_[coarseFacei];

                coarseCoeffs[targetFacei] = sum;
            }
        );
    });
    #else
    forAll(faceRestrictAddressing_, ffi)
    {
        coarseCoeffs[faceRestrictAddressing_[ffi]] += fineCoeffs[ffi];
    }
    #endif

    return tcoarseCoeffs;
}


void Foam::GAMGInterface::write(Ostream& os) const
{
    os  << faceCells_ << token::SPACE << faceRestrictAddressing_;
}


// ************************************************************************* //
