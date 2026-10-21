/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

    Some parts here added have been inspired by the work done at the following
    repository https://github.com/SimFlowCFD/RapidCFD-dev

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "SubField.H"
#include "calculatedFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "coupledFvPatchFields.H"
#include "IndirectList.H"
#include "UniformList.H"

#include "cyclicFvPatchField.H"
#include "cyclicAMIFvPatchField.H"
#include "cyclicACMIFvPatchField.H"

#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    #ifdef FOAM_OFFLOAD
    const labelUList& sort,
    const labelUList& sortStart,
    #endif
    const UList<Type2>& pf,
    Field<Type2>& intf
) const
{
    #ifdef FOAM_OFFLOAD
    if (sort.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << sort.size()
            << ") and field (" << pf.size() << ") are different sizes" << endl
            << abort(FatalError);
    }

    labelRange indices(sortStart.size() - 1);

    auto neiFun = [&pf](const label, const label facei)
    {
        return pf[facei];
    };

    auto neiOp = [](Type2& result, const Type2& value)
    {
        result += value;
    };

    apply_policy(indices.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            indices.begin(),
            indices.end(),
            patchCellReduceFunctor<Type2, decltype(neiFun), decltype(neiOp)>
            (
                intf,
                sortStart,
                sort,
                addr,
                neiFun,
                neiOp
            )
        );
    });
    #else
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << addr.size()
            << ") and field (" << pf.size() << ") are different sizes" << endl
            << abort(FatalError);
    }

    forAll(addr, facei)
    {
        intf[addr[facei]] += pf[facei];
    }
    #endif
}


#ifdef FOAM_OFFLOAD
template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const labelUList& sort,
    const labelUList& sortStart,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, sort, sortStart, tpf(), intf);
    tpf.clear();
}
#else
template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}
#endif

template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    #ifdef FOAM_OFFLOAD
    const labelUList& sort,
    const labelUList& sortStart,
    #endif
    const UList<Type2>& pf,
    Field<Type2>& intf
) const
{
    #ifdef FOAM_OFFLOAD
    if (sort.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << sort.size()
            << ") and field (" << pf.size() << ") are different sizes" << endl
            << abort(FatalError);
    }

    auto neiFun = [&pf](const label, const label facei)
    {
        return pf[facei];
    };

    auto neiOp = [](Type2& result, const Type2& value)
    {
        result -= value;
    };

    labelRange indices(pf.size());

    apply_policy(pf.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            indices.begin(),
            indices.end(),
            patchCellReduceFunctor<Type2, decltype(neiFun), decltype(neiOp)>
            (
                intf,
                sortStart,
                sort,
                addr,
                neiFun,
                neiOp
            )
        );
    });
    #else
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << addr.size()
            << ") and field (" << pf.size() << ") are different sizes" << endl
            << abort(FatalError);
    }

    forAll(addr, facei)
    {
        intf[addr[facei]] -= pf[facei];
    }
    #endif
}


#ifdef FOAM_OFFLOAD
template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const labelUList& sort,
    const labelUList& sortStart,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, sort, sortStart, tpf(), intf);
    tpf.clear();
}
#else
template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}
#endif

template<class Type>
void Foam::fvMatrix<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt
) const
{
    for (label fieldi = 0; fieldi < nMatrices(); ++fieldi)
    {
        const auto& bpsi = this->psi(fieldi).boundaryField();

        forAll(bpsi, ptfi)
        {
            const label patchi = globalPatchID(fieldi, ptfi);

            if (patchi != -1)
            {
                addToInternalField
                (
                    #ifdef FOAM_OFFLOAD
                    lduAddr().patchSortCells(patchi),
                    lduAddr().patchSortAddr(patchi),
                    lduAddr().patchSortStartAddr(patchi),
                    #else
                    lduAddr().patchAddr(patchi),
                    #endif
                    internalCoeffs(patchi).component(solveCmpt),
                    diag
                );
            }
        }
    }
}


template<class Type>
void Foam::fvMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
    {
        const auto& bpsi = this->psi(fieldi).boundaryField();

        forAll(bpsi, ptfi)
        {
            const label patchi = globalPatchID(fieldi, ptfi);
            if (patchi != -1)
            {
                addToInternalField
                (
                    #ifdef FOAM_OFFLOAD
                    lduAddr().patchSortCells(patchi),
                    lduAddr().patchSortAddr(patchi),
                    lduAddr().patchSortStartAddr(patchi),
                    #else
                    lduAddr().patchAddr(patchi),
                    #endif
                    cmptAv(internalCoeffs(patchi)),
                    diag
                );
            }
        }
    }
}


template<class Type>
void Foam::fvMatrix<Type>::addBoundarySource
(
    Field<Type>& source,
    const bool couples
) const
{
    for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
    {
        const auto& bpsi = this->psi(fieldi).boundaryField();

        forAll(bpsi, ptfi)
        {
            const fvPatchField<Type>& ptf = bpsi[ptfi];

            const label patchi = globalPatchID(fieldi, ptfi);

            if (patchi != -1)
            {
                const SubField<Type> pbc = boundaryCoeffs(patchi);

                if (!ptf.coupled())
                {
                    addToInternalField
                    (
                        #ifdef FOAM_OFFLOAD
                        lduAddr().patchSortCells(patchi),
                        lduAddr().patchSortAddr(patchi),
                        lduAddr().patchSortStartAddr(patchi),
                        #else
                        lduAddr().patchAddr(patchi),
                        #endif
                        pbc,
                        source
                    );
                }
                else if (couples)
                {
                    const tmp<Field<Type>> tpnf = ptf.patchNeighbourField();
                    const Field<Type>& pnf = tpnf();

                    #ifdef FOAM_OFFLAOD
                    const labelUList& patchCells =
                        lduAddr().patchSortCells(patchi);

                    const labelUList& sort =
                        lduAddr().patchSortAddr(patchi);

                    const labelUList& sortStart =
                        lduAddr().patchSortStartAddr(patchi);

                    auto neiFun = [&pbc, &pnf](const label, const label facei)
                    {
                        return cmptMultiply(pbc[facei], pnf[facei]);
                    };

                    auto neiOp = [](Type& result, const Type& value)
                    {
                        result += value;
                    };

                    labelRange indices(patchCells.size());

                    apply_policy(indices.size(), [&](auto policy)
                    {
                        std::for_each
                        (
                            policy,
                            indices.begin(),
                            indices.end(),
                            patchCellReduceFunctor<Type, decltype(neiFun),
                                decltype(neiOp)>
                            (
                                source,
                                sortStart,
                                sort,
                                patchCells,
                                neiFun,
                                neiOp
                            )
                        );
                    });
                    #else
                    const labelUList& addr = lduAddr().patchAddr(patchi);

                    forAll(addr, facei)
                    {
                        source[addr[facei]] +=
                            cmptMultiply(pbc[facei], pnf[facei]);
                    }
                    #endif
                }
            }
        }
    }
}


template<class Type>
template<template<class> class ListType>
void Foam::fvMatrix<Type>::setValuesFromList
(
    const labelUList& cellLabels,
    const ListType<Type>& values
)
{
    const fvMesh& mesh = psi_.mesh();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    scalarField& Diag = diag();
    Field<Type>& psi =
        const_cast
        <
            GeometricField<Type, fvPatchField, volMesh>&
        >(psi_).primitiveFieldRef();

    // Following actions:
    // - adjust local field psi
    // - set local matrix to be diagonal (so adjust source)
    //      - cut connections to neighbours
    // - make (on non-adjusted cells) contribution explicit

    #ifdef FOAM_OFFLOAD
    // Expand values into a full cell-indexed field for O(1) neighbour lookup.
    // cellMask uses label (not bool) to avoid std::vector<bool> bit-packing
    // issues under par_unseq / GPU offload.
    const label nCells = psi_.size();
    const label nInternalFaces = mesh.nInternalFaces();

    Field<Type>   cellValues(nCells, Zero);
    List<label>   cellMask(nCells, 0);

    labelRange indices(cellLabels.size());

    apply_policy(cellLabels.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            indices.begin(),
            indices.end(),
            [&](const label i)
            {
                const label celli = cellLabels[i];
                cellMask[celli]   = 1;
                cellValues[celli] = values[i];
                psi[celli]        = values[i];
                source_[celli]    = values[i] * Diag[celli];
            }
        );
    });

    // Following actions for the matrix connectivity:
    // - adjust source on non-fixed cells that neighbour a fixed cell
    // - zero face coefficients on faces touching any fixed cell
    // - zero boundary coefficients for patches adjacent to fixed cells

    if (symmetric() || asymmetric())
    {
        // Per-face masks: 1 if that face's owner/neighbour is a fixed cell.
        // Used to (a) steer the source adjustment and (b) zero coefficients.
        List<label> ownMask(nInternalFaces, 0);
        List<label> neiMask(nInternalFaces, 0);

        labelRange faceIndices(nInternalFaces);

        apply_policy(nInternalFaces, [&](auto policy)
        {
            std::for_each
            (
                policy,
                faceIndices.begin(),
                faceIndices.end(),
                [&](const label facei)
                {
                    ownMask[facei] = cellMask[own[facei]];
                    neiMask[facei] = cellMask[nei[facei]];
                }
            );
        });

        // Adjust source for every non-fixed cell that has at least one
        // fixed neighbour.  Each cell writes only to source_[celli], so
        // there are no write races.
        //
        // LDU addressing (direct port of the RapidCFD functor):
        //   ownerStart[celli]   .. ownerStart[celli+1]   - faces where celli
        //                                                   is the owner; the
        //                                                   other cell is nei.
        //   losortStart[celli]  .. losortStart[celli+1]  - faces (via losort)
        //                                                   where celli is the
        //                                                   neighbour; other
        //                                                   cell is own.
        //
        //   symmetric:  lower() == upper()
        //   asymmetric: owner side uses lower, neighbour side uses upper

        const lduAddressing& addr = lduAddr();
        const labelUList& ownStart    = addr.ownerStartAddr();
        const labelUList& losortStart = addr.losortStartAddr();
        const labelUList& losort      = addr.losortAddr();

        const scalarField& Upper = upper();
        const scalarField& Lower = symmetric() ? upper() : lower();

        labelRange cellIndices(nCells);

        apply_policy(nCells, [&](auto policy)
        {
            std::for_each
            (
                policy,
                cellIndices.begin(),
                cellIndices.end(),
                [&](const label celli)
                {
                    if (cellMask[celli])
                        return;

                    // Faces where celli is the owner: fixed cell would be nei.
                    const label oStart = ownStart[celli];
                    const label oEnd   = ownStart[celli + 1];

                    for (label facei = oStart; facei < oEnd; ++facei)
                    {
                        if (neiMask[facei])
                        {
                            source_[celli] -=
                                Upper[facei] * cellValues[nei[facei]];
                        }
                    }

                    // Faces where celli is the neighbour: fixed cell is own.
                    const label nStart = losortStart[celli];
                    const label nEnd   = losortStart[celli + 1];

                    for (label i = nStart; i < nEnd; ++i)
                    {
                        const label facei = losort[i];

                        if (ownMask[facei])
                        {
                            source_[celli] -=
                                Lower[facei] * cellValues[own[facei]];
                        }
                    }
                }
            );
        });

        // Zero upper (and lower for asymmetric) on every face that touches
        // a fixed cell.  Each face has its own slot - no races.

        scalarField& Upper_rw = upper();

        apply_policy(nInternalFaces, [&](auto policy)
        {
            std::for_each
            (
                policy,
                faceIndices.begin(),
                faceIndices.end(),
                [&](const label facei)
                {
                    if (ownMask[facei] || neiMask[facei])
                    {
                        Upper_rw[facei] = 0.0;
                    }
                }
            );
        });

        if (asymmetric())
        {
            scalarField& Lower_rw = lower();

            apply_policy(nInternalFaces, [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    faceIndices.begin(),
                    faceIndices.end(),
                    [&](const label facei)
                    {
                        if (ownMask[facei] || neiMask[facei])
                        {
                            Lower_rw[facei] = 0.0;
                        }
                    }
                );
            });
        }

        // Boundary patches: zero internalCoeffs / boundaryCoeffs for faces
        // whose owner cell is in the fixed set.

        forAll(mesh.boundaryMesh(), patchi)
        {
            SubField<Type> iCoeffs = internalCoeffs(patchi);
            SubField<Type> bCoeffs = boundaryCoeffs(patchi);

            if (!iCoeffs.size())
            {
                continue;
            }

            const labelUList& pcells = mesh.boundary()[patchi].faceCells();
            const label nPatchFaces = pcells.size();

            const labelRange patchIndices(nPatchFaces);

            // Resolve the PtrList element references on the host. Indexing
            // internalCoeffs_[patchi] inside the device lambda would invoke
            // UPtrList::operator[] -> at() -> FatalError, which nvc++ cannot
            // compile for the device. Capture raw pointers instead.
            Type* const intPtr = iCoeffs.data();
            Type* const bndPtr = bCoeffs.data();
            const label* const pcellsPtr = pcells.cdata();
            const label* const maskPtr = cellMask.cdata();

            apply_policy(nPatchFaces, [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    patchIndices.begin(),
                    patchIndices.end(),
                    [=](const label patchFacei)
                    {
                        if (maskPtr[pcellsPtr[patchFacei]])
                        {
                            intPtr[patchFacei] = Zero;
                            bndPtr[patchFacei] = Zero;
                        }
                    }
                );
            });
        }
    }
    #else
    // Keep cpu code since does not require extra memory.
    if (symmetric() || asymmetric())
    {
        const cellList& cells = mesh.cells();
        forAll(cellLabels, i)
        {
            const label celli = cellLabels[i];
            const Type& value = values[i];

            for (const label facei : cells[celli])
            {
                const label patchi = mesh.boundaryMesh().patchID(facei);

                if (patchi == -1)
                {
                    if (symmetric())
                    {
                        if (celli == own[facei])
                        {
                            source_[nei[facei]] -= upper()[facei]*value;
                        }
                        else
                        {
                            source_[own[facei]] -= upper()[facei]*value;
                        }

                        upper()[facei] = 0.0;
                    }
                    else
                    {
                        if (celli == own[facei])
                        {
                            source_[nei[facei]] -= lower()[facei]*value;
                        }
                        else
                        {
                            source_[own[facei]] -= upper()[facei]*value;
                        }

                        upper()[facei] = 0.0;
                        lower()[facei] = 0.0;
                    }
                }
                else
                {
                    SubField<Type> iCoeffs = internalCoeffs(patchi);
                    if (iCoeffs.size())
                    {
                        const label patchFacei =
                            mesh.boundaryMesh()[patchi].whichFace(facei);

                        iCoeffs[patchFacei] = Zero;
                        boundaryCoeffs(patchi)[patchFacei] = Zero;
                    }
                }
            }
        }
    }

    // Note: above loop might have affected source terms on adjusted cells
    // so make sure to adjust them afterwards
    forAll(cellLabels, i)
    {
        const label celli = cellLabels[i];
        const Type& value = values[i];

        psi[celli] = value;
        source_[celli] = value*Diag[celli];
    }
    #endif
}



template<class Type>
bool Foam::fvMatrix<Type>::checkImplicit(const label fieldi)
{
    const auto& bpsi = this->psi(fieldi).boundaryField();

    word idName;
    forAll(bpsi, patchi)
    {
        if (bpsi[patchi].useImplicit())
        {
            if (debug)
            {
                Pout<< "fvMatrix<Type>::checkImplicit "
                    << " field:" << this->psi(fieldi).name()
                    << " on mesh:"
                    << this->psi(fieldi).mesh().name()
                    << " patch:" << bpsi[patchi].patch().name()
                    << endl;
            }

            idName += Foam::name(patchi);
            useImplicit_ = true;
        }
    }

    if (useImplicit_)
    {
        lduAssemblyName_ = word("lduAssembly") + idName;
    }

    return !idName.empty();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvMatrix<Type>::fvMatrix
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    useImplicit_(false),
    lduAssemblyName_(),
    nMatrix_(0),
    dimensions_(ds),
    source_(psi.size(), Zero),
    internalCoeffs_(psi.mesh().nBoundaryFaces(), Zero),
    boundaryCoeffs_(psi.mesh().nBoundaryFaces(), Zero),
    coeffsOffset_(psi.mesh().boundary().size() + 1, 0)
{
    DebugInFunction
        << "Constructing fvMatrix<Type> for field " << psi_.name() << endl;

    checkImplicit();

    // Populate offset array for per-patch coefficient access
    forAll(psi.mesh().boundary(), patchi)
    {
        const label patchSize = psi.mesh().boundary()[patchi].size();
        coeffsOffset_[patchi+1] = coeffsOffset_[patchi] + patchSize;
    }

    auto& psiRef = this->psi(0);
    const label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryFieldRef().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const fvMatrix<Type>& fvm)
:
    lduMatrix(fvm),
    psi_(fvm.psi_),
    useImplicit_(fvm.useImplicit_),
    lduAssemblyName_(fvm.lduAssemblyName_),
    nMatrix_(fvm.nMatrix_),
    dimensions_(fvm.dimensions_),
    source_(fvm.source_),
    internalCoeffs_(fvm.internalCoeffs_),
    boundaryCoeffs_(fvm.boundaryCoeffs_),
    coeffsOffset_(fvm.coeffsOffset_)
{
    DebugInFunction
        << "Copying fvMatrix<Type> for field " << psi_.name() << endl;

    if (fvm.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            *(fvm.faceFluxCorrectionPtr_)
        );
    }
}


template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const tmp<fvMatrix<Type>>& tmat)
:
    lduMatrix(tmat.constCast(), tmat.movable()),
    psi_(tmat().psi_),
    useImplicit_(tmat().useImplicit_),
    lduAssemblyName_(tmat().lduAssemblyName_),
    nMatrix_(tmat().nMatrix_),
    dimensions_(tmat().dimensions_),
    source_(tmat.constCast().source_, tmat.movable()),
    internalCoeffs_(tmat.constCast().internalCoeffs_, tmat.movable()),
    boundaryCoeffs_(tmat.constCast().boundaryCoeffs_, tmat.movable()),
    coeffsOffset_(tmat.constCast().coeffsOffset_, tmat.movable())
{
    DebugInFunction
        << "Copy/move fvMatrix<Type> for field " << psi_.name() << endl;

    if (tmat().faceFluxCorrectionPtr_)
    {
        if (tmat.movable())
        {
            faceFluxCorrectionPtr_ =
                std::move(tmat.constCast().faceFluxCorrectionPtr_);
        }
        else if (tmat().faceFluxCorrectionPtr_)
        {
            faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
            (
                *(tmat().faceFluxCorrectionPtr_)
            );
        }
    }

    tmat.clear();
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fvMatrix<Type>::~fvMatrix()
{
    DebugInFunction
        << "Destroying fvMatrix<Type> for field " << psi_.name() << endl;

    subMatrices_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::setInterfaces
(
    lduInterfaceFieldPtrsList& interfaces,
    PtrDynList<lduInterfaceField>& newInterfaces
)
{
    interfaces.setSize(coeffsOffset_.size() - 1);
    for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
    {
        const auto& bpsi = this->psi(fieldi).boundaryField();
        lduInterfaceFieldPtrsList fieldInterfaces(bpsi.scalarInterfaces());

        forAll (fieldInterfaces, patchi)
        {
            label globalPatchID = lduMeshPtr()->patchMap()[fieldi][patchi];

            if (globalPatchID != -1)
            {
                if (fieldInterfaces.set(patchi))
                {
                    if (isA<cyclicLduInterfaceField>(bpsi[patchi]))
                    {
                        newInterfaces.append
                        (
                            new cyclicFvPatchField<Type>
                            (
                                refCast<const fvPatch>
                                (
                                    lduMeshPtr()->interfaces()[globalPatchID]
                                ),
                                bpsi[patchi].internalField()
                            )
                        );
                        interfaces.set(globalPatchID, &newInterfaces.last());

                    }
                    else if (isA<cyclicAMILduInterfaceField>(bpsi[patchi]))
                    {
                        newInterfaces.append
                        (
                            new cyclicAMIFvPatchField<Type>
                            (
                                refCast<const fvPatch>
                                (
                                    lduMeshPtr()->interfaces()[globalPatchID]
                                ),
                                bpsi[patchi].internalField()
                            )
                        );
                        interfaces.set(globalPatchID, &newInterfaces.last());
                    }
                    else if (isA<cyclicACMILduInterfaceField>(bpsi[patchi]))
                    {
                        newInterfaces.append
                        (
                            new cyclicACMIFvPatchField<Type>
                            (
                                refCast<const fvPatch>
                                (
                                    lduMeshPtr()->interfaces()[globalPatchID]
                                ),
                                bpsi[patchi].internalField()
                            )
                        );
                        interfaces.set(globalPatchID, &newInterfaces.last());
                    }
                    else
                    {
                        interfaces.set(globalPatchID, &fieldInterfaces[patchi]);
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMatrix<Type>::mapContributions
(
    label fieldi,
    const FieldField<Field, Type>& fluxContrib,
    FieldField<Field, Type>& contrib,
    bool internal
) const
{
    const lduPrimitiveMeshAssembly* ptr = lduMeshPtr();

    const labelList& patchMap = ptr->patchMap()[fieldi];

    forAll(contrib, patchi)
    {
        const label globalPtchId = patchMap[patchi];

        if (globalPtchId != -1)
        {
            // Cache contrib before overwriting
            const Field<Type> saveContrib(fluxContrib[globalPtchId]);
            contrib[patchi].setSize(psi_.boundaryField()[patchi].size()),
            contrib[patchi] = pTraits<Type>::zero;

            if (internal)
            {
                contrib[patchi] =
                    cmptMultiply
                    (
                        saveContrib,
                        psi_.boundaryField()[patchi].patchInternalField()
                    );
            }
            else
            {
                if (this->psi(fieldi).boundaryField()[patchi].coupled())
                {
                    contrib[patchi] =
                        cmptMultiply
                        (
                            saveContrib,
                            psi_.boundaryField()[patchi].patchNeighbourField()
                        );
                }
            }
        }
        else if (globalPtchId == -1)
        {
            const polyPatch& pp =
                this->psi(fieldi).mesh().boundaryMesh()[patchi];

            if (pp.masterImplicit())
            {
                label virtualPatch =
                    ptr->patchLocalToGlobalMap()[fieldi][patchi];

                const label nbrPatchId = pp.neighbPolyPatchID();

                // Copy contrib before overwriting
                const Field<Type> saveContrib(fluxContrib[virtualPatch]);

                Field<Type>& coeffs = contrib[patchi];
                Field<Type>& nbrCoeffs = contrib[nbrPatchId];

                coeffs.setSize(psi_.boundaryField()[patchi].size());
                nbrCoeffs.setSize(psi_.boundaryField()[nbrPatchId].size());

                coeffs = pTraits<Type>::zero;
                nbrCoeffs = pTraits<Type>::zero;

                // nrb cells
                const labelList& nbrCellIds =
                    ptr->cellBoundMap()[fieldi][patchi];

                const labelList& cellIds =
                    ptr->cellBoundMap()[fieldi][nbrPatchId];

                const GeometricField<Type, fvPatchField, volMesh>& psi =
                    this->psi(fieldi);

                forAll(saveContrib, subFaceI)
                {
                    const label faceId =
                        ptr->facePatchFaceMap()[fieldi][patchi][subFaceI];
                    const label nbrFaceId =
                        ptr->facePatchFaceMap()[fieldi][nbrPatchId][subFaceI];

                    const label nbrCellId = nbrCellIds[subFaceI];
                    const label cellId = cellIds[subFaceI];

                    if (internal)
                    {
                        coeffs[faceId] +=
                            cmptMultiply(saveContrib[subFaceI], psi[cellId]);

                        nbrCoeffs[nbrFaceId] +=
                            cmptMultiply(saveContrib[subFaceI], psi[nbrCellId]);
                    }
                    else //boundary
                    {
                        coeffs[faceId] +=
                            cmptMultiply(saveContrib[subFaceI], psi[nbrCellId]);

                        nbrCoeffs[nbrFaceId] +=
                            cmptMultiply(saveContrib[subFaceI], psi[cellId]);
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMatrix<Type>::setBounAndInterCoeffs()
{
    // Count number of interfaces and offsets
    label interfaceI = 0;
    for (label i = 0; i < nMatrices(); i++)
    {
        interfaceI += this->psi(i).mesh().boundary().size();
    }

    // Sizes in terms of new ordering
    labelList sizes(interfaceI, 0);
    for (label i = 0; i < nMatrices(); i++)
    {
        const auto& boundary = this->psi(i).mesh().boundary();
        forAll(boundary, patchI)
        {
            label globalPatchId = lduMeshPtr()->patchMap()[i][patchI];
            if (globalPatchId != -1)
            {
                sizes[globalPatchId] = boundary[patchI].size();
            }
            else
            {
                globalPatchId =
                    lduMeshPtr()->patchLocalToGlobalMap()[i][patchI];
                sizes[globalPatchId] = boundary[patchI].size();
            }
        }
    }

    // Sum up to get offsets
    labelList offsets(sizes.size()+1);
    offsets[0] = 0;
    for (label i = 1; i <= sizes.size(); i++)
    {
        offsets[i] = offsets[i-1] + sizes[i-1];
    }

    // Construct new internal/boundary in new order:
    // - global patches first (according to patchMap())
    // - implicit patches next (according to patchLocalToGlobalMap())

    Field<Type> newInternalCoeffs(offsets.last());
    Field<Type> newBoundaryCoeffs(offsets.last());

    // Preserved patches first
    for (label i = 0; i < nMatrices(); i++)
    {
        const auto& boundary = this->psi(i).mesh().boundary();
        const auto& matrix = this->matrix(i);
        forAll(boundary, patchI)
        {
            const label globalPatchId = lduMeshPtr()->patchMap()[i][patchI];
            if (globalPatchId != -1)
            {
                const auto srcInt = matrix.internalCoeffs(patchI);
                SubList<Type>
                (
                    newInternalCoeffs,
                    srcInt.size(),
                    offsets[globalPatchId]
                ) = srcInt;

                const auto srcBnd = matrix.boundaryCoeffs(patchI);
                SubList<Type>
                (
                    newBoundaryCoeffs,
                    srcBnd.size(),
                    offsets[globalPatchId]
                ) = srcBnd;
            }
        }
    }

    // Implicit patches next
    for (label i = 0; i < nMatrices(); i++)
    {
        const auto& boundary = this->psi(i).mesh().boundary();
        const auto& matrix = this->matrix(i);
        forAll(boundary, patchI)
        {
            if (lduMeshPtr()->patchMap()[i][patchI] == -1)
            {
                const label globalPatchId =
                    lduMeshPtr()->patchLocalToGlobalMap()[i][patchI];

                const auto srcInt = matrix.internalCoeffs(patchI);
                SubList<Type>
                (
                    newInternalCoeffs,
                    srcInt.size(),
                    offsets[globalPatchId]
                ) = srcInt;

                const auto srcBnd = matrix.boundaryCoeffs(patchI);
                SubList<Type>
                (
                    newBoundaryCoeffs,
                    srcBnd.size(),
                    offsets[globalPatchId]
                ) = srcBnd;
            }
        }
    }

    internalCoeffs_ = std::move(newInternalCoeffs);
    boundaryCoeffs_ = std::move(newBoundaryCoeffs);
    coeffsOffset_ = std::move(offsets);

    // Pout<< "source:" << flatOutput(source_) << endl;
    // for (label patchi = 0; patchi < coeffsOffset_.size()-1; ++patchi)
    // {
    //     Pout<< "patch:" << patchi
    //         << " offset: " << coeffsOffset_[patchi]
    //         << " to " << coeffsOffset_[patchi+1] << endl;
    //     Pout<< "    internalCoeffs: " << flatOutput(internalCoeffs(patchi)) << endl;
    //     Pout<< "    boundaryCoeffs: " << flatOutput(boundaryCoeffs(patchi)) << endl;
    // }
}


template<class Type>
void Foam::fvMatrix<Type>::manipulateMatrix(direction cmp)
{
    for (label i=0; i < nMatrices(); ++i)
    {
        forAll(psi(i).boundaryField(), patchI)
        {
            label globalPatchId = lduMeshPtr()->patchMap()[i][patchI];

            if (globalPatchId == -1)
            {
                psi(i).boundaryFieldRef()[patchI].manipulateMatrix
                (
                    *this,
                    i,
                    cmp
                );
            }
        }
    }
}


template<class Type>
void Foam::fvMatrix<Type>::transferFvMatrixCoeffs()
{
    const labelListList& faceMap = lduMeshPtr()->faceMap();
    const labelList& cellMap = lduMeshPtr()->cellOffsets();

    label newFaces = lduMeshPtr()->lduAddr().upperAddr().size();
    label newCells = lduMeshPtr()->lduAddr().size();

    scalarField lowerAssemb(newFaces, Zero);
    scalarField upperAssemb(newFaces, Zero);
    scalarField diagAssemb(newCells, Zero);
    Field<Type> sourceAssemb(newCells, Zero);

    bool asymmetricAssemby = false;
    for (label i=0; i < nMatrices(); ++i)
    {
        if (matrix(i).asymmetric())
        {
            asymmetricAssemby = true;
        }
    }
    // Move append contents into intermediate list
    for (label i=0; i < nMatrices(); ++i)
    {
        if (asymmetricAssemby)
        {
            const scalarField lowerSub(matrix(i).lower());
            forAll(lowerSub, facei)
            {
                lowerAssemb[faceMap[i][facei]] = lowerSub[facei];
            }
        }

        const scalarField upperSub(matrix(i).upper());
        const scalarField diagSub(matrix(i).diag());
        const Field<Type> sourceSub(matrix(i).source());

        forAll(upperSub, facei)
        {
            upperAssemb[faceMap[i][facei]] = upperSub[facei];
        }

        forAll(diagSub, celli)
        {
            const label globalCelli = cellMap[i] + celli;
            diagAssemb[globalCelli] = diagSub[celli];
            sourceAssemb[globalCelli] = sourceSub[celli];
        }
    }

    if (asymmetricAssemby)
    {
        lower().setSize(newFaces, Zero);
        lower() = lowerAssemb;
    }
    upper().setSize(newFaces, Zero);
    upper() = upperAssemb;

    diag().setSize(newCells, Zero);
    diag() = diagAssemb;

    source().setSize(newCells, Zero);
    source() = sourceAssemb;
}


template<class Type>
Foam::lduPrimitiveMeshAssembly* Foam::fvMatrix<Type>::lduMeshPtr()
{
    return
    (
        psi_.mesh().thisDb().objectRegistry::template getObjectPtr
        <
            lduPrimitiveMeshAssembly
        > (lduAssemblyName_)
    );
}


template<class Type>
const Foam::lduPrimitiveMeshAssembly* Foam::fvMatrix<Type>::lduMeshPtr() const
{
    return
    (
        psi_.mesh().thisDb().objectRegistry::template cfindObject
        <
            lduPrimitiveMeshAssembly
        > (lduAssemblyName_)
    );
}


template<class Type>
void Foam::fvMatrix<Type>::createOrUpdateLduPrimitiveAssembly()
{
    lduPrimitiveMeshAssembly* ptr = lduMeshPtr();

    IOobject io
    (
        lduAssemblyName_,
        psi_.mesh().time().timeName(),
        psi_.mesh().thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::REGISTER
    );

    UPtrList<lduMesh> uMeshPtr(nMatrices());

    UPtrList<GeometricField<Type, fvPatchField, volMesh>>
        uFieldPtr(nMatrices());

    for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
    {
        const fvMesh& meshi = this->psi(fieldi).mesh();
        uMeshPtr.set
        (
            fieldi,
            &const_cast<fvMesh&>(meshi)
        );
        uFieldPtr.set(fieldi, &this->psi(fieldi));
    }

    if (!ptr)
    {
        lduPrimitiveMeshAssembly* lduAssemMeshPtr =
            new lduPrimitiveMeshAssembly(io, uMeshPtr);

        lduAssemMeshPtr->store();
        lduAssemMeshPtr->update(uFieldPtr);

        Info
            << "Creating lduPrimitiveAssembly: " << lduAssemblyName_ << endl;
    }
    else if
    (
        psi_.mesh().changing() && !psi_.mesh().upToDatePoints(*ptr)
    )
    {
        // Clear losortPtr_, ownerStartPtr_, losortStartPtr_
        ptr->lduAddr().clearOut();
        ptr->update(uFieldPtr);
        psi_.mesh().setUpToDatePoints(*ptr);

        Info
            << "Updating lduPrimitiveAssembly: " << lduAssemblyName_ << endl;
    }
    else
    {
        Info
            << "Using lduPrimitiveAssembly: " << lduAssemblyName_ << endl;
    }
}


template<class Type>
void Foam::fvMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const Type& value
)
{
    this->setValuesFromList(cellLabels, UniformList<Type>(value));
}


template<class Type>
void Foam::fvMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const UList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}


template<class Type>
void Foam::fvMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const UIndirectList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}


template<class Type>
void Foam::fvMatrix<Type>::setReference
(
    const label celli,
    const Type& value,
    const bool forceReference
)
{
    if ((forceReference || psi_.needReference()) && celli >= 0)
    {
        source()[celli] += diag()[celli]*value;
        diag()[celli] += diag()[celli];
    }
}


template<class Type>
void Foam::fvMatrix<Type>::setReferences
(
    const labelUList& cellLabels,
    const Type& value,
    const bool forceReference
)
{
    if (forceReference || psi_.needReference())
    {
        labelRange cellRange(0, cellLabels.size());

        Field<Type>& src = source();
        scalarField& d = diag();

        Type* srcPtr = src.begin();
        scalar* dPtr = d.begin();

        apply_policy(cellLabels.size(), [&](auto policy)
        {
            std::for_each
            (
                policy,
                cellRange.begin(),
                cellRange.end(),
                [=, &cellLabels](const label celli)
                {
                    const label cellId = cellLabels[celli];

                    if (cellId >= 0)
                    {
                        srcPtr[cellId] += dPtr[cellId]*value;
                        dPtr[cellId] += dPtr[cellId];
                    }
                }
            );
        });
    }
}


template<class Type>
void Foam::fvMatrix<Type>::setReferences
(
    const labelUList& cellLabels,
    const UList<Type>& values,
    const bool forceReference
)
{
    if (forceReference || psi_.needReference())
    {
        labelRange cellRange(0, cellLabels.size());

        Field<Type>& src = source();
        scalarField& d = diag();

        Type* srcPtr = src.begin();
        scalar* dPtr = d.begin();

        const Type* valuesPtr = values.begin();

        apply_policy(cellLabels.size(), [&](auto policy)
        {
            std::for_each
            (
                policy,
                cellRange.begin(),
                cellRange.end(),
                [=, &cellLabels](const label celli)
                {
                    const label cellId = cellLabels[celli];

                    if (cellId >= 0)
                    {
                        srcPtr[cellId] += dPtr[cellId]*valuesPtr[celli];
                        dPtr[cellId] += dPtr[cellId];
                    }
                }
            );
        });
    }
}


template<class Type>
void Foam::fvMatrix<Type>::addFvMatrix(fvMatrix& matrix)
{
    subMatrices_.append(matrix.clone());
    ++nMatrix_;

    if (dimensions_ != matrix.dimensions())
    {
        FatalErrorInFunction
            << "incompatible dimensions for matrix addition "
            << endl << "    "
            << "[" << dimensions_  << " ] "
            << " [" << matrix.dimensions() << " ]"
            << abort(FatalError);
    }

    for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
    {
        if (checkImplicit(fieldi))
        {
            break;
        }
    }

    internalCoeffs_.clear();
    boundaryCoeffs_.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    DebugInFunction
        << "Relaxing " << psi_.name() << " by " << alpha << endl;

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), Zero);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (ptf.size())
        {
            const SubField<Type> iCoeffs = internalCoeffs(patchi);

            if (ptf.coupled())
            {
                const SubField<Type> pCoeffs = boundaryCoeffs(patchi);

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                #ifdef FOAM_OFFLOAD
                const labelUList& patchSortCells =
                    lduAddr().patchSortCells(patchi);

                const labelUList& patchSortStart =
                    lduAddr().patchSortStartAddr(patchi);

                const labelUList& patchSort =
                    lduAddr().patchSortAddr(patchi);

                auto dFun = [iCoeffs](const label, const label facei)
                {
                    return component(iCoeffs[facei], 0);
                };

                auto sumOffFun = [pCoeffs](const label, const label facei)
                {
                    return mag(component(pCoeffs[facei], 0));
                };

                labelRange cells(patchSortCells.size());

                apply_policy(cells.size(), [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        cells.begin(),
                        cells.end(),
                        patchCellReduceFunctor
                        <
                            scalar,
                            decltype(dFun),
                            plusEqOp<scalar>
                        >
                        (
                            D,
                            patchSortStart,
                            patchSort,
                            patchSortCells,
                            dFun,
                            plusEqOp<scalar>()
                        )
                    );

                    std::for_each
                    (
                        policy,
                        cells.begin(),
                        cells.end(),
                        patchCellReduceFunctor
                        <
                            scalar,
                            decltype(sumOffFun),
                            plusEqOp<scalar>
                        >
                        (
                            sumOff,
                            patchSortStart,
                            patchSort,
                            patchSortCells,
                            sumOffFun,
                            plusEqOp<scalar>()
                        )
                    );
                });
                #else
                const labelUList& pa = lduAddr().patchAddr(patchi);

                forAll(pa, face)
                {
                    D[pa[face]] += component(iCoeffs[face], 0);
                    sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
                }
                #endif
            }
            else
            {
                // For non-coupled boundaries add the maximum magnitude diagonal
                // contribution to ensure stability
                #ifdef FOAM_OFFLOAD
                const labelUList& patchSortCells =
                    lduAddr().patchSortCells(patchi);

                const labelUList& patchSortStart =
                    lduAddr().patchSortStartAddr(patchi);

                const labelUList& patchSort =
                    lduAddr().patchSortAddr(patchi);

                auto faceOp = [iCoeffs](const label, const label facei)
                {
                    return cmptMax(cmptMag(iCoeffs[facei]));
                };

                labelRange cells(patchSortCells.size());

                apply_policy(cells.size(), [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        cells.begin(),
                        cells.end(),
                        patchCellReduceFunctor
                        <
                            scalar,
                            decltype(faceOp),
                            plusEqOp<scalar>
                        >
                        (
                            D,
                            patchSortStart,
                            patchSort,
                            patchSortCells,
                            faceOp,
                            plusEqOp<scalar>()
                        )
                    );
                });
                #else
                const labelUList& pa = lduAddr().patchAddr(patchi);

                forAll(pa, face)
                {
                    D[pa[face]] += cmptMax(cmptMag(iCoeffs[face]));
                }
                #endif
            }
        }
    }


    if (debug)
    {
        // Calculate amount of non-dominance.
        label nNon = 0;
        scalar maxNon = 0.0;
        scalar sumNon = 0.0;
        forAll(D, celli)
        {
            scalar d = (sumOff[celli] - D[celli])/mag(D[celli]);

            if (d > 0)
            {
                nNon++;
                maxNon = max(maxNon, d);
                sumNon += d;
            }
        }

        reduce(nNon, sumOp<label>(), UPstream::msgType(), psi_.mesh().comm());
        reduce
        (
            maxNon,
            maxOp<scalar>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );
        reduce
        (
            sumNon,
            sumOp<scalar>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );
        sumNon /= returnReduce
        (
            D.size(),
            sumOp<label>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );

        InfoInFunction
            << "Matrix dominance test for " << psi_.name() << nl
            << "    number of non-dominant cells   : " << nNon << nl
            << "    maximum relative non-dominance : " << maxNon << nl
            << "    average relative non-dominance : " << sumNon << nl
            << endl;
    }


    // Ensure the matrix is diagonally dominant...
    // Assumes that the central coefficient is positive and ensures it is
    labelRange cells(D.size());

    apply_policy(cells.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            cells.begin(),
            cells.end(),
            [&](const label celli)
            {
                D[celli] = max(mag(D[celli]), sumOff[celli]);
            }
        );
    });

    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (ptf.size())
        {
            const SubField<Type> iCoeffs = internalCoeffs(patchi);

            if (ptf.coupled())
            {
                #ifdef FOAM_OFFLOAD
                const labelUList& patchSortCells =
                    lduAddr().patchSortCells(patchi);

                const labelUList& patchSortStart =
                    lduAddr().patchSortStartAddr(patchi);

                const labelUList& patchSort =
                    lduAddr().patchSortAddr(patchi);

                auto faceOp = [iCoeffs](const label, const label facei)
                {
                    return component(iCoeffs[facei], 0);
                };

                labelRange cells(patchSortCells.size());

                apply_policy(cells.size(), [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        cells.begin(),
                        cells.end(),
                        patchCellReduceFunctor
                        <
                            scalar,
                            decltype(faceOp),
                            minusEqOp<scalar>
                        >
                        (
                            D,
                            patchSortStart,
                            patchSort,
                            patchSortCells,
                            faceOp,
                            minusEqOp<scalar>()
                        )
                    );
                });
                #else
                const labelUList& pa = lduAddr().patchAddr(patchi);

                forAll(pa, face)
                {
                    D[pa[face]] -= component(iCoeffs[face], 0);
                }
                #endif
            }
            else
            {
                #ifdef FOAM_OFFLOAD
                const labelUList& patchSortCells =
                    lduAddr().patchSortCells(patchi);

                const labelUList& patchSortStart =
                    lduAddr().patchSortStartAddr(patchi);

                const labelUList& patchSort =
                    lduAddr().patchSortAddr(patchi);

                auto faceOp = [iCoeffs](const label, const label facei)
                {
                    return cmptMin(iCoeffs[facei]);
                };

                labelRange cells(patchSortCells.size());

                apply_policy(cells.size(), [&](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        cells.begin(),
                        cells.end(),
                        patchCellReduceFunctor
                        <
                            scalar,
                            decltype(faceOp),
                            minusEqOp<scalar>
                        >
                        (
                            D,
                            patchSortStart,
                            patchSort,
                            patchSortCells,
                            faceOp,
                            minusEqOp<scalar>()
                        )
                    );
                });
                #else
                const labelUList& pa = lduAddr().patchAddr(patchi);

                forAll(pa, face)
                {
                    D[pa[face]] -= cmptMin(iCoeffs[face]);
                }
                #endif
            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.primitiveField();
}


template<class Type>
void Foam::fvMatrix<Type>::relax()
{
    word name = psi_.select
    (
        psi_.mesh().data().isFinalIteration()
    );

    scalar relaxCoeff = 0;

    if (psi_.mesh().relaxEquation(name, relaxCoeff))
    {
        relax(relaxCoeff);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::boundaryManipulate
(
    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& bFields
)
{
    forAll(bFields, patchi)
    {
        bFields[patchi].manipulateMatrix(*this);
    }
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::fvMatrix<Type>::D() const
{
    auto tdiag = tmp<scalarField>::New(diag());
    addCmptAvBoundaryDiag(tdiag.ref());
    return tdiag;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvMatrix<Type>::DD() const
{
    tmp<Field<Type>> tdiag(pTraits<Type>::one*diag());

    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (!ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                #ifdef FOAM_OFFLOAD
                lduAddr().patchSortCells(patchi),
                lduAddr().patchSortAddr(patchi),
                lduAddr().patchSortStartAddr(patchi),
                #else
                lduAddr().patchAddr(patchi),
                #endif
                internalCoeffs(patchi),
                tdiag.ref()
            );
        }
    }

    return tdiag;
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Type>::A() const
{
    auto tAphi = volScalarField::New
    (
        "A(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions_/psi_.dimensions()/dimVol,
        fvPatchFieldBase::extrapolatedCalculatedType()
    );

    tAphi.ref().primitiveFieldRef() = D()/psi_.mesh().V();
    tAphi.ref().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvMatrix<Type>::H() const
{
    auto tHphi = GeometricField<Type, fvPatchField, volMesh>::New
    (
        "H(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions_/dimVol,
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& Hphi = tHphi.ref();

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt(psi_.primitiveField().component(cmpt));

        scalarField boundaryDiagCmpt(psi_.size(), Zero);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        Hphi.primitiveFieldRef().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    Hphi.primitiveFieldRef() += lduMatrix::H(psi_.primitiveField()) + source_;
    addBoundarySource(Hphi.primitiveFieldRef());

    Hphi.primitiveFieldRef() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    typename Type::labelType validComponents
    (
        psi_.mesh().template validComponents<Type>()
    );

    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1)
        {
            Hphi.replace
            (
                cmpt,
                dimensionedScalar(Hphi.dimensions(), Zero)
            );
        }
    }

    return tHphi;
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Type>::H1() const
{
    auto tH1 = volScalarField::New
    (
        "H(1)",
        psi_.mesh(),
        dimensions_/(dimVol*psi_.dimensions()),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& H1_ = tH1.ref();

    H1_.primitiveFieldRef() = lduMatrix::H1();

    forAll(psi_.boundaryField(), patchi)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchi];

        if (ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                #ifdef FOAM_OFFLOAD
                lduAddr().patchSortCells(patchi),
                lduAddr().patchSortAddr(patchi),
                lduAddr().patchSortStartAddr(patchi),
                #else
                lduAddr().patchAddr(patchi),
                #endif
                boundaryCoeffs(patchi).component(0),
                H1_
            );
        }
    }

    H1_.primitiveFieldRef() /= psi_.mesh().V();
    H1_.correctBoundaryConditions();

    return tH1;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvMatrix<Type>::
flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorInFunction
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary"
               " of fvSchemes."
            << abort(FatalError);
    }

    if (nMatrices() > 1)
    {
        FatalErrorInFunction
            << "Flux requested but " << psi_.name()
            << " can't handle multiple fvMatrix."
            << abort(FatalError);
    }

    auto tfieldFlux = GeometricField<Type, fvsPatchField, surfaceMesh>::New
    (
        "flux(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions(),
        lduMatrix::faceH(psi_.primitiveField())
    );
    auto& fieldFlux = tfieldFlux.ref();
    fieldFlux.setOriented();

    auto& ffbf = fieldFlux.boundaryFieldRef();

    const label nPatches = psi_.boundaryField().size();
    label fieldi = 0;

    if (!useImplicit_)
    {
        // Direct computation without intermediate FieldField storage
        forAll(ffbf, patchi)
        {
            // Internal contribution: internalCoeffs * patchInternalField
            ffbf[patchi] =
                cmptMultiply
                (
                    Field<Type>(internalCoeffs(patchi)),
                    psi_.boundaryField()[patchi].patchInternalField()
                );

            // Neighbour contribution: boundaryCoeffs * patchNeighbourField
            if (psi_.boundaryField()[patchi].coupled())
            {
                ffbf[patchi] -=
                    cmptMultiply
                    (
                        Field<Type>(boundaryCoeffs(patchi)),
                        psi_.boundaryField()[patchi].patchNeighbourField()
                    );
            }
            else
            {
                ffbf[patchi] -= Field<Type>(boundaryCoeffs(patchi));
            }
            //DebugVar(gSum(ffbf[patchi]))
        }
    }
    else
    {
        // Implicit case - needs intermediate storage for mapContributions
        FieldField<Field, Type> InternalContrib(nPatches);
        FieldField<Field, Type> fluxInternalContrib(nPatches);
        forAll(fluxInternalContrib, patchi)
        {
            InternalContrib.set(patchi, new Field<Type>(internalCoeffs(patchi)));
            fluxInternalContrib.set(patchi, new Field<Type>(internalCoeffs(patchi)));
        }
        mapContributions(fieldi, fluxInternalContrib, InternalContrib, true);

        FieldField<Field, Type> NeighbourContrib(nPatches);
        FieldField<Field, Type> fluxBoundaryContrib(nPatches);
        forAll(fluxBoundaryContrib, patchi)
        {
            NeighbourContrib.set(patchi, new Field<Type>(boundaryCoeffs(patchi)));
            fluxBoundaryContrib.set(patchi, new Field<Type>(boundaryCoeffs(patchi)));
        }
        mapContributions(fieldi, fluxBoundaryContrib, NeighbourContrib, false);

        forAll(ffbf, patchi)
        {
            ffbf[patchi] = InternalContrib[patchi] - NeighbourContrib[patchi];
        }
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


template<class Type>
const Foam::dictionary& Foam::fvMatrix<Type>::solverDict
(
    const word& name
) const
{
    return psi_.mesh().solverDict(name);
}


template<class Type>
const Foam::dictionary& Foam::fvMatrix<Type>::solverDict() const
{
    return psi_.mesh().solverDict
    (
        psi_.select(psi_.mesh().data().isFinalIteration())
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::operator=(const fvMatrix<Type>& fvmv)
{
    if (this == &fvmv)
    {
        return;  // Self-assignment is a no-op
    }

    if (&psi_ != &(fvmv.psi_))
    {
        FatalErrorInFunction
            << "different fields"
            << abort(FatalError);
    }

    dimensions_ = fvmv.dimensions_;
    lduMatrix::operator=(fvmv);
    source_ = fvmv.source_;
    internalCoeffs_ = fvmv.internalCoeffs_;
    boundaryCoeffs_ = fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            *fvmv.faceFluxCorrectionPtr_
        );
    }

    useImplicit_ = fvmv.useImplicit_;
    lduAssemblyName_ = fvmv.lduAssemblyName_;
}


template<class Type>
void Foam::fvMatrix<Type>::operator=(const tmp<fvMatrix<Type>>& tfvmv)
{
    operator=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ += fvmv.dimensions_;
    lduMatrix::operator+=(fvmv);
    source_ += fvmv.source_;
    internalCoeffs_ += fvmv.internalCoeffs_;
    boundaryCoeffs_ += fvmv.boundaryCoeffs_;

    useImplicit_ = fvmv.useImplicit_;
    lduAssemblyName_ = fvmv.lduAssemblyName_;
    nMatrix_ = fvmv.nMatrix_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            *fvmv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=(const tmp<fvMatrix<Type>>& tfvmv)
{
    operator+=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "-=");

    dimensions_ -= fvmv.dimensions_;
    lduMatrix::operator-=(fvmv);
    source_ -= fvmv.source_;
    internalCoeffs_ -= fvmv.internalCoeffs_;
    boundaryCoeffs_ -= fvmv.boundaryCoeffs_;

    useImplicit_ = fvmv.useImplicit_;
    lduAssemblyName_ = fvmv.lduAssemblyName_;
    nMatrix_ = fvmv.nMatrix_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = std::make_unique<faceFluxFieldType>
        (
            -*fvmv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=(const tmp<fvMatrix<Type>>& tfvmv)
{
    operator-=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().V()*su.field();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().V()*su.field();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= psi().mesh().V()*su;
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += psi().mesh().V()*su;
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const volScalarField::Internal& dsf
)
{
    dimensions_ *= dsf.dimensions();
    lduMatrix::operator*=(dsf.field());
    source_ *= dsf.field();

    forAll(psi_.boundaryField(), patchi)
    {
        scalarField pisf
        (
            dsf.mesh().boundary()[patchi].patchInternalField(dsf.field())
        );

        internalCoeffs(patchi) *= pisf;
        boundaryCoeffs(patchi) *= pisf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorInFunction
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const tmp<volScalarField::Internal>& tfld
)
{
    operator*=(tfld());
    tfld.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const tmp<volScalarField>& tfld
)
{
    operator*=(tfld());
    tfld.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& mat1,
    const fvMatrix<Type>& mat2,
    const char* op
)
{
    if (&mat1.psi() != &mat2.psi())
    {
        FatalErrorInFunction
            << "Incompatible fields for operation\n    "
            << "[" << mat1.psi().name() << "] "
            << op
            << " [" << mat2.psi().name() << "]"
            << abort(FatalError);
    }

    if
    (
        dimensionSet::checking()
     && mat1.dimensions() != mat2.dimensions()
    )
    {
        FatalErrorInFunction
            << "Incompatible dimensions for operation\n    "
            << "[" << mat1.psi().name() << mat1.dimensions()/dimVolume << " ] "
            << op
            << " [" << mat2.psi().name() << mat2.dimensions()/dimVolume << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& mat,
    const DimensionedField<Type, volMesh>& fld,
    const char* op
)
{
    if
    (
        dimensionSet::checking()
     && mat.dimensions()/dimVolume != fld.dimensions()
    )
    {
        FatalErrorInFunction
            << "Incompatible dimensions for operation\n    "
            << "[" << mat.psi().name() << mat.dimensions()/dimVolume << " ] "
            << op
            << " [" << fld.name() << fld.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& mat,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if
    (
        dimensionSet::checking()
     && mat.dimensions()/dimVolume != dt.dimensions()
    )
    {
        FatalErrorInFunction
            << "Incompatible dimensions for operation\n    "
            << "[" << mat.psi().name() << mat.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    fvMatrix<Type>& mat,
    const dictionary& solverControls
)
{
    return mat.solve(solverControls);
}

template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    const tmp<fvMatrix<Type>>& tmat,
    const dictionary& solverControls
)
{
    SolverPerformance<Type> solverPerf(tmat.constCast().solve(solverControls));

    tmat.clear();

    return solverPerf;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    fvMatrix<Type>& mat,
    const word& name
)
{
    return mat.solve(name);
}

template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    const tmp<fvMatrix<Type>>& tmat,
    const word& name
)
{
    SolverPerformance<Type> solverPerf(tmat.constCast().solve(name));

    tmat.clear();

    return solverPerf;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve(fvMatrix<Type>& mat)
{
    return mat.solve();
}

template<class Type>
Foam::SolverPerformance<Type> Foam::solve(const tmp<fvMatrix<Type>>& tmat)
{
    SolverPerformance<Type> solverPerf(tmat.constCast().solve());

    tmat.clear();

    return solverPerf;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::correction
(
    const fvMatrix<Type>& A
)
{
    tmp<Foam::fvMatrix<Type>> tAcorr = A - (A & A.psi());

    // Delete the faceFluxCorrection from the correction matrix
    // as it does not have a clear meaning or purpose
    tAcorr.ref().faceFluxCorrectionPtr(nullptr);

    return tAcorr;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::correction
(
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<Foam::fvMatrix<Type>> tAcorr = tA - (tA() & tA().psi());

    // Delete the faceFluxCorrection from the correction matrix
    // as it does not have a clear meaning or purpose
    tAcorr.ref().faceFluxCorrectionPtr(nullptr);

    return tAcorr;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "==");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += A.psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tC().psi().mesh().V()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const fvMatrix<Type>& A,
    const Foam::zero
)
{
    return A;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator==
(
    const tmp<fvMatrix<Type>>& tA,
    const Foam::zero
)
{
    return tA;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A
)
{
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<fvMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<fvMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<fvMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<DimensionedField<Type, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().V()*su.field();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().field();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().V()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const tmp<fvMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref().negate();
    tC.ref().source() -= su.value()*A.psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const volScalarField::Internal& dsf,
    const fvMatrix<Type>& A
)
{
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField::Internal>& tdsf,
    const fvMatrix<Type>& A
)
{
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const fvMatrix<Type>& A
)
{
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const volScalarField::Internal& dsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField::Internal>& tdsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const fvMatrix<Type>& A
)
{
    auto tC = tmp<fvMatrix<Type>>::New(A);
    tC.ref() *= ds;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const tmp<fvMatrix<Type>>& tA
)
{
    tmp<fvMatrix<Type>> tC(tA.ptr());
    tC.ref() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const DimensionedField<Type, volMesh>& psi
)
{
    auto tMphi = GeometricField<Type, fvPatchField, volMesh>::New
    (
        "M&" + psi.name(),
        psi.mesh(),
        M.dimensions()/dimVol,
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& Mphi = tMphi.ref();

    // Loop over field components
    if (M.hasDiag())
    {
        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            scalarField psiCmpt(psi.field().component(cmpt));
            scalarField boundaryDiagCmpt(M.diag());
            M.addBoundaryDiag(boundaryDiagCmpt, cmpt);
            Mphi.primitiveFieldRef().replace(cmpt, -boundaryDiagCmpt*psiCmpt);
        }
    }
    else
    {
        Mphi.primitiveFieldRef() = Zero;
    }

    Mphi.primitiveFieldRef() += M.lduMatrix::H(psi.field()) + M.source();
    M.addBoundarySource(Mphi.primitiveFieldRef());

    Mphi.primitiveFieldRef() /= -psi.mesh().V();
    Mphi.correctBoundaryConditions();

    return tMphi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<DimensionedField<Type, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & psi;
    tM.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const tmp<DimensionedField<Type, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::operator&
(
    const tmp<fvMatrix<Type>>& tM,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvMatrix<Type>& fvm)
{
    os  << static_cast<const lduMatrix&>(fvm) << nl
        << fvm.dimensions_ << nl
        << fvm.source_ << nl
        << fvm.internalCoeffs_ << nl
        << fvm.boundaryCoeffs_ << endl;

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * Expression Templates  * * * * * * * * * * * * * //

template<class Type>
template<typename E>
Foam::fvMatrix<Type>::fvMatrix
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const Expression::fvMatrixExpression
    <
        E,
        typename E::DiagExpr,
        typename E::UpperExpr,
        typename E::LowerExpr,
        typename E::FaceFluxExpr,
        typename E::SourceExpr
   >& expr
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    useImplicit_(false),
    lduAssemblyName_(),
    nMatrix_(0),
    dimensions_(expr.dimensions()),
    source_(psi.size(), Zero),
    internalCoeffs_(psi.mesh().nBoundaryFaces(), Zero),
    boundaryCoeffs_(psi.mesh().nBoundaryFaces(), Zero),
    coeffsOffset_(psi.mesh().boundary().size() + 1, 0)
{
    DebugInFunction
        << "Constructing fvMatrix<Type> from expression for field "
        << psi_.name() << endl;

    checkImplicit();

    // Populate offset array for per-patch coefficient access
    forAll(psi.mesh().boundary(), patchi)
    {
        const label patchSize = psi.mesh().boundary()[patchi].size();
        coeffsOffset_[patchi+1] = coeffsOffset_[patchi] + patchSize;
    }

    // Fill in diag,upper,lower etc.
    expr.evaluate(*this);

    auto& psiRef = this->psi(0);
    const label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryFieldRef().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
Foam::Expression::fvMatrixConstRefWrap<Foam::fvMatrix<Type>>
Foam::fvMatrix<Type>::expr() const
{
    return Expression::fvMatrixConstRefWrap<Foam::fvMatrix<Type>>(*this);
}


template<class Type>
template<typename E>
void Foam::fvMatrix<Type>::operator=
(
    const Expression::fvMatrixExpression
    <
        E,
        typename E::DiagExpr,
        typename E::UpperExpr,
        typename E::LowerExpr,
        typename E::FaceFluxExpr,
        typename E::SourceExpr
   >& expr
)
{
    expr.evaluate(*this);
}


// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "fvMatrixSolve.C"

// ************************************************************************* //
