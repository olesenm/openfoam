/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "fvcReconstruct.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"

#ifdef FOAM_OFFLOAD
#include "lduAddressingFunctors.H"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> reconstructMag(const surfaceScalarField& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    tmp<volScalarField> treconField
    (
        new volScalarField
        (
            IOobject
            (
                "reconstruct("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(ssf.dimensions()/dimArea, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    scalarField& rf = treconField.ref();
    #ifdef FOAM_OFFLOAD
    const labelUList& ownStart = mesh.lduAddr().ownerStartAddr();
    const labelUList& neiStart = mesh.lduAddr().losortStartAddr();
    const labelUList& losort = mesh.lduAddr().losortAddr();

    auto ownFun = [&Sf, &Cf, &C, &ssf, &magSf]
        (const label celli, const label facei)
    {
        return
        (
            (Sf[facei] & (Cf[facei] - C[celli]))
        *ssf[facei]/magSf[facei]
        );
    };

    auto neiFun =[&Sf, &Cf, &C, &ssf, &magSf]
        (const label celli, const label facei)
    {
        return
        (
        -((Sf[facei] & (Cf[facei] - C[celli]))
            *ssf[facei]/magSf[facei])
        );
    };

    auto plusOp = [](scalar& result, const scalar value)
    {
        result += value;
    };

    labelRange cells(rf.size());

    apply_policy(cells.size(), [&](auto policy)
    {
        std::transform
        (
            policy,
            cells.begin(),
            cells.end(),
            rf.begin(),
            cellReduceFunctor
            <
                scalar,
                decltype(ownFun),
                decltype(neiFun),
                decltype(plusOp),
                decltype(plusOp)
            >
            (
                pTraits<scalar>::zero,
                ownStart,
                neiStart,
                losort,
                ownFun,
                neiFun,
                plusOp,
                plusOp
            )
        );
    });
    #else
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        rf[own] += (Sf[facei] & (Cf[facei] - C[own]))*ssf[facei]/magSf[facei];
        rf[nei] -= (Sf[facei] & (Cf[facei] - C[nei]))*ssf[facei]/magSf[facei];
    }
    #endif

    const surfaceScalarField::Boundary& bsf = ssf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvsPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelUList& patchCells = mesh.lduAddr().patchSortCells(patchi);
        const labelUList& sort = mesh.lduAddr().patchSortAddr(patchi);
        const labelUList& sortStart = mesh.lduAddr().patchSortStartAddr(patchi);

        // patch-local views (avoid repeated boundaryField() inside kernel)
        const auto& pssf = ssf.boundaryField()[patchi];

        labelRange indices(pOwner.size());

        // device-callable lambda
        auto neiFun = [&](const label id, const label facei) -> scalar
        {
            const label celli = pOwner[facei];

            return (pSf[facei] & (pCf[facei] - C[celli]))
                *pssf[facei]/pMagSf[facei];
        };

        apply_policy(pOwner.size(), [&]( auto policy )
        {
            std::for_each
            (
                policy,
                indices.begin(),
                indices.end(),
                patchCellReduceFunctor<scalar, decltype(neiFun),
                    plusEqOp<scalar>>
                (
                    rf,
                    sortStart,
                    sort,
                    patchCells,
                    neiFun,
                    plusEqOp<scalar>()
                )
            );
        });
        #else
        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];
            rf[own] +=
                (pSf[pFacei] & (pCf[pFacei] - C[own]))
               *psf[pFacei]/pMagSf[pFacei];
        }
        #endif
    }

    rf /= mesh.V();

    treconField.ref().correctBoundaryConditions();

    return treconField;
}


tmp<volScalarField> reconstructMag(const tmp<surfaceScalarField>& tssf)
{
    tmp<volScalarField> tvf
    (
        fvc::reconstructMag(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
