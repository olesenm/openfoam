/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "gaussGrad.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GradFieldType> tgGrad
    (
        new GradFieldType
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>(ssf.dimensions()/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GradFieldType& gGrad = tgGrad.ref();

    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    #ifdef FOAM_OFFLOAD
    auto fun = [&] (const label, const label facei) -> GradType
    {
        return Sf[facei]*issf[facei];
    };

    auto ownOp = [](GradType& result, const GradType& value)
    {
        result += value;
    };

    auto neiOp = [](GradType& result, const GradType& value)
    {
        result -= value;
    };

    labelRange cellRange(igGrad.size());

    apply_policy(igGrad.size(), [&](auto policy)
    {
        std::transform
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            igGrad.begin(),   // d = old igGrad[celli]
            cellReduceFunctor
            <
                GradType,
                decltype(fun),
                decltype(fun),
                decltype(ownOp),
                decltype(neiOp)
            >
            (
                pTraits<GradType>::zero,
                mesh.lduAddr().ownerStartAddr(),
                mesh.lduAddr().losortStartAddr(),
                mesh.lduAddr().losortAddr(),
                fun,
                fun,
                ownOp,
                neiOp
            )
        );
    });
    #else
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }
    #endif

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelUList& patchCells = mesh.lduAddr().patchSortCells(patchi);
        const labelUList& sort = mesh.lduAddr().patchSortAddr(patchi);
        const labelUList& sortStart = mesh.lduAddr().patchSortStartAddr(patchi);

        const labelRange cellRange(patchCells.size());

        auto neiFun = [&](const label, const label facei) -> GradType
        {
            return pSf[facei]*pssf[facei];
        };

        auto neiOp = [](GradType& result, const GradType& value)
        {
            result += value;
        };

        apply_policy(cellRange.size(), [&](auto policy)
        {
            std::for_each
            (
                policy,
                cellRange.begin(),
                cellRange.end(),
                patchCellReduceFunctor
                <
                    GradType,
                    decltype(neiFun),
                    decltype(neiOp)
                >
                (
                    igGrad,
                    sortStart,
                    sort,
                    patchCells,
                    neiFun,
                    neiOp
                )
            );
        });
        #else
        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
        #endif
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    tmp<GradFieldType> tgGrad
    (
        gradf(tinterpScheme_().interpolate(vsf), name)
    );
    GradFieldType& gGrad = tgGrad.ref();

    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void Foam::fv::gaussGrad<Type>::calcGrad
(
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& gGrad,
    const GeometricField<Type, fvPatchField, volMesh>& vsf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    DebugPout<< "gaussGrad<Type>::calcGrad on " << vsf.name()
        << " into " << gGrad.name() << endl;

    const fvMesh& mesh = vsf.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    const auto tssf(tinterpScheme_().interpolate(vsf));
    const auto& ssf = tssf();

    gGrad = dimensioned<GradType>(vsf.dimensions()/dimLength, Zero);

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    forAll(owner, facei)
    {
        const GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    correctBoundaryConditions(vsf, gGrad);
}


template<class Type>
void Foam::fv::gaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    auto& gGradbf = gGrad.boundaryFieldRef();

    forAll(vsf.boundaryField(), patchi)
    {
        if (!vsf.boundaryField()[patchi].coupled())
        {
            const vectorField n
            (
                vsf.mesh().Sf().boundaryField()[patchi]
              / vsf.mesh().magSf().boundaryField()[patchi]
            );

            gGradbf[patchi] += n *
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
     }
}


// ************************************************************************* //
