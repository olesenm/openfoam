/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "cellLimitedGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::limitGradient
(
    const Field<scalar>& limiter,
    Field<vector>& gIf
) const
{
    gIf *= limiter;
}


template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::limitGradient
(
    const Field<vector>& limiter,
    Field<tensor>& gIf
) const
{
    labelRange cellRange(gIf.size());

    apply_policy(gIf.size(), [&](auto policy)
    {   
        std::for_each
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            [&](label celli)
            {
                gIf[celli] = tensor
                (
                    cmptMultiply(limiter[celli], gIf[celli].x()),
                    cmptMultiply(limiter[celli], gIf[celli].y()),
                    cmptMultiply(limiter[celli], gIf[celli].z())
                );
            }
        );
    });
}


template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::calcGrad
(
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& g,
    const GeometricField<Type, fvPatchField, volMesh>& vsf
) const
{
    const fvMesh& mesh = vsf.mesh();

    basicGradScheme_().calcGrad(g, vsf);

    if (k_ < SMALL)
    {
        return;
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    Field<Type> maxVsf(vsf.primitiveField());
    Field<Type> minVsf(vsf.primitiveField());

    #ifdef FOAM_OFFLOAD
    const labelUList& ownStart = mesh.lduAddr().ownerStartAddr();
    const labelUList& neiStart = mesh.lduAddr().losortStartAddr();
    const labelUList& losort = mesh.lduAddr().losortAddr();

    labelRange ownerRange(vsf.size());

    apply_policy(ownerRange.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            ownerRange.begin(),
            ownerRange.end(),
            [&](const label celli)
            {
                Type maxVal(vsf[celli]);
                Type minVal(vsf[celli]);

                for
                (
                    label facei = ownStart[celli];
                    facei < ownStart[celli + 1];
                    ++facei
                )
                {
                    const label nei = neighbour[facei];

                    maxVal = max(maxVal, vsf[nei]);
                    minVal = min(minVal, vsf[nei]);
                }

                for
                (
                    label i = neiStart[celli];
                    i < neiStart[celli + 1];
                    ++i
                )
                {
                    const label facei = losort[i];
                    const label own = owner[facei];

                    maxVal = max(maxVal, vsf[own]);
                    minVal = min(minVal, vsf[own]);
                }

                maxVsf[celli] = maxVal;
                minVsf[celli] = minVal;
            }
        );
    });
    #else
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const Type& vsfOwn = vsf[own];
        const Type& vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }
    #endif

    const auto& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchField<Type>& psf = bsf[patchi];
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const Field<Type> psfNei(psf.patchNeighbourField());

            labelRange pOwnerRange(pOwner.size());

            apply_policy(pOwnerRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    pOwnerRange.begin(),
                    pOwnerRange.end(),
                    [&](label pFacei)
                    {
                        const label own = pOwner[pFacei];
                        const Type& vsfNei = psfNei[pFacei];

                        maxVsf[own] = max(maxVsf[own], vsfNei);
                        minVsf[own] = min(minVsf[own], vsfNei);
                    }
                );
            });
        }
        else
        {
            labelRange pOwnerRange(pOwner.size());

            apply_policy(pOwnerRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    pOwnerRange.begin(),
                    pOwnerRange.end(),
                    [&](label pFacei)
                    {
                        const label own = pOwner[pFacei];
                        const Type& vsfNei = psf[pFacei];

                        maxVsf[own] = max(maxVsf[own], vsfNei);
                        minVsf[own] = min(minVsf[own], vsfNei);
                    }
                );
            });
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const Field<Type> maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;
    }


    // Create limiter initialized to 1
    // Note: the limiter is not permitted to be > 1
    Field<Type> limiter(vsf.primitiveField().size(), pTraits<Type>::one);

    #ifdef FOAM_OFFLOAD
    labelRange cellRange(limiter.size());

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        apply_policy(cellRange.size(), [&](auto policy)
        {
            std::for_each
            (
                policy,
                cellRange.begin(),
                cellRange.end(),
                [=, &limiter, &maxVsf, &minVsf, &Cf, &C, &g,
                &ownStart, &neiStart, &losort]
                (const label celli)
                {
                    Type lim(limiter[celli]);

                    for
                    (
                        label facei = ownStart[celli];
                        facei < ownStart[celli + 1];
                        ++facei
                    )
                    {
                        limitFaceCmpt
                        (
                            setComponent(lim, cmpt),
                            component(maxVsf[celli], cmpt),
                            component(minVsf[celli], cmpt),
                            component((Cf[facei] - C[celli]) & g[celli], cmpt)
                        );
                    }

                    for
                    (
                        label i = neiStart[celli];
                        i < neiStart[celli + 1];
                        ++i
                    )
                    {
                        const label facei = losort[i];

                        limitFaceCmpt
                        (
                            setComponent(lim, cmpt),
                            component(maxVsf[celli], cmpt),
                            component(minVsf[celli], cmpt),
                            component((Cf[facei] - C[celli]) & g[celli], cmpt)
                        );
                    }

                    limiter[celli] = lim;
                }
            );
        });
    }
    #else
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }
    #endif

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD

        const labelUList& patchCells = mesh.lduAddr().patchSortCells(patchi);
        const labelUList& sort = mesh.lduAddr().patchSortAddr(patchi);
        const labelUList& sortStart = mesh.lduAddr().patchSortStartAddr(patchi);

        labelRange patchCellsRange(patchCells.size());

        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
        {
            apply_policy(patchCellsRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    patchCellsRange.begin(),
                    patchCellsRange.end(),
                    [=, &limiter, &maxVsf, &minVsf, &pCf, &C, &g,
                    &patchCells, &sort, &sortStart]
                    (const label patchCelli)
                    {
                        const label own = patchCells[patchCelli];

                        Type lim(limiter[own]);

                        for
                        (
                            label i = sortStart[patchCelli];
                            i < sortStart[patchCelli + 1];
                            ++i
                        )
                        {
                            const label pFacei = sort[i];

                            limitFaceCmpt
                            (
                                setComponent(lim, cmpt),
                                component(maxVsf[own], cmpt),
                                component(minVsf[own], cmpt),
                                component((pCf[pFacei] - C[own]) & g[own], cmpt)
                            );
                        }

                        limiter[own] = lim;
                    }
                );
            });
        }
        #else
        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                ((pCf[pFacei] - C[own]) & g[own])
            );
        }
        #endif
    }

    if (fv::debug)
    {
        auto limits = gMinMax(limiter);
        auto avg = gAverage(limiter);

        Info<< "gradient limiter for: " << vsf.name()
            << " min = " << limits.min()
            << " max = " << limits.max()
            << " average: " << avg << endl;
    }

    limitGradient(limiter, g);
    g.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, g);
}


template<class Type, class Limiter>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::cellLimitedGrad<Type, Limiter>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GradFieldType> tGrad
    (
        new GradFieldType
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vsf.dimensions()/dimLength,
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    
    calcGrad(tGrad.ref(), vsf);

    return tGrad;
}


// ************************************************************************* //
