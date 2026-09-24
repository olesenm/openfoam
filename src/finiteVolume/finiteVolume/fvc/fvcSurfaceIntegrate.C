/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

template<class Type>
void surfaceIntegrate
(
    Field<Type>& ivf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const Field<Type>& issf = ssf;

    #ifdef FOAM_OFFLOAD
    const lduAddressing& addr = mesh.lduAddr();

    const labelUList& ownStart = addr.ownerStartAddr();
    const labelUList& neiStart = addr.losortStartAddr();
    const labelUList& losort = addr.losortAddr();

    auto faceOp = [&](const label, const label facei)
    {
        return issf[facei];
    };

    labelRange cellRange(0, ivf.size());

    apply_policy(cellRange.size(), [&](auto policy)
    {
        std::transform
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            ivf.begin(),
            cellReduceFunctor
            <
                Type,
                decltype(faceOp),
                decltype(faceOp),
                plusEqOp<Type>,
                minusEqOp<Type>
            >
            (
                pTraits<Type>::zero,
                ownStart,
                neiStart,
                losort,
                faceOp,
                faceOp,
                plusEqOp<Type>(),
                minusEqOp<Type>()
            )
        );
    });
    #else
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf[facei];
        ivf[neighbour[facei]] -= issf[facei];
    }
    #endif

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelUList& patchCells = mesh.lduAddr().patchSortCells(patchi);
        const labelUList& neiStart = mesh.lduAddr().patchSortStartAddr(patchi);
        const labelUList& losort = mesh.lduAddr().patchSortAddr(patchi);

        labelRange indices(patchCells.size());

        auto neiFun = [&pssf](const label, const label facei)
        {
            return pssf[facei];
        };

        auto neiOp = [](Type& result, const Type& value)
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
                patchCellReduceFunctor<Type, decltype(neiFun),
                    decltype(neiOp)>
                (
                    ivf,
                    neiStart,
                    losort,
                    patchCells,
                    neiFun,
                    neiOp
                )
            );
        });
        #else
        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
        #endif
    }

    ivf /= mesh.Vsc();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions()/dimVol, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    surfaceIntegrate(vf.primitiveFieldRef(), ssf);
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        fvc::surfaceIntegrate(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions(), Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    #ifdef FOAM_OFFLOAD
    const lduAddressing& addr = mesh.lduAddr();

    const labelUList& ownStart = addr.ownerStartAddr();
    const labelUList& neiStart = addr.losortStartAddr();
    const labelUList& losort = addr.losortAddr();

    auto faceOp = [&](const label, const label facei)
    {
        return ssf[facei];
    };

    labelRange cellRange(0, vf.size());

    apply_policy(cellRange.size(), [&](auto policy)
    {
        std::transform
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            vf.begin(),
            cellReduceFunctor
            <
                Type,
                decltype(faceOp),
                decltype(faceOp),
                plusEqOp<Type>,
                plusEqOp<Type>
            >
            (
                pTraits<Type>::zero,
                ownStart,
                neiStart,
                losort,
                faceOp,
                faceOp,
                plusEqOp<Type>(),
                plusEqOp<Type>()
            )
        );
    });
    #else
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += ssf[facei];
        vf[neighbour[facei]] += ssf[facei];
    }
    #endif

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelUList& patchCells =
            mesh.lduAddr().patchSortCells(patchi);

        const labelUList& neiStart =
            mesh.lduAddr().patchSortStartAddr(patchi);

        const labelUList& losort =
            mesh.lduAddr().patchSortAddr(patchi);

        auto neiFun = [&pssf](const label, const label facei)
        {
            return pssf[facei];
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
                    vf,
                    neiStart,
                    losort,
                    patchCells,
                    neiFun,
                    neiOp
                )
            );
        });
        #else
        forAll(mesh.boundary()[patchi], facei)
        {
            vf[pFaceCells[facei]] += pssf[facei];
        }
        #endif
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> surfaceSum
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf = surfaceSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
