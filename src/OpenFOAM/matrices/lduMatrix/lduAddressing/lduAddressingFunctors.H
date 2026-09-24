/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    functors for LDU-level offloading

\*---------------------------------------------------------------------------*/

#ifdef FOAM_OFFLOAD

#ifndef lduAddressingFunctors_H
#define lduAddressingFunctors_H

#include "lduAddressing.H"
#include "ops.H"

namespace Foam
{
    template<class Type, class OwnFun, class NeiFun, class OwnOp, class NeiOp>
    struct lduAddressingFunctor
    {
        const labelUList& ownStart_;
        const labelUList& neiStart_;
        const labelUList& losort_;

        OwnFun ownFun_;
        NeiFun neiFun_;
        OwnOp ownOp_;
        NeiOp neiOp_;

        lduAddressingFunctor
        (
            const labelUList& ownStart,
            const labelUList& neiStart,
            const labelUList& losort,
            const OwnFun ownFun,
            const NeiFun neiFun,
            const OwnOp ownOp,
            const NeiOp neiOp
        )
        :
            ownStart_(ownStart),
            neiStart_(neiStart),
            losort_(losort),
            ownFun_(ownFun),
            neiFun_(neiFun),
            ownOp_(ownOp),
            neiOp_(neiOp)
        {}

        Type operator()(const label& id,const Type& d)
        {
            Type out = d;
            label oStart = ownStart_[id];
            label oSize = ownStart_[id+1] - oStart;
    
            label nStart = neiStart_[id];
            label nSize = neiStart_[id+1] - nStart;
    
            for(label i = 0; i<oSize; i++)
            {
                label face = oStart + i;
                out = ownOp_(out, ownFun_(id, face));
            }
    
            for(label i = 0; i<nSize; i++)
            {
                label face = losort_[nStart + i];
                out = neiOp_(out, neiFun_(id, face));
            }
            return out;
        }
    };


    template<class Type, class LUType, class Op>
    struct matrixCoeffsMultiplyFunctor
    {
        const Field<Type>& psi_;
        const scalarField& coeffs_;
        const labelUList& addr_;
        Op op_;

        matrixCoeffsMultiplyFunctor
        (
            const Field<Type>& psi,
            const scalarField& coeffs,
            const labelUList& addr,
            const Op& op
        )
        :
            psi_(psi),
            coeffs_(coeffs),
            addr_(addr),
            op_(op)
        {}

        Type operator()(const label& cell, const label& face)
        {
            return op_(coeffs_[face]*psi_[addr_[face]]);
        }
    };


    template<class Type, class Op>
    struct matrixCoeffsFunctor
    {
        const scalarField& coeffs_;
        Op op_;

        matrixCoeffsFunctor
        (
            const scalarField& coeffs,
            const Op& op
        )
        :
            coeffs_(coeffs),
            op_(op)
        {}

        Type operator()(const label& cell, const label& face)
        {
            return op_(coeffs_[face]);
        }
    };


    template<class Type, class Input, class OwnFun, class NeiFun, class OwnOp, 
        class NeiOp>
    inline void matrixOperation
    (
        Input in,
        List<Type>& out,
        const lduAddressing& addr,
        OwnFun o,
        NeiFun n,
        OwnOp oo,
        NeiOp no
    )
    {
        const labelRange indices(out.size());
        apply_policy(indices.size(), [&](auto policy)
        {
            std::transform
            (
                policy,
                indices.begin(),
                indices.end(),
                in,
                out.begin(),
                lduAddressingFunctor<Type, OwnFun, NeiFun, OwnOp, NeiOp>
                (
                    addr.ownerStartAddr(),
                    addr.losortStartAddr(),
                    addr.losortAddr(),
                    o,
                    n,
                    oo,
                    no
                )
            );
        });
    }


    template<class Type, class Input, class OwnFun, class NeiFun>
    inline void matrixOperation
    (
        Input in,
        List<Type>& out,
        const lduAddressing& addr,
        OwnFun o,
        NeiFun n
    )
    {
        matrixOperation
        (
            in,
            out,
            addr,
            o,
            n,
            sumOp<Type>(),
            sumOp<Type>()
        );
    }
}

#endif

#endif  // FOAM_OFFLOAD

// ************************************************************************* //
