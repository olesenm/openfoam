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
    Functors for GPU offloading of GAMG matrix agglomeration.
\*---------------------------------------------------------------------------*/

#ifndef GAMGAgglomerationF_H
#define GAMGAgglomerationF_H

namespace Foam
{
namespace GAMG
{

template<class Type, bool nonNegative>
struct restrictT
{
    Field<Type>& cf_;

    const Field<Type>& ff_;

    const labelList& sort_;
    const labelList& target_;
    const labelList& targetStart_;

    const Type zero_;

    restrictT
    (
        Field<Type>& cf,
        const Field<Type>& ff,
        const labelList& sort,
        const labelList& target,
        const labelList& targetStart
    )
    :
        cf_(cf),
        ff_(ff),
        sort_(sort),
        target_(target),
        targetStart_(targetStart),
        zero_(pTraits<Type>::zero)
    {}

    void operator()(const label index) const
    {
        if (nonNegative && target_[index] < 0)
        {
            return;
        }

        const label coarseCell = target_[index];
        const label start = targetStart_[index];
        const label end = targetStart_[index + 1];

        Type result(zero_);

        for (label i = start; i < end; ++i)
        {
            result += ff_[sort_[i]];
        }

        cf_[coarseCell] = result;
    }
};


} // End namespace GAMG
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
