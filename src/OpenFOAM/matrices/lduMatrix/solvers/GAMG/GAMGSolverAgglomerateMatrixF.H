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

Description
    Functors for GPU offloading of GAMG matrix agglomeration.

\*---------------------------------------------------------------------------*/

#ifndef GAMGSolverAgglomerateMatrixF_H
#define GAMGSolverAgglomerateMatrixF_H

namespace Foam
{
namespace GAMG
{

struct asymAgglomerate
{
    scalar* coarseUpper_;
    scalar* coarseLower_;
    const scalar* uf_;
    const scalar* lf_;
    const bool* flip_;
    const label* sort_;
    const label* target_;
    const label* targetStart_;

    asymAgglomerate
    (
        scalarField& coarseUpper,
        scalarField& coarseLower,
        const scalarField& uf,
        const scalarField& lf,
        const boolList& flip,
        const labelList& sort,
        const labelList& target,
        const labelList& targetStart
    )
    :
        coarseUpper_(coarseUpper.begin()),
        coarseLower_(coarseLower.begin()),
        uf_(uf.begin()),
        lf_(lf.begin()),
        flip_(flip.begin()),
        sort_(sort.begin()),
        target_(target.begin()),
        targetStart_(targetStart.begin())
    {}

    void operator()(const label index) const
    {
        const label coarseFace = target_[index];

        scalar uc = coarseUpper_[coarseFace];
        scalar lc = coarseLower_[coarseFace];

        const label start = targetStart_[index];
        const label end = targetStart_[index + 1];

        for (label i = start; i < end; ++i)
        {
            const label idx = sort_[i];
            const bool flipVal = flip_[idx];

            const scalar ufVal = uf_[idx];
            const scalar lfVal = lf_[idx];

            if (!flipVal)
            {
                uc += ufVal;
                lc += lfVal;
            }
            else
            {
                uc += lfVal;
                lc += ufVal;
            }
        }

        coarseUpper_[coarseFace] = uc;
        coarseLower_[coarseFace] = lc;
    }
};


struct diagAsymAgglomerate
{
    scalar* coarseDiag_;
    const scalar* uf_;
    const scalar* lf_;
    const label* sort_;
    const label* target_;
    const label* targetStart_;

    diagAsymAgglomerate
    (
        scalarField& coarseDiag,
        const scalarField& uf,
        const scalarField& lf,
        const labelList& sort,
        const labelList& target,
        const labelList& targetStart
    )
    :
        coarseDiag_(coarseDiag.begin()),
        uf_(uf.begin()),
        lf_(lf.begin()),
        sort_(sort.begin()),
        target_(target.begin()),
        targetStart_(targetStart.begin())
    {}

    void operator()(const label index) const
    {
        const label targetFace = target_[index];

        if (targetFace >= 0)
        {
            return;
        }

        const label id = -1 - targetFace;

        scalar out = coarseDiag_[id];

        const label start = targetStart_[index];
        const label end = targetStart_[index + 1];

        for (label i = start; i < end; ++i)
        {
            const label idx = sort_[i];

            out += uf_[idx] + lf_[idx];
        }

        coarseDiag_[id] = out;
    }
};

struct symAgglomerate
{
    scalar* coarseUpper_;
    const scalar* ff_;
    const label* sort_;
    const label* target_;
    const label* targetStart_;

    symAgglomerate
    (
        scalarField& coarseUpper,
        const scalarField& ff,
        const labelList& sort,
        const labelList& target,
        const labelList& targetStart
    )
    :
        coarseUpper_(coarseUpper.begin()),
        ff_(ff.begin()),
        sort_(sort.begin()),
        target_(target.begin()),
        targetStart_(targetStart.begin())
    {}

    void operator()(const label index) const
    {
        const label coarseFace = target_[index];

        if (coarseFace < 0)
        {
            return;
        }

        scalar out = coarseUpper_[coarseFace];

        const label start = targetStart_[index];
        const label end = targetStart_[index + 1];

        for (label i = start; i < end; ++i)
        {
            out += ff_[sort_[i]];
        }

        coarseUpper_[coarseFace] = out;
    }
};


struct diagSymAgglomerate
{
    scalar* coarseDiag_;
    const scalar* ff_;
    const label* sort_;
    const label* target_;
    const label* targetStart_;

    diagSymAgglomerate
    (
        scalarField& coarseDiag,
        const scalarField& ff,
        const labelList& sort,
        const labelList& target,
        const labelList& targetStart
    )
    :
        coarseDiag_(coarseDiag.begin()),
        ff_(ff.begin()),
        sort_(sort.begin()),
        target_(target.begin()),
        targetStart_(targetStart.begin())
    {}

    void operator()(const label index) const
    {
        const label targetFace = target_[index];

        if (targetFace >= 0)
        {
            return;
        }

        const label id = -1 - targetFace;

        scalar out = coarseDiag_[id];

        const label start = targetStart_[index];
        const label end = targetStart_[index + 1];

        for (label i = start; i < end; ++i)
        {
            out += 2*ff_[sort_[i]];
        }

        coarseDiag_[id] = out;
    }
};


} // End namespace GAMG
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //