/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2025 OpenCFD Ltd.
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

Class
    Foam::VectorSpaceOps

Description
    Operator functions for VectorSpace.

\*---------------------------------------------------------------------------*/

#ifndef Foam_VectorSpaceOps_H
#define Foam_VectorSpaceOps_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Recursive execution. Terminating at \<N\>, starting at index \<I\>
template<direction N, direction I=0>
struct VectorSpaceOps
{
    //- Like std::copy_n() but with templated loops (constexpr).
    //  \param [in] input indexable input data
    //  \param [out] result indexable output data
    template<class Input, class Output>
    static inline constexpr void copy_n
    (
        [[maybe_unused]] Input input,
        [[maybe_unused]] Output result
    )
    {
        if constexpr (I < N)
        {
            result[I] = input[I];
            VectorSpaceOps<N, I+1>::copy_n(input, result);
        }
    }

    //- Like std::fill_n() but with templated loops (constexpr).
    //  \param [out] result indexable output data
    //  \param val the value to assign for each entry
    template<class Output, class T>
    static inline constexpr void fill_n
    (
        [[maybe_unused]] Output result,
        [[maybe_unused]] const T& val
    )
    {
        if constexpr (I < N)
        {
            result[I] = val;
            VectorSpaceOps<N, I+1>::fill_n(result, val);
        }
    }

    //- Apply binary assignment operation to each vector-space component.
    //  \param [in,out] vs vector-space (indexed) data
    //  \param s scalar/component data (non-indexed)
    template<class V, class S, class EqOp>
    static inline constexpr void eqOpS
    (
        [[maybe_unused]] V& vs,
        [[maybe_unused]] const S& s,
        [[maybe_unused]] EqOp eo
    )
    {
        if constexpr (I < N)
        {
            eo(vs.v_[I], s);
            VectorSpaceOps<N, I+1>::eqOpS(vs, s, eo);
        }
    }

    //- Apply the inplace binary reduction operation.
    //  \param [in,out] s scalar or component data (non-indexed)
    //  \param [in] vs input vector-space (indexed) data
    //  \param eo binary combine/assign operation
    template<class S, class V, class EqOp>
    static inline constexpr void SeqOp
    (
        [[maybe_unused]] S& s,
        [[maybe_unused]] const V& vs,
        [[maybe_unused]] EqOp eo
    )
    {
        if constexpr (I < N)
        {
            eo(s, vs.v_[I]);
            VectorSpaceOps<N, I+1>::SeqOp(s, vs, eo);
        }
    }

    //- Apply the inplace binary assignment operation to the components.
    //  \param [in,out] vs1 vector-space (indexed) data
    //  \param [in] vs2 second vector-space (indexed) data
    //  \param eo binary combine/assign operation
    template<class V1, class V2, class EqOp>
    static inline constexpr void eqOp
    (
        [[maybe_unused]] V1& vs1,
        [[maybe_unused]] const V2& vs2,
        [[maybe_unused]] EqOp eo
    )
    {
        if constexpr (I < N)
        {
            eo(vs1.v_[I], vs2.v_[I]);
            VectorSpaceOps<N, I+1>::eqOp(vs1, vs2, eo);
        }
    }

    //- Apply the binary operation between vector-space and scalar data
    //- and assign the result.
    //  \param [out] vs vector-space (indexed) data
    //  \param [in] vs1 vector-space (indexed) data operand
    //  \param [in] s scalar operand
    //  \param bop binary operation
    template<class V, class V1, class S, class BinaryOp>
    static inline constexpr void opVS
    (
        [[maybe_unused]] V& vs,
        [[maybe_unused]] const V1& vs1,
        [[maybe_unused]] const S& s,
        [[maybe_unused]] BinaryOp bop
    )
    {
        if constexpr (I < N)
        {
            vs.v_[I] = bop(vs1.v_[I], s);
            VectorSpaceOps<N, I+1>::opVS(vs, vs1, s, bop);
        }
    }

    //- Apply the binary operation between scalar and vector-space data
    //- and assign the result.
    //  \param [out] vs vector-space (indexed) data
    //  \param [in] s scalar operand
    //  \param [in] vs1 vector-space (indexed) data operand
    //  \param bop binary operation
    template<class V, class S, class V1, class BinaryOp>
    static inline constexpr void opSV
    (
        [[maybe_unused]] V& vs,
        [[maybe_unused]] const S& s,
        [[maybe_unused]] const V1& vs1,
        [[maybe_unused]] BinaryOp bop
    )
    {
        if constexpr (I < N)
        {
            vs.v_[I] = bop(s, vs1.v_[I]);
            VectorSpaceOps<N, I+1>::opSV(vs, s, vs1, bop);
        }
    }

    //- Apply the binary operation between two vector-space data
    //- and assign the result.
    //  \param [out] vs vector-space (indexed) data
    //  \param [in] vs1 vector-space (indexed) data operand
    //  \param [in] vs2 vector-space (indexed) data operand
    //  \param bop binary operation
    template<class V, class V1, class BinaryOp>
    static inline constexpr void op
    (
        [[maybe_unused]] V& vs,
        [[maybe_unused]] const V1& vs1,
        [[maybe_unused]] const V1& vs2,
        [[maybe_unused]] BinaryOp bop
    )
    {
        if constexpr (I < N)
        {
            vs.v_[I] = bop(vs1.v_[I], vs2.v_[I]);
            VectorSpaceOps<N, I+1>::op(vs, vs1, vs2, bop);
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
