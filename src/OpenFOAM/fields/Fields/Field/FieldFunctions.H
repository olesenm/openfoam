/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#define TEMPLATE template<class Type>
#include "FieldFunctionsM.H"
#include "UPstream.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
inline T deviceZero()
{
    #ifdef FOAM_OFFLOAD

    if constexpr (std::is_arithmetic_v<T>)
    {
        return T(0);
    }
    else
    {
        T t;
        for (direction d = 0; d < pTraits<T>::nComponents; ++d)
        {
            t[d] = typename pTraits<T>::cmptType(0);
        }
        return t;
    }

    #else

    return pTraits<T>::zero;

    #endif
}

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

#ifdef FOAM_OFFLOAD

template<class Type, class OwnFun, class NeiFun, class OwnOp, class NeiOp>
struct cellReduceFunctor
{
    Type zero_;

    const labelUList& ownStart_;
    const labelUList& neiStart_;
    const labelUList& losort_;

    OwnFun ownFun_;
    NeiFun neiFun_;
    OwnOp ownOp_;
    NeiOp neiOp_;

    cellReduceFunctor
    (
        const Type& zero,
        const labelUList& ownStart,
        const labelUList& neiStart,
        const labelUList& losort,
        const OwnFun& ownFun,
        const NeiFun& neiFun,
        const OwnOp& ownOp,
        const NeiOp& neiOp
    )
    :
        zero_(zero),
        ownStart_(ownStart),
        neiStart_(neiStart),
        losort_(losort),
        ownFun_(ownFun),
        neiFun_(neiFun),
        ownOp_(ownOp),
        neiOp_(neiOp)
    {}

    Type operator()(const label celli) const
    {
        Type result(zero_);

        const label ownStart = ownStart_[celli];
        const label ownSize = ownStart_[celli + 1] - ownStart;

        for (label i = 0; i < ownSize; ++i)
        {
            const label facei = ownStart + i;
            ownOp_(result, ownFun_(celli, facei));
        }

        const label neiStart = neiStart_[celli];
        const label neiSize = neiStart_[celli + 1] - neiStart;

        for (label i = 0; i < neiSize; ++i)
        {
            const label facei = losort_[neiStart + i];
            neiOp_(result, neiFun_(celli, facei));
        }

        return result;
    }
};

template<class Type, class NeiFun, class NeiOp>
struct patchCellReduceFunctor
{
    Field<Type>& field_;

    const labelUList& neiStart_;
    const labelUList& losort_;
    const labelUList& patchCells_;

    NeiFun neiFun_;
    NeiOp neiOp_;

    patchCellReduceFunctor
    (
        Field<Type>& field,
        const labelUList& neiStart,
        const labelUList& losort,
        const labelUList& patchCells,
        const NeiFun& neiFun,
        const NeiOp& neiOp
    )
    :
        field_(field),
        neiStart_(neiStart),
        losort_(losort),
        patchCells_(patchCells),
        neiFun_(neiFun),
        neiOp_(neiOp)
    {}

    bool operator()(const label id) const
    {
        const label celli = patchCells_[id];

        Type result(field_[celli]);

        const label neiStart = neiStart_[id];
        const label neiSize = neiStart_[id + 1] - neiStart;

        for (label i = 0; i < neiSize; ++i)
        {
            const label facei = losort_[neiStart + i];
            neiOp_(result, neiFun_(id, facei));
        }

        field_[celli] = result;

        return true;
    }
};
#endif

template<class Type>
void component
(
    Field<typename Field<Type>::cmptType>& res,
    const UList<Type>& f,
    const direction d
);


template<class Type>
void T(Field<Type>& res, const UList<Type>& f);


template<class Type, direction r>
void pow
(
    Field<typename powProduct<Type, r>::type>& res,
    const UList<Type>& vf
);


template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const UList<Type>& f,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const tmp<Field<Type>>& tf,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);


template<class Type>
void sqr
(
    Field<typename outerProduct<Type, Type>::type>& res,
    const UList<Type>& vf
);

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type>>
sqr(const UList<Type>& f);

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type>>
sqr(const tmp<Field<Type>>& tf);


template<class Type>
void magSqr
(
    Field<typename typeOfMag<Type>::type>& res,
    const UList<Type>& f
);

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
magSqr(const UList<Type>& f);

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
magSqr(const tmp<Field<Type>>& tf);


template<class Type>
void mag
(
    Field<typename typeOfMag<Type>::type>& res,
    const UList<Type>& f
);

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
mag(const UList<Type>& f);

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
mag(const tmp<Field<Type>>& tf);


template<class Type>
void cmptMax(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMax(const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType>>
cmptMax(const tmp<Field<Type>>& tf);


template<class Type>
void cmptMin(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMin(const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType>>
cmptMin(const tmp<Field<Type>>& tf);


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptAv(const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptAv(const tmp<Field<Type>>& tf);


template<class Type>
void cmptMag(Field<Type>& res, const UList<Type>& f);

template<class Type>
tmp<Field<Type>> cmptMag(const UList<Type>& f);

template<class Type>
tmp<Field<Type>> cmptMag(const tmp<Field<Type>>& tf);


template<class Type>
void cmptMagSqr(Field<Type>& res, const UList<Type>& f);

template<class Type>
tmp<Field<Type>> cmptMagSqr(const UList<Type>& f);

template<class Type>
tmp<Field<Type>> cmptMagSqr(const tmp<Field<Type>>& tf);


#define TMP_UNARY_FUNCTION(ReturnType, Func)                                   \
                                                                               \
/*! \brief Apply the \c Func() function on the tmp field */                    \
template<class Type>                                                           \
ReturnType Func(const tmp<Field<Type>>& tf1);

template<class Type>
Type max(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, sum)


// From MinMaxOps.H:
//   - Foam::minMax(const UList<Type>&)
//   - Foam::minMaxMag(const UList<Type>&)

TMP_UNARY_FUNCTION(MinMax<Type>, minMax)
TMP_UNARY_FUNCTION(scalarMinMax, minMaxMag)


template<class Type>
Type maxMagSqr(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, maxMagSqr)

template<class Type>
Type minMagSqr(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, minMagSqr)


template<class Type>
typename scalarProduct<Type, Type>::type
sumProd
(
    const UList<Type>& f1,
    const UList<Type>& f2
);

template<class Type>
Type sumCmptProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
typename outerProduct1<Type>::type sumSqr(const UList<Type>& f);

template<class Type>
typename outerProduct1<Type>::type sumSqr(const tmp<Field<Type>>& tf);

template<class Type>
typename typeOfMag<Type>::type sumMag(const UList<Type>& f);

TMP_UNARY_FUNCTION(typename typeOfMag<Type>::type, sumMag)

template<class Type>
Type sumCmptMag(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, sumCmptMag)

//- The local arithmetic average of a field
template<class Type>
Type average(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, average)


// With reduction on ReturnType
#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                       \
                                                                               \
template<class Type>                                                           \
ReturnType gFunc(const UList<Type>& f, const label comm = UPstream::worldComm);\
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

G_UNARY_FUNCTION(MinMax<Type>, gMinMax, minMax, sum)
G_UNARY_FUNCTION(scalarMinMax, gMinMaxMag, minMaxMag, sum)

G_UNARY_FUNCTION(typename outerProduct1<Type>::type, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(typename typeOfMag<Type>::type, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION


template<class Type>
typename scalarProduct<Type, Type>::type gSumProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm = UPstream::worldComm
);

template<class Type>
Type gSumCmptProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm = UPstream::worldComm
);

//- The global arithmetic average of a field
template<class Type>
Type gAverage
(
    const UList<Type>& f,
    //! The communicator. A negative value disables the reduction
    const label comm = UPstream::worldComm
);

TMP_UNARY_FUNCTION(Type, gAverage)

//- The global weighted average of a field,
//- using the mag() of the weights.
template<class Type>
Type gWeightedAverage
(
    const UList<scalar>& weights,
    const UList<Type>& fld,
    //! The communicator. A negative value disables the reduction
    const label comm = UPstream::worldComm
);

//- The local weighted average of a field,
//- using the mag() of the weights.
template<class Type>
Type weightedAverage
(
    const UList<scalar>& weights,
    const UList<Type>& fld
)
{
    // Use the global version, but without reduction
    return gWeightedAverage(weights, fld, -1);
}

//- The global weighted sum (integral) of a field,
//- using the mag() of the weights.
template<class Type>
Type gWeightedSum
(
    const UList<scalar>& weights,
    const UList<Type>& fld,
    //! The communicator. A negative value disables the reduction
    const label comm = UPstream::worldComm
);

//- The local weighted sum (integral) of a field,
//- using the mag() of the weights.
template<class Type>
Type weightedSum
(
    const UList<scalar>& weights,
    const UList<Type>& fld
)
{
    // Use the global version, but without reduction
    return gWeightedSum(weights, fld, -1);
}


#undef TMP_UNARY_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION_FS(Type, Type, MinMax<Type>, clamp)
BINARY_TYPE_FUNCTION_FS(Type, Type, Foam::zero_one, clamp)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TERNARY_FUNCTION(Type, Type, Type, scalar, lerp)
TERNARY_TYPE_FUNCTION_FFS(Type, Type, Type, scalar, lerp)


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type1, Type2>::type>& res,                          \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
);                                                                             \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const UList<Type1>& f1, const UList<Type2>& f2);                   \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const UList<Type1>& f1, const tmp<Field<Type2>>& tf2);             \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<Field<Type1>>& tf1, const UList<Type2>& f2);             \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<Field<Type1>>& tf1, const tmp<Field<Type2>>& tf2);       \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type, Form>::type>& res,                            \
    const UList<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
tmp<Field<typename product<Type, Form>::type>>                                 \
operator Op(const UList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs);    \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
tmp<Field<typename product<Type, Form>::type>>                                 \
operator Op(const tmp<Field<Type>>&tf1,const VectorSpace<Form,Cmpt,nCmpt>&vs); \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Form, Type>::type>& res,                            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const UList<Type>& f1                                                      \
);                                                                             \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
tmp<Field<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const UList<Type>& f1);    \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
tmp<Field<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&vs,const tmp<Field<Type>>&tf1);

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"
#include "scalarField.H"

// ************************************************************************* //
