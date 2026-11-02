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

Description
    Template-like definitions for Field operations that can be offloaded
    to GPU.

SourceFiles
    FieldOffload.H

\*---------------------------------------------------------------------------*/

#ifdef FOAM_OFFLOAD

#include "FieldOffload.H"

#include "scalar.H"
#include "bool.H"
#include "label.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"

#include <execution>
#include <numeric>

namespace Foam
{
    // -------------------- max
    template<class Type>
    Type maxOffload(const UList<Type>& f1)
    {
        Type result = f1[0];

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                f1[0],
                [] __host__ __device__ (const Type& a, const Type& b) -> Type
                {
                    return max(a, b);
                },
                [] __host__ __device__ (const Type& x) -> Type
                {
                    return x;
                }
            );
        });

        return result;
    }

    // explicit instantiations as needed
    template scalar maxOffload<scalar>(const UList<scalar>&);
    template label  maxOffload<label>(const UList<label>&);
    template vector maxOffload<vector>(const UList<vector>&);
    template sphericalTensor maxOffload<sphericalTensor>(const UList<sphericalTensor>&);
    template symmTensor maxOffload<symmTensor>(const UList<symmTensor>&);
    template tensor maxOffload<tensor>(const UList<tensor>&);


    //------------------- min
    template<class Type>
    Type minOffload(const UList<Type>& f1)
    {
        Type result = f1[0];

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                f1[0],
                [] __host__ __device__ (const Type& a, const Type& b) -> Type
                {
                    return min(a, b);
                },
                [] __host__ __device__ (const Type& x) -> Type
                {
                    return x;
                }
            );
        });

        return result;
    }

    // explicit instantiations as needed
    template scalar minOffload<scalar>(const UList<scalar>&);
    template label  minOffload<label>(const UList<label>&);
    template vector minOffload<vector>(const UList<vector>&);
    template sphericalTensor minOffload<sphericalTensor>(const UList<sphericalTensor>&);
    template symmTensor minOffload<symmTensor>(const UList<symmTensor>&);
    template tensor minOffload<tensor>(const UList<tensor>&);


    // ------------------- sum
    template<class Type>
    typename typeOfSolve<Type>::type
    sumOffload(const UList<Type>& f)
    {
        using resultType = typename typeOfSolve<Type>::type;

        resultType result = resultType(Zero);

        apply_policy(f.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f.begin(),
                f.end(),
                resultType(Zero),
                [] __host__ __device__ (const resultType& a, const resultType& b) -> resultType
                {
                    return a + b;
                },
                [] __host__ __device__ (const Type& x) -> resultType
                {
                    return x;
                }
            );
        });

        return result;
    }

    // explicit instantiations as needed
    template typeOfSolve<float>::type
    sumOffload<float>(const UList<float>&);

    template typeOfSolve<double>::type
    sumOffload<double>(const UList<double>&);

    template typeOfSolve<int>::type
    sumOffload<int>(const UList<int>&);

    template typeOfSolve<bool>::type
    sumOffload<bool>(const UList<bool>&);

    template typeOfSolve<vector>::type
    sumOffload<vector>(const UList<vector>&);

    template typeOfSolve<sphericalTensor>::type
    sumOffload<sphericalTensor>(const UList<sphericalTensor>&);

    template typeOfSolve<symmTensor>::type
    sumOffload<symmTensor>(const UList<symmTensor>&);

    template typeOfSolve<tensor>::type
    sumOffload<tensor>(const UList<tensor>&);


    // ------------------- sumProd
    template<class Type>
    typename scalarProduct<Type, Type>::type
    sumProdOffload(const UList<Type>& f1, const UList<Type>& f2)
    {
        using resultType = typename scalarProduct<Type, Type>::type;
        resultType result = resultType(Zero);

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                f2.begin(),
                resultType(Zero),
                [] __host__ __device__ (const resultType& a, const resultType& b) -> resultType
                {
                    return a + b;
                },
                [] __host__ __device__ (const Type& a, const Type& b) -> resultType
                {
                    // For float and double this is not needed because scalarField.C contains their definition
                    // if constexpr (std::is_arithmetic_v<Type>)
                    // {
                    //     return a * b;        // scalar case
                    // }
                    // else
                    {
                        return (a && b);     // OpenFOAM inner product
                    }
                }
            );
        });

        return result;
    }


    // explicit instantiations as needed
    template scalar sumProdOffload<scalar>(const UList<scalar>&, const UList<scalar>&);
    template scalar sumProdOffload<vector>(const UList<vector>&, const UList<vector>&);
    template scalar sumProdOffload<sphericalTensor>(const UList<sphericalTensor>&, const UList<sphericalTensor>&);
    template scalar sumProdOffload<symmTensor>(const UList<symmTensor>&, const UList<symmTensor>&);
    template scalar sumProdOffload<tensor>(const UList<tensor>&, const UList<tensor>&);


    // ------------------- sumProdFloatOffload
    float Foam::sumProdScalarOffload(const UList<float>& f1, const UList<float>& f2)
    {
        float result = 0.0;

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                f2.begin(),
                0.0f,
                std::plus<float>{},
                std::multiplies<float>{}
            );
        });

        return result;
    }


    // ------------------- sumProdDoubleOffload
    double Foam::sumProdScalarOffload(const UList<double>& f1, const UList<double>& f2)
    {
        double result = 0.0;

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                f2.begin(),
                0.0,
                std::plus<double>{},
                std::multiplies<double>{}
            );
        });

        return result;
    }


    // ------------------- maxMagSqr
    template<class Type>
    Type maxMagSqrOffload(const UList<Type>& f1)
    {
        using magType = decltype(magSqr(std::declval<Type>()));

        struct MaxMagSqrData
        {
            Type value;
            magType metric;
        };

        MaxMagSqrData result{Zero, magType(-1)};

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                MaxMagSqrData{Zero, magType(-1)},
                [] __host__ __device__ (const MaxMagSqrData& a, const MaxMagSqrData& b) -> MaxMagSqrData
                {
                    return (a.metric < b.metric) ? b : a;
                },
                [] __host__ __device__ (const Type& x) -> MaxMagSqrData
                {
                    return MaxMagSqrData{x, magSqr(x)};
                }
            );
        });

        return result.value;
    }

    // explicit instantiations as needed
    template scalar
    maxMagSqrOffload<scalar>(const UList<scalar>&);

    template vector
    maxMagSqrOffload<vector>(const UList<vector>&);


    // ------------------- minMagSqr
    template<class Type>
    Type minMagSqrOffload(const UList<Type>& f1)
    {
        using magType = decltype(magSqr(std::declval<Type>()));

        struct MinMagSqrData
        {
            Type value;
            magType metric;
        };

        MinMagSqrData result{Zero, magType(-1)};

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                MinMagSqrData{Zero, magType(-1)},
                [] __host__ __device__ (const MinMagSqrData& a, const MinMagSqrData& b) -> MinMagSqrData
                {
                    return (a.metric > b.metric) ? b : a;
                },
                [] __host__ __device__ (const Type& x) -> MinMagSqrData
                {
                    return MinMagSqrData{x, magSqr(x)};
                }
            );
        });

        return result.value;
    }

    // explicit instantiations as needed
    template scalar
    minMagSqrOffload<scalar>(const UList<scalar>&);

    template vector
    minMagSqrOffload<vector>(const UList<vector>&);


    // ------------------- sumCmptProd
    template<class Type>
    Type sumCmptProdOffload(const UList<Type>& f1, const UList<Type>& f2)
    {
        Type result = Type(Zero);

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                f2.begin(),
                Type(Zero),
                [] __host__ __device__ (const Type& a, const Type& b) -> Type
                {
                    return a + b;
                },
                [] __host__ __device__ (const Type& a, const Type& b) -> Type
                {
                    return cmptMultiply(a, b);
                }
            );
        });

        return result;
    }

    // explicit instantiations as needed
    template scalar sumCmptProdOffload<scalar>(const UList<scalar>&, const UList<scalar>&);
    template vector sumCmptProdOffload<vector>(const UList<vector>&, const UList<vector>&);
    template sphericalTensor sumCmptProdOffload<sphericalTensor>(const UList<sphericalTensor>&, const UList<sphericalTensor>&);
    template symmTensor sumCmptProdOffload<symmTensor>(const UList<symmTensor>&, const UList<symmTensor>&);
    template tensor sumCmptProdOffload<tensor>(const UList<tensor>&, const UList<tensor>&);


    // ------------------- sumSqr
    template<class Type>
    typename outerProduct1<Type>::type
    sumSqrOffload(const UList<Type>& f1)
    {
        using resultType = typename outerProduct1<Type>::type;

        resultType result = resultType(Zero);

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                resultType(Zero),
                [] __host__ __device__ (const resultType& a, const resultType& b) -> resultType
                {
                    return a + b;
                },
                [] __host__ __device__ (const Type& x) -> resultType
                {
                    return sqr(x);
                }
            );
        });

        return result;
    }

    // explicit instantiations as needed
    template outerProduct1<float>::type
    sumSqrOffload<float>(const UList<float>&);

    template outerProduct1<double>::type
    sumSqrOffload<double>(const UList<double>&);

    // template outerProduct1<vector>::type
    // sumSqrOffload<vector>(const UList<vector>&);

    // template outerProduct1<sphericalTensor>::type
    // sumSqrOffload<sphericalTensor>(const UList<sphericalTensor>&);

    // template outerProduct1<symmTensor>::type
    // sumSqrOffload<symmTensor>(const UList<symmTensor>&);

    // template outerProduct1<tensor>::type
    // sumSqrOffload<tensor>(const UList<tensor>&);


    // ------------------- sumMag
    template<class Type>
    typename typeOfMag<Type>::type
    sumMagOffload(const UList<Type>& f1)
    {
        using resultType = typename typeOfMag<Type>::type;

        resultType result = resultType(Zero);

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                resultType(Zero),
                [] __host__ __device__ (const resultType& a, const resultType& b) -> resultType
                {
                    return a + b;
                },
                [] __host__ __device__ (const Type& x) -> resultType
                {
                    return mag(x);
                }
            );
        });

        return result;
    }

    // explicit instantiations as needed
    template typeOfMag<float>::type
    sumMagOffload<float>(const UList<float>&);

    template typeOfMag<double>::type
    sumMagOffload<double>(const UList<double>&);

    template typeOfMag<vector>::type
    sumMagOffload<vector>(const UList<vector>&);

    template typeOfMag<sphericalTensor>::type
    sumMagOffload<sphericalTensor>(const UList<sphericalTensor>&);

    template typeOfMag<symmTensor>::type
    sumMagOffload<symmTensor>(const UList<symmTensor>&);

    template typeOfMag<tensor>::type
    sumMagOffload<tensor>(const UList<tensor>&);


    // ------------------- sumCmptMag
    template<class Type>
    Type sumCmptMagOffload(const UList<Type>& f1)
    {
        Type result = Type(Zero);

        apply_policy(f1.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                f1.begin(),
                f1.end(),
                Type(Zero),
                [] __host__ __device__ (const Type& a, const Type& b) -> Type
                {
                    return a + b;
                },
                [] __host__ __device__ (const Type& x) -> Type
                {
                    return cmptMag(x);
                }
            );
        });

        return result;
    }

    // explicit instantiations as needed
    template scalar sumCmptMagOffload<scalar>(const UList<scalar>&);
    template vector sumCmptMagOffload<vector>(const UList<vector>&);
    template sphericalTensor sumCmptMagOffload<sphericalTensor>(const UList<sphericalTensor>&);
    template symmTensor sumCmptMagOffload<symmTensor>(const UList<symmTensor>&);
    template tensor sumCmptMagOffload<tensor>(const UList<tensor>&);


    // ------------------- gWeightedAverage
    template<class Type>
    void gWeightedAverageOffload
    (
        const UList<scalar>& weights,
        const UList<Type>& fld,
        scalar& weight,
        Type& result
    )
    {
        apply_policy(weights.size(), [&](auto policy)
        {
            weight = std::transform_reduce
            (
                policy,
                weights.begin(),
                weights.end(),
                scalar(0),
                [] __host__ __device__ (const scalar& a, const scalar& b) -> scalar
                {
                    return a + b;
                },
                [] __host__ __device__ (const scalar& x) -> scalar
                {
                    return mag(x);
                }
            );
        });

        apply_policy(weights.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                weights.begin(),
                weights.end(),
                fld.begin(),
                Type(Zero),
                [] __host__ __device__ (const Type& a, const Type& b) -> Type
                {
                    return a + b;
                },
                [] __host__ __device__ (const scalar& w, const Type& f) -> Type
                {
                    return mag(w)*f;
                }
            );
        });
    }

    // explicit instantiations as needed
    template void gWeightedAverageOffload<scalar>
    (const UList<scalar>&, const UList<scalar>&, scalar&, scalar&);

    template void gWeightedAverageOffload<vector>
    (const UList<scalar>&, const UList<vector>&, scalar&, vector&);

    template void gWeightedAverageOffload<sphericalTensor>
    (const UList<scalar>&, const UList<sphericalTensor>&, scalar&, sphericalTensor&);

    template void gWeightedAverageOffload<symmTensor>
    (const UList<scalar>&, const UList<symmTensor>&, scalar&, symmTensor&);

    template void gWeightedAverageOffload<tensor>
    (const UList<scalar>&, const UList<tensor>&, scalar&, tensor&);


    // ------------------- gWeightedSum
    template<class Type>
    void gWeightedSumOffload
    (
        const UList<scalar>& weights,
        const UList<Type>& fld,
        Type& result
    )
    {
        apply_policy(weights.size(), [&](auto policy)
        {
            result = std::transform_reduce
            (
                policy,
                weights.begin(),
                weights.end(),
                fld.begin(),
                Type(Zero),
                [] __host__ __device__ (const Type& a, const Type& b) -> Type
                {
                    return a + b;
                },
                [] __host__ __device__ (const scalar& w, const Type& f) -> Type
                {
                    return mag(w) * f;
                }
            );
        });
    }

    // explicit instantiations as needed
    template void gWeightedSumOffload<scalar>
    (const UList<scalar>&, const UList<scalar>&, scalar&);

    template void gWeightedSumOffload<vector>
    (const UList<scalar>&, const UList<vector>&, vector&);

    template void gWeightedSumOffload<sphericalTensor>
    (const UList<scalar>&, const UList<sphericalTensor>&, sphericalTensor&);

    template void gWeightedSumOffload<symmTensor>
    (const UList<scalar>&, const UList<symmTensor>&, symmTensor&);

    template void gWeightedSumOffload<tensor>
    (const UList<scalar>&, const UList<tensor>&, tensor&);
}

#endif // FOAM_OFFLOAD

// ************************************************************************* //

