/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application
    Test-contiguous2

Description
    Test contiguous type traits (primarily for vector-space types)

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "wordRes.H"
#include "contiguous.H"
#include "error.H"

#include "IOstreams.H"
#include "scalar.H"
#include "vector.H"
#include "vector2D.H"
#include "tensor.H"
#include "tensor2D.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "diagTensor.H"
#include "barycentric.H"
#include "barycentricTensor.H"
#include "spatialVector.H"
#include "spatialTensor.H"
#include "triad.H"
#include "treeBoundBox.H"
#include "Switch.H"

#include "labelRange.H"
#include "scalarList.H"
#include "HashOps.H"
#include "FixedList.H"
#include "Pair.H"
#include "Tuple2.H"

namespace Foam
{

// Wrong, but interesting to test
template<> struct is_contiguous<Pair<word>> : std::true_type {};

} // end namespace Foam

using namespace Foam;


template<class T>
void printInfo(const char* const name = nullptr)
{
    if (name == nullptr)
    {
        Info<< error::demangle<T>();
    }
    else
    {
        Info<< name;
    }

    if constexpr (is_contiguous_v<T>)
    {
        Info<< " contiguous";
    }
    else
    {
        Info<< " NO";
    }

    if constexpr (is_contiguous_label<T>::value)
    {
        Info<< "(label)";
    }
    if constexpr (is_contiguous_scalar<T>::value)
    {
        Info<< "(scalar)";
    }

    {
        using elemType = contiguous_element_t<T>;

        Info<< " element=" << error::demangle<elemType>();
    }

    Info<< nl;
}


template<class cmptType>
void test_vectorspace()
{
    // General VectorSpace / MatrixSpace derived types,
    // all parameterized on the same component (cmptType)

    printInfo<Vector<cmptType>>();
    printInfo<Vector2D<cmptType>>();
    printInfo<Tensor<cmptType>>();
    printInfo<Tensor2D<cmptType>>();
    printInfo<SymmTensor<cmptType>>();
    printInfo<SphericalTensor<cmptType>>();
    printInfo<DiagTensor<cmptType>>();
    printInfo<Barycentric<cmptType>>();
    printInfo<BarycentricTensor<cmptType>>();
    printInfo<SpatialVector<cmptType>>();
    printInfo<SpatialTensor<cmptType>>();

    // This only exists in scalar (double) form
    if constexpr (std::is_same_v<scalar, cmptType>)
    {
        printInfo<triad>();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();

    Info<< "report [type, contiguous, ..., element]" << nl
        << endl;

    printInfo<bool>();
    printInfo<int32_t>();
    printInfo<int64_t>();
    printInfo<float>();
    printInfo<double>();
    printInfo<FixedList<double, 4>>();
    printInfo<const Pair<long>>();

    printInfo<FixedList<word, 2>>();
    printInfo<Pair<word>>();

    printInfo<std::pair<int, vector>>();

    printInfo<FixedList<FixedList<int, 2>, 2>>();

    printInfo<treeBoundBox>();
    printInfo<labelPair>();
    printInfo<Tuple2<labelPair, labelPair>>();

    Info<< nl << "VectorSpace / MatrixSpace derived types (float)" << nl
        << endl;

    test_vectorspace<float>();

    Info<< nl << "VectorSpace / MatrixSpace derived types (double)" << nl
        << endl;

    test_vectorspace<double>();

    Info<< nl << "Nested/mixed VectorSpace compositions" << nl << endl;

    // A Vector of Vectors: cmptType is itself a VectorSpace type
    printInfo<Vector<vector>>();

    // Mixing a fundamental type with a VectorSpace type of the same
    // underlying (float) element - contiguous_element should resolve
    // to a common 'float'
    printInfo<Tuple2<float, floatVector>>();

    // Mixing two different VectorSpace types (vector, tensor) that both
    // ultimately resolve to 'scalar'
    printInfo<Tuple2<vector, tensor>>();

    // Mismatched underlying element types (scalar vs label) - the
    // Tuple2 contiguous_element falls back to a pass-through (non-useful)
    // type since there is no common element
    printInfo<Tuple2<vector, label>>();

    // std::tuple counterparts to the Tuple2 cases above - same common
    // element resolution (or fallback), but for an arbitrary arity
    printInfo<std::tuple<float, float, floatVector>>();
    printInfo<std::tuple<vector, tensor>>();
    printInfo<std::tuple<vector, label>>();

    return 0;
}

// ************************************************************************* //
