/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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
    Test-vector

Description
    Some simple tests for vector

\*---------------------------------------------------------------------------*/

#include "vectorField.H"
#include "boolVector.H"
#include "complexVector.H"
#include "labelVector.H"
#include "vector2D.H"
#include "IOstreams.H"
#include "FixedList.H"
#include "Random.H"
#include <algorithm>
#include <random>

using namespace Foam;

void printInfo(const vector& vec)
{
    Info<< "vector : " << vec << nl
        << "magSqr : " << magSqr(vec) << nl;

    Info<< "component"
        << " max:" << cmptMax(vec)
        << " sum:" << cmptSum(vec)
        << " prod:" << cmptProduct(vec)
        << " mag:" << cmptMag(vec)
        << " magSqr:" << cmptMagSqr(vec)
        << nl << nl;
}


void doTest(vector& vec1, vector& vec2)
{
    Info<<"normalised(vector1): " << normalised(vec1) << nl;
    Info<<"vector1: " << vec1 << nl;
    vec1.normalise();
    Info<<"normalised: " << vec1 << nl;

    vector vecsmall((1e-100 * vec1));

    Info<<"small: " << vecsmall << nl;
    Info<<"small normalised: " << normalised(vecsmall) << nl;
    Info<<"small diff: " << (vecsmall.normalise() - vec1) << nl;

    vec1 *= 4.0;

    Info<< "scalar mult: " << vec1 << nl;
    Info<< "addition   : " << (vec1 + vec1) << nl;

    printInfo(vec1);
    printInfo(vec2);

    Info<< "vector: " << vec1 << nl
        << "vector: " << vec2 << nl
        << "   min: " << min(vec1, vec2) << nl
        << "  dist: " << vec1.dist(vec2) << ' ' << mag(vec1 - vec2) << nl
        << "dist^2: " << vec1.distSqr(vec2) << ' ' << magSqr(vec1 - vec2) << nl
        << nl;
}


template<class VecSpace>
void testIterator(const VecSpace& vs)
{
    Info<< "size: " << vs.size() << " for:";
    for (const auto& val : vs)
    {
        Info<< " " << val;
    }
    Info<< nl;
}


template<class VecSpace>
void testData(const VecSpace& vs)
{
    Info<< "size: " << vs.size() << " for:";

    const auto* data = vs.cdata();
    for
    (
        const auto* endData = data + VecSpace::nComponents;
        data != endData;
        ++data
    )
    {
        Info<< " " << *data;
    }
    Info<< nl;
}


template<class Type>
void testNormalise(Field<Type>& fld)
{
    Info<< nl << pTraits<Type>::typeName << " Field" << nl
        << "  orig: " << fld << nl;
    fld.normalise();
    Info<< "  norm: " << fld << nl;
}


// Transcribe vectorspace information into a FixedList
template<class Type>
void testTranscribe(Type& input)
{
    if constexpr
    (
        is_vectorspace_v<Type>
     && std::is_floating_point_v<typename pTraits_cmptType<Type>::type>
    )
    {
        constexpr auto nCmpts = pTraits_nComponents_v<Type>;
        using cmpt = typename pTraits_cmptType<Type>::type;

        FixedList<cmpt, nCmpts+1> values;
        values.back() = 100;  // some additional data

        VectorSpaceOps<nCmpts>::copy_n(input.cdata(), values.data());

        Info<< "Transcribed " << input << " => " << values << nl;

        for (auto& val : values)
        {
            val *= -1;
        }

        VectorSpaceOps<nCmpts>::copy_n(values.cdata(), input.data());
        Info<< "  copied back (-1) as " << input
            << " from " << values << nl;
    }
    else
    {
        Info<< "Did not transcribe " << input << nl;
    }
}


// Test of adding to vectors, possible of dissimilar types
template<class Target, class Source>
void testAdding(const Target& a, const Source& b)
{
    typedef typename Target::cmptType cmptType1;
    typedef typename Source::cmptType cmptType2;

    Info<< "    "
        << pTraits<Target>::typeName << " += "
        << pTraits<Source>::typeName << " : ";

    if constexpr (std::is_convertible_v<cmptType2, cmptType1>)
    {
        Target res(a);
        res += b;
        Info<< a << " += " << b << " == " << res << nl;
    }
    else
    {
        Info<< "unsupported" << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<<"normalised: " << vector(1,2,3).normalise() << nl;
    Info<<"normalised: " << vector::uniform(VSMALL).normalise() << nl;
    Info<<"normalised: " << vector::uniform(ROOTVSMALL).normalise() << nl;

    {
        // VectorSpaceOps::copy_args
        vector vec;
        VectorSpaceOps<3>::copy_args(vec.data(), 1, 2, 3);
        Info<< "copy_args: " << vec << nl;

        VectorSpaceOps<3, 1>::copy_args(vec.data(), 4, 5);
        Info<< "copy_args: " << vec << nl;

        // completely ignored:
        // VectorSpaceOps<3, 10>::copy_args(vec.data(), 4, 5);
        // Info<< "copy_args: " << vec << nl;

        double a[2]{};
        VectorSpaceOps<2>::copy_args(a, 4, 5.5);
        Info<< "copy_args: " << a[0] << ' ' << a[1] << nl;
    }

    {
        complexVector cvec(complex(1), complex(0), complex(0));
        doubleVector dvec(1, 0, 0);
        floatVector  fvec(1, 0, 0);
        labelVector  lvec(1, 0, 0);

        Info<< nl << "additions:" << nl;

        testAdding(cvec, cvec);
        testAdding(cvec, dvec);
        testAdding(cvec, fvec);
        testAdding(cvec, lvec);

        testAdding(dvec, cvec);
        testAdding(dvec, dvec);
        testAdding(dvec, fvec);
        testAdding(dvec, lvec);

        testAdding(fvec, cvec);
        testAdding(fvec, dvec);
        testAdding(fvec, fvec);
        testAdding(fvec, lvec);

        testAdding(lvec, cvec);
        testAdding(lvec, dvec);
        testAdding(lvec, fvec);
        testAdding(lvec, lvec);

        Info<< nl;
    }

    {
        vector vec1(0.5, 0.5, 0.5);
        vector vec2(0.5, 0.51, -0.5);

        doTest(vec1, vec2);

        testIterator(vec1);
        testIterator(vec2);
        testData(vec2);

        // Use STL algorithm(s)

        std::sort(vec2.begin(), vec2.end());
        Info<< "sorted: " << vec2 << nl;

        std::shuffle(vec2.begin(), vec2.end(), std::default_random_engine());
        Info<< "shuffled: " << vec2 << nl;

        // Vectors with some identical components
        vectorField vectors
        ({
            {1.1, 2.2, 3.3 },
            {2.2, 3.3, 4.4 },
            {-1.1, 2.2, 3.3 },
            {-2.2, 3.3, 4.4 },

            {-1.1, -2.2, 3.3 },
            {-2.2, -3.3, 4.4 },

            {-1.1, -2.2, -3.3 },
            {-2.2, -3.3, -4.4 },
            {-3.3, 2.1, 12 },
            {3.145, 1.6, 2 },

            {0, 0, 0}
        });

        shuffle(vectors);

        Info<< "initial vectors: ";
        vectors.writeList(Info, 1) << nl;

        Foam::sort(vectors);
        Info<< "regular sort:";
        vectors.writeList(Info, 1) << nl;

        std::sort(vectors.begin(), vectors.end(), vector::less_xyz);
        Info<< "sorted xyz:";
        vectors.writeList(Info, 1) << nl;

        std::sort(vectors.begin(), vectors.end(), vector::less_yzx);
        Info<< "sorted yzx:";
        vectors.writeList(Info, 1) << nl;

        std::sort(vectors.begin(), vectors.end(), vector::less_zxy);
        Info<< "sorted zxy:";
        vectors.writeList(Info, 1) << nl;

        vectorField vecCmptMag1(cmptMag(vectors));
        Info<< "cmptMag:";
        vecCmptMag1.writeList(Info, 1) << nl;

        vectorField vecCmptMag2(vectors.size());
        vecCmptMag2.replace(vector::X, mag(vectors.component(vector::X)));
        vecCmptMag2.replace(vector::Y, mag(vectors.component(vector::Y)));
        vecCmptMag2.replace(vector::Z, mag(vectors.component(vector::Z)));

        Info<< "cmptMag:";
        vecCmptMag2.writeList(Info, 1) << nl;

        Info<< "mult:";
        cmptMultiply(vecCmptMag2, vecCmptMag2, vector(2,3,4));
        vecCmptMag2.writeList(Info, 1) << nl;
    }
    // Basic tests for fields
    {
        scalarField sfld1
        ({
            0.0,
            1.0,
            2.0
        });

        testNormalise(sfld1);

        Field<floatVector> vfld1
        ({
            floatVector(0.0, 1.0, 2.0),
            floatVector(0, 0, floatScalarSMALL),
            floatVector(0, 0, floatScalarVSMALL),
            floatVector(0, 2, 1),
        });

        testNormalise(vfld1);

        Field<doubleVector> vfld2
        ({
            doubleVector(0.0, 1.0, 2.0),
            doubleVector(0, 0, doubleScalarSMALL),
            doubleVector(0, 0, doubleScalarVSMALL),
            doubleVector(0, 2, 1),
        });

        testNormalise(vfld2);
    }

    Info<< nl
        << "Test constexpr evaluation" << nl;
    {
        // constexpr: vector magSqr (free function)
        constexpr scalar magSqr123 = magSqr(vector(1, 2, 3));
        static_assert(magSqr123 == 14, "magSqr(vector) should be constexpr");
        Info<< "constexpr magSqr(vector(1,2,3)) = " << magSqr123 << nl;

        // constexpr: Vector2D magSqr()/distSqr() member functions
        constexpr vector2D v2a(3, 4);
        constexpr vector2D v2b(0, 0);

        constexpr scalar magSqr2D = v2a.magSqr();
        constexpr scalar distSqr2D = v2a.distSqr(v2b);

        static_assert(magSqr2D == 25, "Vector2D::magSqr() should be constexpr");
        static_assert
        (
            distSqr2D == 25,
            "Vector2D::distSqr() should be constexpr"
        );

        Info<< "constexpr vector2D(3,4).magSqr() = " << magSqr2D << nl
            << "constexpr vector2D(3,4).distSqr(vector2D(0,0)) = "
            << distSqr2D << nl;
    }

    Info<< nl
        << "Test transcribing components" << nl;
    {
        vector vec1(1.1, 2.2, 3.3);
        testTranscribe(vec1);

        labelVector vec2(10, 20, 30);
        testTranscribe(vec2);
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
