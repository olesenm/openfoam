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
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-AABBTree

Description
    Deterministic smoke tests for AABBTree.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "AABBTree.H"
#include "clockTime.H"
#include "faceList.H"
#include "pointField.H"
#include "StringStream.H"
#include <limits>

using namespace Foam;


void check(const bool condition, const char* message)
{
    if (!condition)
    {
        FatalErrorInFunction
            << message
            << abort(FatalError);
    }
}


void testNodeAddressing(const pointField& points)
{
    AABBTree<point> oneLeafTree
    (
        points,
        true,   // median split
        2,      // max level
        8       // min leaf size
    );

    check(oneLeafTree.nLeafs() == 1, "Expected one leaf");
    oneLeafTree.reindex();

    const auto& oneLeafNodes = oneLeafTree.nodes();
    check(oneLeafNodes.size() == 1, "Expected one node for one leaf");
    check(oneLeafNodes[0].is_leaf(), "Expected one-leaf root to be a leaf");

    AABBTree<point> multiLeafTree
    (
        points,
        true,   // median split
        1,      // max level
        1       // min leaf size
    );

    check(multiLeafTree.nLeafs() == 2, "Expected two leaves");

    const auto& multiLeafNodes = multiLeafTree.nodes();
    check(multiLeafNodes.size() == 3, "Expected three nodes for two leaves");
    check
    (
        !multiLeafNodes[0].is_leaf()
     && multiLeafNodes[1].is_leaf()
     && multiLeafNodes[2].is_leaf(),
        "Expected a binary root with two leaf children"
    );
}


void testQueries()
{
    pointField points(2);
    points[0] = point(10, 0, 0);
    points[1] = point(0, 0, 0);

    AABBTree<point> tree
    (
        points,
        true,   // median split
        1,      // max level
        1       // min leaf size
    );

    check(tree.nLeafs() == 2, "Expected two query leaves");
    check(tree.pointInside(points[0]), "Expected high point query to match");
    check(tree.pointInside(points[1]), "Expected low point query to match");
    check
    (
        !tree.pointInside(point(5, 0, 0)),
        "Expected point in the gap between leaves not to match"
    );

    const boundBox lowQuery
    (
        point(-0.1, -0.1, -0.1),
        point(0.1, 0.1, 0.1)
    );
    const boundBox highQuery
    (
        point(9.9, -0.1, -0.1),
        point(10.1, 0.1, 0.1)
    );
    const boundBox gapQuery
    (
        point(4.9, -0.1, -0.1),
        point(5.1, 0.1, 0.1)
    );

    check(tree.overlaps(lowQuery), "Expected low overlap query to match");
    check(tree.overlaps(highQuery), "Expected high overlap query to match");
    check(!tree.overlaps(gapQuery), "Expected gap overlap query not to match");
}


void testObjectCentreConstruction()
{
    pointField points(5);
    points[0] = point(0, 0, 0);
    points[1] = point(1, 0, 0);
    points[2] = point(2, 0, 0);
    points[3] = point(10, 0, 0);
    points[4] = point(11, 0, 0);

    faceList faces(2);
    faces[0] = face({0, 1, 2, 2});
    faces[1] = face({2, 3, 4});

    AABBTree<face> tree
    (
        faces,
        points,
        true,   // median split
        1,      // max level
        1       // min leaf size
    );

    check(tree.nLeafs() == 2, "Expected two object leaves");
    check
    (
        tree.pointInside(points[0]) && tree.pointInside(points[4]),
        "Expected object leaves to include their referenced points"
    );
}


void testLargeMinLeafSize(const pointField& points)
{
    AABBTree<point> tree
    (
        points,
        true,   // median split
        1,      // max level
        std::numeric_limits<int>::max()
    );

    check(tree.nLeafs() == 1, "Large min leaf size should stop splitting");
}


void testStreamReadInvalidatesNodes()
{
    pointField oldPoints(2);
    oldPoints[0] = point(0, 0, 0);
    oldPoints[1] = point(1, 0, 0);

    AABBTree<point> tree
    (
        oldPoints,
        true,   // median split
        1,      // max level
        1       // min leaf size
    );

    check(tree.pointInside(oldPoints[0]), "Expected initial query to match");
    check(tree.hasNodes(), "Expected the initial query to cache nodes");

    pointField newPoints(2);
    newPoints[0] = point(100, 0, 0);
    newPoints[1] = point(101, 0, 0);

    AABBTree<point> replacement
    (
        newPoints,
        true,   // median split
        1,      // max level
        1       // min leaf size
    );

    OStringStream output;
    output << replacement;

    IStringStream input(output.str());
    input >> tree;

    check(!tree.hasNodes(), "Expected stream extraction to invalidate nodes");
    check
    (
        !tree.pointInside(oldPoints[0]),
        "Expected the old tree's point not to match after extraction"
    );
    check
    (
        tree.pointInside(newPoints[0]),
        "Expected the extracted tree's point to match"
    );
}


void runCentreScalingBenchmark()
{
    const label objectCounts[] = {2000, 4000, 8000, 16000};

    for (const label nObjects : objectCounts)
    {
        pointField points(3*nObjects);
        faceList faces(nObjects);

        for (label objecti = 0; objecti < nObjects; ++objecti)
        {
            const label firstPoint = 3*objecti;
            const scalar x = 2*scalar(objecti);

            points[firstPoint] = point(x, 0, 0);
            points[firstPoint + 1] = point(x + 1, 0, 0);
            points[firstPoint + 2] = point(x, 1, 0);

            faces[objecti] =
                face({firstPoint, firstPoint + 1, firstPoint + 2});
        }

        clockTime timing;
        AABBTree<face> tree
        (
            faces,
            points,
            true,   // median split
            12,     // max level
            4       // min leaf size
        );

        Info<< "Objects: " << nObjects
            << " leaves: " << tree.nLeafs()
            << " construction time: " << timing.elapsedTime() << " s"
            << nl;
    }
}


int main(int argc, char* argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addBoolOption
    (
        "benchmark",
        "run the centre-estimation scaling benchmark"
    );
    argList args(argc, argv);

    pointField points(4);
    points[0] = point(0, 0, 0);
    points[1] = point(1, 0, 0);
    points[2] = point(0, 1, 0);
    points[3] = point(1, 1, 0);

    AABBTree<point> tree
    (
        points,
        true,   // median split
        2,      // max level
        8       // min leaf size
    );

    check(tree.nLeafs() == 1, "Expected one AABBTree leaf");
    check(tree.leafBoxes().size() == 1, "Expected one leaf box");
    check
    (
        tree.indices(0).size() == points.size(),
        "Expected all points in the leaf"
    );

    testNodeAddressing(points);
    testQueries();
    testObjectCentreConstruction();
    testLargeMinLeafSize(points);
    testStreamReadInvalidatesNodes();

    if (args.found("benchmark"))
    {
        runCentreScalingBenchmark();
    }

    Info<< "AABBTree tests passed" << endl;
    return 0;
}


// ************************************************************************* //
