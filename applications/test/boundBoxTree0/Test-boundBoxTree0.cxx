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

Description
    Test boundBoxTree traversal, nearest searches, and stream extraction

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "boundBoxTree.H"
#include "Random.H"
#include "StringStream.H"

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


void testFindNearest()
{
    List<treeBoundBox> boxes(1);
    boxes[0] = treeBoundBox(point(1, 1, 1), point(2, 2, 2));

    boundBoxTree tree(boxes);
    DynamicList<label> hits;

    const auto overlap = tree.findNearest(point(0, 1.5, 1.5), 1.1, hits);

    check
    (
        overlap.second() == volumeType::MIXED,
        "Expected the sphere to overlap the box"
    );
    check
    (
        hits.size() == 1 && hits[0] == 0,
        "Expected the overlapping box in the nearest hits"
    );

    const auto tangent = tree.findNearest(point(0, 1.5, 1.5), 1, hits);

    check
    (
        tangent.second() == volumeType::MIXED,
        "Expected an exact tangent to count as an overlap"
    );
    check
    (
        hits.size() == 1 && hits[0] == 0,
        "Expected the tangent box in the nearest hits"
    );

    List<treeBoundBox> orderBoxes(2);
    orderBoxes[0] =
        treeBoundBox(point(1, -10, -10), point(2, 10, 10));
    orderBoxes[1] =
        treeBoundBox(point(2, 0, 0), point(3, 0.1, 0.1));

    boundBoxTree orderTree(orderBoxes);
    const auto ordered =
        orderTree.findNearest(point(0, 0, 0), 4.1, hits);

    check
    (
        ordered.second() == volumeType::MIXED,
        "Expected both boxes to overlap the search sphere"
    );
    check
    (
        hits.size() == 2 && hits[0] == 0 && hits[1] == 1,
        "Expected boxes to be sorted by point-to-box distance"
    );
}


// void testUnsupportedExtraction()
// {
//     List<treeBoundBox> boxes(1);
//     boxes[0] = treeBoundBox(point(1, 1, 1), point(2, 2, 2));
//
//     boundBoxTree tree(boxes);
//     IStringStream input("unsupported input");
//     input >> tree;
//
//     check(input.fail(), "Expected unsupported extraction to fail");
//     check
//     (
//         tree.size() == 1
//      && tree[0].min() == boxes[0].min()
//      && tree[0].max() == boxes[0].max(),
//         "Expected failed extraction to leave the tree unchanged"
//     );
// }


// Visitor tallying leaf/interior nodes and total indexed boxes
struct treeStatsVisitor
{
    label nLeaves = 0;
    label nInterior = 0;
    label nIndices = 0;

    void leaf(const treeBoundBox&, const labelUList& indices)
    {
        ++nLeaves;
        nIndices += indices.size();
    }

    bool interior(const treeBoundBox&)
    {
        ++nInterior;
        return true;  // always descend
    }
};


// Visitor collecting indices of all leaves overlapping a query box
struct leafOverlapVisitor
{
    const treeBoundBox& queryBb;
    DynamicList<label> hits;

    void leaf(const treeBoundBox&, const labelUList& indices)
    {
        hits.push_back(indices);
    }

    bool interior(const treeBoundBox& bb)
    {
        return bb.overlaps(queryBb);
    }
};


int main(int argc, char* argv[])
{
    argList::noBanner();
    argList::noParallel();
    #include "setRootCase.H"

    testFindNearest();
    // testUnsupportedExtraction();

    Random rndGen(1234);

    List<treeBoundBox> boxes(200);
    for (auto& bb : boxes)
    {
        const point lo
        (
            rndGen.sample01<scalar>(),
            rndGen.sample01<scalar>(),
            rndGen.sample01<scalar>()
        );
        bb = treeBoundBox(lo, lo + vector(0.01, 0.01, 0.01));
    }

    boundBoxTree tree
    (
        boundBoxTree::algorithm::median_split,
        List<treeBoundBox>(boxes),
        8,   // maxLevel
        6    // minLeafSize
    );

    Info<< "nodes: " << tree.nNodes() << nl;

    treeStatsVisitor stats;
    tree.visit(stats);

    Info<< "leaves: " << stats.nLeaves
        << " interior: " << stats.nInterior << nl
        << "total indexed boxes: " << stats.nIndices
        << " (expect " << boxes.size() << ")" << nl;

    const treeBoundBox queryBb(point(0, 0, 0), point(0.3, 0.3, 0.3));
    leafOverlapVisitor collector{queryBb, DynamicList<label>()};
    tree.visit(collector);

    Info<< "boxes in leaves overlapping " << queryBb << ": "
        << collector.hits.size() << nl
        << "overlaps(query): " << tree.overlaps(queryBb) << nl;

    return 0;
}


// ************************************************************************* //
