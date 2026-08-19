/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    surfaceRedistributePar

Group
    grpSurfaceUtilities

Description
    (Re)distribution of triSurface. Either takes an undecomposed surface
    or an already decomposed surface and redistributes it so that each
    processor has all triangles that overlap its mesh.

Note
    - best decomposition option is hierarchical since it guarantees
      square decompositions.
    - triangles might be present on multiple processors.
    - merging uses geometric tolerance so take care with writing precision.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "distributedTriSurfaceMesh.H"
#include "mapDistribute.H"
#include "decompositionModel.H"
#include "foamVtkSurfaceWriter.H"
#include "foamVtuWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Print on master all the per-processor surface stats.
// Determine surface bounding boxes, faces, points
void writeProcStats
(
    const triSurface& s,
    const UList<List<treeBoundBox>>& meshBb
)
{
    const auto communicator = UPstream::worldComm;

    // Need to gather data separately since we don't have a corresponding
    // is_contiguous data type...

    typedef Pair<point> bb_type;

    List<bb_type> surfBoxes
    (
        UPstream::listGatherValues<bb_type>(s.box(), communicator)
    );

    List<labelPair> surfSizes
    (
        UPstream::listGatherValues<labelPair>
        (
            labelPair(s.points().size(), s.size()),
            communicator
        )
    );

    if (UPstream::master(communicator))
    {
        for (int proci : UPstream::allProcs(communicator))
        {
            Info<< "processor" << proci << nl;

            const auto [nPoints, nFaces] = surfSizes[proci];
            const auto& [smin, smax] = surfBoxes[proci];
            const auto& bbs = meshBb[proci];

            forAll(bbs, i)
            {
                if (!i)
                {
                    Info<< "\tMesh bounds          : ";
                }
                else
                {
                    Info<< "\t                       ";
                }
                Info<< bbs[i] << nl;
            }
            Info<< "\tSurface bounding box : " << smin << ' ' << smax << nl
                << "\tTriangles            : " << nFaces << nl
                << "\tVertices             : " << nPoints << endl;
        }
        Info<< endl;
    }
}


// Extract and return the surface bounding boxes
List<treeBoundBox> getSurfaceBoundingBoxes(const triSurface& s)
{
    const auto communicator = UPstream::worldComm;

    treeBoundBox localBb(s.box());

    return UPstream::listGatherValues(localBb, communicator);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Redistribute a triSurface."
        " The specified surface must be located in the constant/triSurface"
        " directory"
    );

    argList::addArgument("triSurfaceMesh");
    argList::addArgument("distributionType");
    argList::addBoolOption
    (
        "keepNonMapped",
        "Preserve surface outside of mesh bounds"
    );
    argList::addDryRunOption
    (
        "Test without writing the decomposition."
        " Writes surfaces and bounding boxes as VTK files."
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();  // Extra safety?

    const auto surfFileName = args.get<fileName>(1);
    const auto distTypeName = args.get<word>(2);
    const label distType =
        distributedTriSurfaceMesh::distributionTypeNames_[distTypeName];

    Info<< "Reading surface from " << surfFileName << nl << nl
        << "Using distribution method "
        << distTypeName << nl << endl;

    const bool keepNonMapped = args.found("keepNonMapped");

    if (keepNonMapped)
    {
        Info<< "Preserving surface outside of mesh bounds." << nl << endl;
    }
    else
    {
        Info<< "Removing surface outside of mesh bounds." << nl << endl;
    }


    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Please run this program on the decomposed case."
            << " It will read surface " << surfFileName
            << " and decompose it such that it overlaps the mesh bounding box."
            << exit(FatalError);
    }


    Random rndGen(653213);

    // For independent decomposition, ensure that distributedTriSurfaceMesh
    // can find the alternative decomposeParDict specified via the
    // -decomposeParDict option.
    if (distType == distributedTriSurfaceMesh::INDEPENDENT)
    {
        // Ensure demand-driven decompositionMethod finds alternative
        // decomposeParDict location properly.

        IOdictionary* dictPtr = new IOdictionary
        (
            IOobject::selectIO
            (
                IOobject
                (
                    decompositionModel::canonicalName,
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                args.getOrDefault<fileName>("decomposeParDict", "")
            )
        );

        // Store it on the object registry, but to be found it must also
        // have the expected "decomposeParDict" name.

        dictPtr->rename(decompositionModel::canonicalName);
        runTime.store(dictPtr);
    }

    // Determine mesh bounding boxes:
    List<List<treeBoundBox>> meshBb(Pstream::nProcs());
    if (distType == distributedTriSurfaceMesh::FOLLOW)
    {
        #include "createPolyMesh.H"

        meshBb[UPstream::myProcNo()].push_back
        (
            treeBoundBox(mesh.points()).extend(rndGen, 1e-3)
        );
        Pstream::allGatherList(meshBb);
    }

    IOobject io
    (
        surfFileName,         // name
        //runTime.findInstance("triSurface", surfFileName),   // instance
        runTime.constant(),   // instance
        "triSurface",         // local
        runTime,              // registry
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    // Look for file (using searchableSurface rules)
    const fileName actualPath(io.typeFilePath<searchableSurface>());
    fileName localPath(actualPath);
    localPath.replace(runTime.rootPath() + '/', "");


    autoPtr<distributedTriSurfaceMesh> surfMeshPtr;

    if (actualPath == io.objectPath())
    {
        Info<< "Loading local (decomposed) surface " << localPath << nl <<endl;
        surfMeshPtr.reset(new distributedTriSurfaceMesh(io));
    }
    else
    {
        Info<< "Loading undecomposed surface " << localPath
            << " on master only" << endl;

        triSurface s;
        List<treeBoundBox> bbs;
        if (Pstream::master())
        {
            // Actually load the surface
            const bool oldParRun = Pstream::parRun(false);
            triSurfaceMesh surf(io);
            Pstream::parRun(oldParRun);  // Restore parallel state
            s = surf;
            bbs.push_back(treeBoundBox(boundBox::greatBox));
        }
        else
        {
            bbs.push_back(treeBoundBox::null());
        }

        dictionary dict;
        dict.add("distributionType", distTypeName);
        dict.add("mergeDistance", SMALL);
        dict.add("bounds", bbs);

        // Scatter patch information
        Pstream::broadcast(s.patches());

        // Construct distributedTrisurfaceMesh from components
        IOobject notReadIO(io);
        notReadIO.readOpt(IOobject::NO_READ);
        surfMeshPtr.reset(new distributedTriSurfaceMesh(notReadIO, s, dict));
    }

    distributedTriSurfaceMesh& surfMesh = surfMeshPtr();


    // Write per-processor stats
    Info<< "Before redistribution:" << endl;
    writeProcStats(surfMesh, meshBb);


    // Do redistribution
    Info<< "Redistributing surface" << nl << endl;
    autoPtr<mapDistribute> faceMap;
    autoPtr<mapDistribute> pointMap;
    surfMesh.distribute
    (
        meshBb[Pstream::myProcNo()],
        keepNonMapped,
        faceMap,
        pointMap
    );
    faceMap.clear();
    pointMap.clear();

    Info<< endl;


    // Write per-processor stats
    Info<< "After redistribution:" << endl;
    writeProcStats(surfMesh, meshBb);

    if (args.dryRun())
    {
        // Write surface
        {
            const auto& points = surfMesh.points();
            const auto& faces  = surfMesh.surfFaces();

            word outputName(args.executable() + ".vtp");

            Info<< "\ndry-run: writing " << outputName << nl;

            vtk::surfaceWriter writer
            (
                points,
                faces,
                (runTime.globalPath()/outputName),
                UPstream::parRun()
            );

            writer.writeGeometry();

            writer.beginCellData(1 + (writer.parallel() ? 1 : 0));

            writer.writeProcIDs();

            labelList regionIDs(surfMesh.regionIDs());

            writer.writeCellData("regionID", regionIDs);
        }

        // Write surface bounding boxes
        {
            List<treeBoundBox> surfBbs = getSurfaceBoundingBoxes(surfMesh);

            pointField bbCorners
            (
                surfBbs.size()*treeBoundBox::nPoints(),
                point::zero
            );

            List<cellShape> shapes
            (
                surfBbs.size(),
                cellShape(cellModel::HEX, Foam::identity(8))
            );

            {
                forAll(surfBbs, boxi)
                {
                    const auto& bb = surfBbs[boxi];
                    auto& shape = shapes[boxi];

                    const label beg = boxi*treeBoundBox::nPoints();
                    const label len = treeBoundBox::nPoints();

                    Foam::identity(static_cast<labelList&>(shape), beg);

                    if (bb.good())
                    {
                        auto corners = bbCorners.slice(beg, len);
                        bb.hexCorners(corners);
                    }
                }
            }

            word outputName(args.executable() + "_bbs.vtu");

            Info<< "\ndry-run: writing " << outputName << nl;

            if (UPstream::master())
            {
                vtk::vtuCells bbCells(shapes);

                // ASCII format for easier inspection
                vtk::vtuWriter writer
                (
                    bbCorners,
                    bbCells,
                    vtk::formatType::INLINE_ASCII,
                    (runTime.globalPath()/outputName),
                    false  // not parallel
                );

                writer.writeGeometry();

                writer.beginCellData(1);

                labelList procIDs(Foam::identity(shapes.size()));

                writer.writeCellData("procID", procIDs);
            }
        }

        Info<< nl;
    }
    else
    {
        Info<< "Writing surface." << nl << endl;
        surfMesh.objectRegistry::write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
