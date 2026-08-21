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
    Test-triSurfaceDecompose.cxx

Description
    Test tree decompose for triSurface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "cellShape.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "treeBoundBox.H"
#include "treeDecomp.H"
#include "Time.H"
#include "foamVtkSurfaceWriter.H"
#include "foamVtuWriter.H"

using namespace Foam;

static word getExtension(const fileName& name)
{
    return
    (
        name.has_ext("gz")
      ? name.stem().ext()
      : name.ext()
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Test tree decomposer for triSurface");

    argList::addOption
    (
        "domains",
        "int",
        "number of domains (default: 8)"
    );
    argList::addOption
    (
        "stride",
        "int",
        "stride for agglomerating on domains"
    );
    argList::addBoolOption
    (
        "bitset",
        "use bitset for stride selection"
    );
    argList::addBoolOption
    (
        "proc-write",
        "write individual surfaces - per output processor"
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");

    argList::addVerboseOption();

    argList args(argc, argv);

    const int verbose = args.verbose();
    const bool procWrite = args.found("proc-write");

    const auto importName = args.get<fileName>(1);

    const word readFileType
    (
        args.getOrDefault<word>("read-format", getExtension(importName))
    );

    label nDomains(8);
    args.readIfPresent("domains", nDomains);

    label domainStride(1);
    args.readIfPresent("stride", domainStride);

    scalar scaleFactor(0);
    Info<< "Reading : " << importName << endl;
    triSurface surf(importName, readFileType, scaleFactor);

    Info<< "Read surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    labelList faceToProc;

    Info<< "Tree decompose into " << nDomains << " domains";

    if (domainStride < 1 || domainStride >= nDomains)
    {
        Info<< ", invalid stride=" << domainStride << " (ignoring)";
    }

    if (domainStride > 1 && domainStride < nDomains)
    {
        Info<< ", stride=" << domainStride << endl;

        autoPtr<treeDecomp> decomposer;

        if (args.found("bitset"))
        {
            bitSet select(nDomains);

            for (label i = 0; i < nDomains; i += domainStride)
            {
                select.set(i);
            }

            decomposer.reset(new treeDecomp(nDomains, select));
        }
        else
        {
            labelList weights(nDomains, Zero);

            for (label i = 0; i < nDomains; i += domainStride)
            {
                // Weight by number of sub-ranks within the stride
                weights[i] = Foam::min(domainStride, nDomains-i);
            }
            decomposer.reset(new treeDecomp(nDomains, weights));
        }

        if (verbose)
        {
            Info<< "processor weights: ";
            decomposer().weights().writeList(Info) << nl;
        }

        faceToProc = decomposer().decompose(surf.faceCentres());
    }
    else
    {
        Info<< endl;

        treeDecomp decomposer(nDomains);

        faceToProc = decomposer.decompose(surf.faceCentres());
    }


    // Write surface
    {
        const auto& allPoints = surf.points();
        const auto& allFaces  = surf.surfFaces();

        word outputName("triSurfaceDecompose.vtp");

        Info<< "\ndry-run: writing " << outputName << nl;

        vtk::surfaceWriter writer
        (
            allPoints,
            allFaces,
            // (runTime.globalPath()/outputName),
            outputName,
            false   // parallel
        );

        writer.writeGeometry();

        writer.beginCellData(2);

        // writer.writeProcIDs();
        writer.write("procID", faceToProc);
        writer.write("regionID", surf.regionIDs());
    }

    // Bounding box visualization
    List<treeBoundBox> surfBbs(nDomains);

    if (procWrite)
    {
        const auto& allCentres = surf.faceCentres();

        bitSet selectFaces(surf.size());
        treeBoundBox bb_centres, bb_faces;

        for (label proci = 0; proci < nDomains; ++proci)
        {
            selectFaces.reset();
            bb_centres.reset();
            bb_faces.reset();

            forAll(faceToProc, facei)
            {
                if (proci == faceToProc[facei])
                {
                    selectFaces.set(facei);
                    bb_centres.add(allCentres[facei]);
                    bb_faces.add(surf.box(facei));
                }
            }

            triSurface subMesh(surf.subsetMesh(selectFaces));

            const auto [bb_min, bb_max] = subMesh.box();

            surfBbs[proci].reset(bb_min, bb_max);

            const auto& subPoints = subMesh.points();
            const auto& subFaces = subMesh.surfFaces();

            word outputName
            (
                "triSurfaceDecompose_proc" + Foam::name(proci) + ".vtp"
            );

            Info<< "\ndry-run: writing " << outputName << nl
                << "    bb (centres) = " << bb_centres << nl
                << "    bb (faces)   = " << bb_faces << nl;

            vtk::surfaceWriter writer
            (
                subPoints,
                subFaces,
                outputName,
                false   // parallel
            );

            writer.writeGeometry();

            writer.beginCellData(1);
            writer.write("regionID", subMesh.regionIDs());
        }
    }
    else
    {
        for (label proci = 0; proci < nDomains; ++proci)
        {
            auto& bb_faces = surfBbs[proci];

            forAll(faceToProc, facei)
            {
                if (proci == faceToProc[facei])
                {
                    bb_faces.add(surf.box(facei));
                }
            }
        }
    }

    // Write surface bounding boxes
    {
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

        word outputName
        (
            "triSurfaceDecompose_bbs.vtu"
        );

        Info<< "\ndry-run: writing " << outputName << nl;

        {
            vtk::vtuCells bbCells(shapes);

            // ASCII format for easier inspection
            vtk::vtuWriter writer
            (
                bbCorners,
                bbCells,
                vtk::formatType::INLINE_ASCII,
                outputName,
                false  // not parallel
            );

            writer.writeGeometry();

            writer.beginCellData(1);

            labelList procIDs(Foam::identity(shapes.size()));

            writer.writeCellData("procID", procIDs);
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
