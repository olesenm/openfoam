/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2026 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "faMeshTools.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelList Foam::faMeshTools::edgeLabels(const faMesh& mesh)
{
    labelList values(Foam::identity(mesh.nEdges()));

    // Internal edges are in the same order for faMesh and primitive patch
    label start = mesh.nInternalEdges();

    // Placeholder to detect if anything was missed
    values.slice(start) = -1;

    // Boundaries
    for (const auto& p : mesh.boundary())
    {
        const auto& labels = p.edgeLabels();
        const auto count = labels.size();

        values.slice(start, count) = labels;
        start += count;
    }

    return values;
}


void Foam::faMeshTools::printMeshChecks
(
    const faMesh& mesh,
    const int verbose
)
{
    const faBoundaryMesh& patches = mesh.boundary();
    const label nNonProcessor = patches.nNonProcessor();
    const label nPatches = patches.size();

    label nLocalProcEdges = 0;
    if (UPstream::parRun())
    {
        for (const faPatch& fap : patches)
        {
            const auto* cpp = isA<processorFaPatch>(fap);

            if (cpp)
            {
                nLocalProcEdges += fap.nEdges();
            }
        }
    }

    const labelList nFaces
    (
        UPstream::listGatherValues<label>(mesh.nFaces())
    );
    const labelList nPoints
    (
        UPstream::listGatherValues<label>(mesh.nPoints())
    );
    const labelList nEdges
    (
        UPstream::listGatherValues<label>(mesh.nEdges())
    );
    const labelList nIntEdges
    (
        UPstream::listGatherValues<label>(mesh.nInternalEdges())
    );

    // The "real" (non-processor) boundary edges
    const labelList nBndEdges
    (
        UPstream::listGatherValues<label>
        (
            mesh.nBoundaryEdges() - nLocalProcEdges
        )
    );
    const labelList nProcEdges
    (
        UPstream::listGatherValues<label>(nLocalProcEdges)
    );


    // Format output as
    //  Number of faces: ...
    //         per-proc: (...)

    const auto reporter =
        [&,verbose](const char* tag, const labelList& list)
        {
            Info<< "  Number of " << tag << ": " << sum(list) << nl;
            if (UPstream::parRun() && verbose)
            {
                int padding = static_cast<int>
                (
                    // strlen("  Number of ") - strlen("per-proc")
                    (12 - 8)
                  + strlen(tag)
                );

                do { Info<< ' '; } while (--padding > 0);

                Info<< "per-proc: " << flatOutput(list) << nl;
            }
        };


    Info<< "----------------" << nl
        << "Mesh Information" << nl
        << "----------------" << nl
        << "  " << "boundingBox: "
        << boundBox(mesh.points(), UPstream::parRun()) << nl;

    if (UPstream::master())
    {
        reporter("faces", nFaces);
        reporter("points", nPoints);
        reporter("edges", nEdges);
        reporter("internal edges", nIntEdges);
        reporter("boundary edges", nBndEdges);

        if (UPstream::parRun())
        {
            reporter("processor edges", nProcEdges);
        }
    }

    Info<< "----------------" << nl
        << "Patches" << nl
        << "----------------" << nl;

    for (label patchi = 0; patchi < nNonProcessor; ++patchi)
    {
        const faPatch& p = patches[patchi];

        // Report physical size (nEdges) not virtual size
        Info<< "  " << "patch " << p.index()
            << " (size: " << returnReduce(p.nEdges(), sumOp<label>())
            << ") name: " << p.name()
            << nl;
    }

    Info<< "----------------" << nl
        << "Used polyPatches: " << flatOutput(mesh.whichPolyPatches()) << nl;


    // Geometry information
    Info<< nl;
    {
        scalarMinMax limit(gMinMax(mesh.S().field()));
        Info<< "Face area:" << nl
            << "    min = " << limit.min() << " max = "  << limit.max() << nl;
    }

    {
        scalarMinMax limit(minMax(mesh.magLe().primitiveField()));

        // Include processor boundaries into 'internal' edges
        if (UPstream::parRun())
        {
            for (label patchi = nNonProcessor; patchi < nPatches; ++patchi)
            {
                limit.add(minMax(mesh.magLe().boundaryField()[patchi]));
            }

            reduce(limit, plusOp<scalarMinMax>{});
        }

        Info<< "Edge length (internal):" << nl
            << "    min = " << limit.min() << " max = "  << limit.max() << nl;


        // Include (non-processor) boundaries
        for (label patchi = 0; patchi < nNonProcessor; ++patchi)
        {
            limit.add(minMax(mesh.magLe().boundaryField()[patchi]));
        }

        if (UPstream::parRun())
        {
            reduce(limit, plusOp<scalarMinMax>{});
        }

        Info<< "Edge length:" << nl
            << "    min = " << limit.min() << " max = "  << limit.max() << nl;
    }

    // Not particularly meaningful
    #if 0
    {
        MinMax<vector> limit(gMinMax(mesh.faceAreaNormals().field()));

        Info<< "Face area normals:" << nl
            << "    min = " << limit.min() << " max = " << limit.max() << nl;
    }
    #endif
}


// ************************************************************************* //
