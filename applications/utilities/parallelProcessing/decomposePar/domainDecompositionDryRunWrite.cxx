/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "domainDecompositionDryRun.H"
#include "foamVtkInternalMeshWriter.H"
#include "volFields.H"
#include "OFstream.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::domainDecompositionDryRun::writeVolField
(
    const word& timeName,
    const labelUList& procIds
) const
{
    // Write decomposition as volScalarField for visualization
    volScalarField cellDist
    (
        IOobject
        (
            "cellDist",
            timeName,
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        this->mesh(),
        dimensionedScalar(word::null, dimless, -1),
        fvPatchFieldBase::zeroGradientType()
    );

    forAll(procIds, celli)
    {
        cellDist[celli] = procIds[celli];
    }

    cellDist.correctBoundaryConditions();
    cellDist.write();

    Info<< nl << "Wrote decomposition to "
        << cellDist.objectRelPath()
        << " (volScalarField) for visualization."
        << endl;
}


void Foam::domainDecompositionDryRun::writeVTK
(
    const fileName& file,
    const labelUList& cellToProc
) const
{
    const vtk::vtuCells cells(this->mesh());

    // not parallel
    vtk::internalMeshWriter writer(this->mesh(), cells, file, false);

    writer.writeGeometry();
    writer.beginCellData();
    writer.writeCellData("procID", cellToProc);
    writer.writeCellIDs();

    Info<< "Wrote decomposition to "
        << this->mesh().time().relativePath(writer.output())
        << " for visualization."
        << endl;
}



// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Write mesh cell-cell connectivity as a SCOTCH graph (.grf) file.
//
// File format
//   Line 1  : version (0)
//   Line 2  : numVertices numEdges
//   Line 3  : baseval (0) and flags (hasVertlabels hasEdgeWeights
//             hasVertWeights)
//   Line 4+ : per vertex: degree followed by the (0-based) neighbour list
static void write_grf
(
    const fileName& file,
    const CompactListList<label>& cellCells
)
{
    OFstream os(file + ".grf");

    const label numVertices = cellCells.size();
    const label numEdges = cellCells.totalSize();

    Info<< "    " << os.relativeName() << " : "
        << numVertices << " vertices, "
        << numEdges << " edges (" << (numEdges/2) << ")" << endl;

    // Version 0 = Graph file (.grf)
    os  << "0" << nl;

    // Number of vertices, number of edges (connections)
    os  << numVertices << ' ' << numEdges << nl;

    // Numbering starts from 0
    // 100*hasVertlabels+10*hasEdgeWeights+1*hasVertWeights
    os  << "0 000" << nl;

    for (label vertex = 0; vertex < numVertices; ++vertex)
    {
        const auto neighbours = cellCells[vertex];

        os  << neighbours.size();

        for (label nbr : neighbours)
        {
            os  << ' ' << nbr;
        }
        os  << nl;
    }
}


// Write the points as a SCOTCH geometry (.xyz) file.
// File format:
//   Line 1  : dimension (3) numVertices
//   Line 2+ : per vertex: (0-based) label, followed by x y z coordinates
static void write_xyz
(
    const fileName& file,
    const pointField& points
)
{
    OFstream os(file + ".xyz");
    os.precision(9);

    const label count = points.size();

    Info<< "    " << os.relativeName()
        << " : " << count << " points" << endl;

    // Dimension (3D) and number of vertices
    os  << "3 " << count << nl;

    label id = 0;

    for (const point& p : points)
    {
        os  << id << ' ' << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        ++id;
    }
}


// Write a (trivial, single-domain) SCOTCH mapping (.map) file.
// File format:
//   Line 1  : numVertices
//   Line 2+ : per vertex: (0-based) label, followed by domain number
static void write_map
(
    const fileName& file,
    const labelUList& cellToProc
)
{
    OFstream os(file + ".map");

    const label count = cellToProc.size();

    label maxDoamin = -1;

    if
    (
        auto max = std::max_element(cellToProc.begin(), cellToProc.end());
        max != cellToProc.end()
    )
    {
        maxDoamin = *max;
    }

    Info<< "    " << os.relativeName()
        << " : " << (maxDoamin+1) << " domains" << endl;

    os  << count << nl;

    label id = 0;

    for (label domain : cellToProc)
    {
        os  << id << ' ' << domain << nl;
        ++id;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::domainDecompositionDryRun::writeScotch
(
    const fileName& file,
    const CompactListList<label>& cellCells,
    const labelUList& cellToProc
) const
{
    Info<< "Writing (scotch)" << endl;

    write_grf(file, cellCells);
    write_xyz(file, mesh_.cellCentres());
    write_map(file, cellToProc);

    Info<< endl;
}


// ************************************************************************* //
