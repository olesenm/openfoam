/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

#include "foamVtkInternalMeshWriter.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::vtk::internalMeshWriter::debug = 0;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::internalMeshWriter::internalMeshWriter
(
    const polyMesh& mesh,
    const vtk::vtuCells& cells,
    const vtk::outputOptions opts
)
:
    vtk::vtuWriter(cells, opts),
    mesh_(mesh)
{}


Foam::vtk::internalMeshWriter::internalMeshWriter
(
    const polyMesh& mesh,
    const vtk::vtuCells& cells,
    const fileName& file,
    bool parallel
)
:
    internalMeshWriter(mesh, cells)
{
    open(file, parallel);
}


Foam::vtk::internalMeshWriter::internalMeshWriter
(
    const polyMesh& mesh,
    const vtk::vtuCells& cells,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    internalMeshWriter(mesh, cells, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::internalMeshWriter::beginFile(std::string title)
{
    if (vtk::vtuWriter::debug < debug)
    {
        vtk::vtuWriter::debug = debug;
    }

    const auto& runTime = mesh_.time();

    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }

    // Provide default title

    DebugInFunction
        << "case=" << runTime.caseName()
        << " region=" << mesh_.name()
        << " time=" << runTime.timeName()
        << " index=" << runTime.timeIndex() << endl;

    if (legacy())
    {
        return vtk::fileWriter::beginFile
        (
            runTime.globalCaseName()
        );
    }

    // XML (inline)

    return vtk::fileWriter::beginFile
    (
        "case='" + runTime.globalCaseName()
      + "' region='" + mesh_.name()
      + "' time='" + runTime.timeName()
      + "' index='" + Foam::name(runTime.timeIndex())
      + "'"
    );
}


bool Foam::vtk::internalMeshWriter::writeGeometry()
{
    return writeGeometry_impl(mesh_.points(), mesh_.cellCentres());
}


// ************************************************************************* //
