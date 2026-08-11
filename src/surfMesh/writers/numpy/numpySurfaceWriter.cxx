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

\*---------------------------------------------------------------------------*/

#include "numpySurfaceWriter.H"
#include "foamNumpyCore.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "primitivePatch.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(numpyWriter);
    addToRunTimeSelectionTable(surfaceWriter, numpyWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, numpyWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::numpyWriter::numpyWriter()
:
    surfaceWriter(),
    writeNormal_(false)
{}


Foam::surfaceWriters::numpyWriter::numpyWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    writeNormal_(options.getOrDefault("normal", false))
{}


Foam::surfaceWriters::numpyWriter::numpyWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    numpyWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::numpyWriter::numpyWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    numpyWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceWriters::numpyWriter::serialWriteGeometry
(
    const fileName& directory,
    const meshedSurf& surf
)
{
    const pointField& points = surf.points();
    const faceList& faces = surf.faces();

    if (verbose_)
    {
        if (this->isPointData())
        {
            Info<< "Writing face points: "
                << directory/"points.npy" << endl;
        }
        else
        {
            Info<< "Writing face centres: "
                << directory/"points.npy" << endl;
        }
    }

    std::unique_ptr<primitivePatch> ppPtr;

    {
        OFstream os
        (
            IOstreamOption::atomicType::ATOMIC,
            directory/"points.npy",
            IOstreamOption::streamFormat::BINARY
        );
        auto& oss = os.stdStream();

        if (this->isPointData() || (this->vertexOutput() && faces.empty()))
        {
            const auto& fld = points;

            fileFormats::numpyCore::writeHeader(oss, fld);
            fileFormats::numpyCore::writeData(oss, fld);
        }
        else
        {
            ppPtr = std::make_unique<primitivePatch>
            (
                SubList<face>(faces),
                points
            );

            const auto& fld = ppPtr->faceCentres();

            fileFormats::numpyCore::writeHeader(oss, fld);
            fileFormats::numpyCore::writeData(oss, fld);
        }
    }

    if (writeNormal_ && bool(ppPtr))
    {
        const auto& pp = *ppPtr;

        OFstream os
        (
            IOstreamOption::atomicType::ATOMIC,
            directory/"normal.npy",
            IOstreamOption::streamFormat::BINARY
        );
        auto& oss = os.stdStream();

        const auto& fld = pp.faceNormals();

        fileFormats::numpyCore::writeHeader(oss, fld);
        fileFormats::numpyCore::writeData(oss, fld);
    }
}


Foam::fileName Foam::surfaceWriters::numpyWriter::write()
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!Foam::isDir(surfaceDir))
        {
            Foam::mkDir(surfaceDir);
        }

        // Write sample locations
        serialWriteGeometry(surfaceDir, surf);
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::numpyWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Field:    rootdir/surfaceName/<TIME>/field.npy

    fileName surfaceDir = outputPath_;
    fileName outputFile(surfaceDir/timeName()/fieldName + ".npy");

    // Implicit geometry merge()
    tmp<Field<Type>> tfield = adjustField(fieldName, mergeField(localValues));

    if (verbose_)
    {
        Info<< " to " << outputFile << endl;
    }


    // const meshedSurf& surf = surface();
    (void) adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!Foam::isDir(outputFile.path()))
        {
            Foam::mkDir(outputFile.path());
        }

        // Write field
        {
            OFstream os
            (
                IOstreamOption::atomicType::ATOMIC,
                outputFile,
                IOstreamOption::streamFormat::BINARY
            );
            auto& oss = os.stdStream();

            const auto& fld = tfield();

            fileFormats::numpyCore::writeHeader(oss, fld);
            fileFormats::numpyCore::writeData(oss, fld);
        }
    }

    return surfaceDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::numpyWriter);


// ************************************************************************* //
