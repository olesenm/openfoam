/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "foamVtkSurfaceWriter.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceType>
void Foam::vtk::surfaceWriter::setFaceReference
(
    const List<FaceType>& faces
)
{
    faces_.reset();
    triFaces_.reset();
    labelledTris_.reset();

    if constexpr (std::is_same_v<Foam::face, FaceType>)
    {
        faces_.cref(faces);
    }
    else if constexpr (std::is_same_v<Foam::triFace, FaceType>)
    {
        triFaces_.cref(faces);
    }
    else if constexpr (std::is_same_v<Foam::labelledTri, FaceType>)
    {
        labelledTris_.cref(faces);
    }
}


template<class FaceType>
Foam::vtk::surfaceWriter::surfaceWriter
(
    const vtk::outputOptions opts,
    const pointField& points,
    const List<FaceType>& faces
)
:
    vtk::polyWriter(opts),

    points_(points),
    instant_(),
    vertexOutput_(false)
{
    setFaceReference(faces);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::surfaceWriter::setTime(const instant& inst)
{
    instant_ = inst;
}


bool Foam::vtk::surfaceWriter::beginFile(std::string title)
{
    if (title.size())
    {
        return vtk::fileWriter::beginFile(title);
    }

    if (!instant_.name().empty())
    {
        return vtk::fileWriter::beginFile
        (
            "time='" + instant_.name() + "'"
        );
    }

    // Provide default title
    return vtk::fileWriter::beginFile("surface");
}


bool Foam::vtk::surfaceWriter::writeGeometry()
{
    if (!points_.good())
    {
        return false;
    }

    if (vertexOutput_)
    {
        return writeVertGeometry(points_.cref());
    }
    else if (const auto* faces_ptr = faces_.get())
    {
        return writePolyGeometry(points_.cref(), *faces_ptr);
    }
    else if (const auto* faces_ptr = triFaces_.get())
    {
        return writePolyGeometry(points_.cref(), *faces_ptr);
    }
    else if (const auto* faces_ptr = labelledTris_.get())
    {
        return writePolyGeometry(points_.cref(), *faces_ptr);
    }
    else
    {
        return false;
    }
}


void Foam::vtk::surfaceWriter::writeTimeValue()
{
    if (!instant_.name().empty())
    {
        vtk::fileWriter::writeTimeValue(instant_.value());
    }
}


// * * * * * * * * * * * * * * * Concrete Types  * * * * * * * * * * * * * * //

#undef  define_methods
#define define_methods(FaceType)                                              \
                                                                              \
Foam::vtk::surfaceWriter::surfaceWriter                                       \
(                                                                             \
    const pointField& points,                                                 \
    const List<FaceType>& faces,                                              \
    const vtk::outputOptions opts                                             \
)                                                                             \
:                                                                             \
    surfaceWriter(opts, points, faces)                                        \
{}                                                                            \
                                                                              \
void Foam::vtk::surfaceWriter::piece                                          \
(                                                                             \
    const pointField& points,                                                 \
    const List<FaceType>& faces                                               \
)                                                                             \
{                                                                             \
    endPiece();                                                               \
    points_.cref(points);                                                     \
    setFaceReference(faces);                                                  \
}


define_methods(Foam::face);
define_methods(Foam::triFace);
define_methods(Foam::labelledTri);

#undef define_methods


// ************************************************************************* //
