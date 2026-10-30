/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2026 OpenCFD Ltd.
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

#include "foamVtkCoordSetWriter.H"
#include "foamVtkOutput.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::coordSetWriter::beginPiece()
{
    // Update sizes, similar to
    // vtk::polyWriter::beginPiece(const pointField&, const edgeList&)

    // Basic sizes
    label nPoints = 0;
    label nLines = 0;

    for (const pointField& pts : points_)
    {
        if (auto npts = pts.size(); npts > 0)
        {
            nPoints += npts;
            ++nLines;
        }
    }

    pointSlab_ = nPoints;
    cellSlab_ = nLines;

    switch (elemOutput_)
    {
        case elemOutputType::NO_ELEMENTS:
        {
            cellSlab_ = 0;
            break;
        }
        case elemOutputType::DEFAULT_ELEMENTS:
        {
            if (points_.size() < 2)
            {
                cellSlab_ = 0;
            }
            break;
        }
        case elemOutputType::VERTEX_ELEMENTS:
        {
            // Vertex cells instead of lines
            cellSlab_ = pointSlab_;
            break;
        }
        case elemOutputType::LINE_ELEMENTS:
        {
            // Already determined
            break;
        }
    }

    // if (parallel_)
    // {
    //     Foam::reduceOffsets
    //     (
    //         UPstream::worldComm,
    //         pointSlab_,
    //         cellSlab_
    //     );
    // }


    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        format().openTag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, nTotalPoints()
        );

        if (nTotalCells())
        {
            if (elemOutput_ == elemOutputType::VERTEX_ELEMENTS)
            {
                format().xmlAttr(vtk::fileAttr::NUMBER_OF_VERTS, nTotalCells());
            }
            else
            {
                format().xmlAttr(vtk::fileAttr::NUMBER_OF_LINES, nTotalCells());
            }
        }
        format().closeTag();
    }
}


void Foam::vtk::coordSetWriter::writePoints()
{
    this->beginPoints(nTotalPoints());

    {
        for (const pointField& pts : points_)
        {
            vtk::writeList(format(), pts);
        }
    }

    this->endPoints();
}


void Foam::vtk::coordSetWriter::writeLines_legacy()
{
    if
    (
        (elemOutput_ == elemOutputType::VERTEX_ELEMENTS)
     || (cellSlab_.total() == 0)
    )
    {
        return;  // Nothing to do
    }

    // connectivity = use each point
    const label nLocalLines = cellSlab_.size();
    const label nLocalConns = pointSlab_.size();

    legacy::beginLines(os_, nLocalLines, nLocalConns);

    labelList vertLabels(nLocalLines + nLocalConns);

    auto iter = vertLabels.begin();

    label localPointi = 0;
    for (const pointField& pts : points_)
    {
        if (label npts = pts.size(); npts > 0)
        {
            *iter++ = npts;
            while (npts--)
            {
                *iter++ = localPointi;
                ++localPointi;
            }
        }
    }

    vtk::writeList(format(), vertLabels);

    if (format_)
    {
        format().flush();
    }
}


void Foam::vtk::coordSetWriter::writeLines()
{
    if
    (
        (elemOutput_ == elemOutputType::VERTEX_ELEMENTS)
     || (cellSlab_.total() == 0)
    )
    {
        return;  // Nothing to do
    }

    // connectivity = use each point
    const label nLocalLines = cellSlab_.size();
    const label nLocalConns = pointSlab_.size();

    if (format_)
    {
        format().tag(vtk::fileTag::LINES);
    }

    //
    // 'offsets'  (connectivity offsets)
    //
    {
        labelList vertOffsets(nLocalLines);
        label nOffs = vertOffsets.size();

        // if (parallel_)
        // {
        //     reduce(nOffs, sumOp<label>());
        // }

        if (format_)
        {
            auto payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }

        // processor-local connectivity offsets
        label off = 0;

        /// label off =
        /// (
        ///     parallel_ ? UPstream::find_offset(nLocalConns) : 0
        /// );

        auto iter = vertOffsets.begin();

        for (const pointField& pts : points_)
        {
            if (auto npts = pts.size(); npts > 0)
            {
                off += npts;  // End offset
                *iter = off;
                ++iter;
            }
        }

        vtk::writeList(format_.ref(), vertOffsets);

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalConns);

        label nConns = nLocalConns;

        // if (parallel_)
        // {
        //     reduce(nConns, sumOp<label>());
        // }

        if (format_)
        {
            auto payLoad = vtk::sizeofData<label>(nConns);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);
        }

        {
            // XML: connectivity only
            // [id1, id2, ..., id1, id2, ...]

            auto iter = vertLabels.begin();

            label localPointi = 0;
            for (const pointField& pts : points_)
            {
                label npts = pts.size();

                while (npts--)
                {
                    *iter++ = localPointi;
                    ++localPointi;
                }
            }
        }


        vtk::writeList(format(), vertLabels);

        if (format_)
        {
            format().flush();
            format().endDataArray();
        }
    }

    if (format_)
    {
        format().endTag(vtk::fileTag::LINES);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::coordSetWriter::coordSetWriter
(
    const UPtrList<const pointField>& points,
    const vtk::outputOptions opts
)
:
    vtk::polyWriter(opts),

    points_(points),
    instant_(),
    elemOutput_(elemOutputType::DEFAULT_ELEMENTS)
{}


Foam::vtk::coordSetWriter::coordSetWriter
(
    const UPtrList<const pointField>& points,
    const fileName& file,
    bool parallel
)
:
    coordSetWriter(points)
{
    open(file, parallel);
}


Foam::vtk::coordSetWriter::coordSetWriter
(
    const UPtrList<const pointField>& points,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    coordSetWriter(points, opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::coordSetWriter::setElementType
(
    const elemOutputType elemOutput
)
{
    elemOutput_ = elemOutput;
}


bool Foam::vtk::coordSetWriter::open
(
    const fileName& file,
    bool parallel
)
{
    return vtk::polyWriter::open(file, false);  // non-parallel only
}


void Foam::vtk::coordSetWriter::setTime(const instant& inst)
{
    instant_ = inst;
}


bool Foam::vtk::coordSetWriter::beginFile(std::string title)
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
    return vtk::fileWriter::beginFile("coord-set");
}


bool Foam::vtk::coordSetWriter::writeGeometry()
{
    enter_Piece();

    beginPiece();

    writePoints();

    if (elemOutput_ == elemOutputType::VERTEX_ELEMENTS)
    {
        if (label nTotalVerts = cellSlab_.total(); nTotalVerts > 0)
        {
            vtk::polyWriter::writeVerts(nTotalVerts);
        }
    }
    else
    {
        if (legacy())
        {
            writeLines_legacy();
        }
        else
        {
            writeLines();
        }
    }

    return true;
}


void Foam::vtk::coordSetWriter::writeTimeValue()
{
    if (!instant_.name().empty())
    {
        vtk::fileWriter::writeTimeValue(instant_.value());
    }
}


void Foam::vtk::coordSetWriter::piece
(
    const UPtrList<const pointField>& points
)
{
    endPiece();

    points_ = points;
}


bool Foam::vtk::coordSetWriter::writeProcIDs()
{
    // Ignore
    return false;
}


// ************************************************************************* //
