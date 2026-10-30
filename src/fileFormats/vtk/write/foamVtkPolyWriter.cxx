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

#include "foamVtkPolyWriter.H"
#include "foamVtkOutput.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// The connectivity count for a list of faces (or edges)
template<class ElemType>
static Foam::label countConnectivity(const Foam::UList<ElemType>& elements)
{
    using namespace Foam;

    if constexpr (std::is_same_v<Foam::edge, ElemType>)
    {
        // Edges have exactly two ends each
        return 2 * elements.size();
    }
    else if constexpr
    (
        std::is_same_v<Foam::triFace, ElemType>
     || std::is_same_v<Foam::labelledTri, ElemType>
    )
    {
        // Triangles have exactly three points each
        return 3 * elements.size();
    }
    else
    {
        // General element type
        label nConnectivity = 0;

        for (const auto& f : elements)
        {
            nConnectivity += f.size();
        }

        return nConnectivity;
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::polyWriter::beginPiece
(
    const pointField& points,
    const UList<edge>& edges
)
{
    // Basic sizes
    pointSlab_ = points.size();
    cellSlab_ = edges.size();

    if (parallel_)
    {
        Foam::reduceOffsets
        (
            UPstream::worldComm,
            pointSlab_,
            cellSlab_
        );
    }


    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        format().tag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, nTotalPoints(),
            vtk::fileAttr::NUMBER_OF_LINES,  nTotalCells()
            // AND: vtk::fileAttr::NUMBER_OF_POLYS,  0
        );
    }
}


template<class FaceType>
void Foam::vtk::polyWriter::beginPiece_impl
(
    const pointField& points,
    const UList<FaceType>& faces
)
{
    // Basic sizes
    pointSlab_ = points.size();
    cellSlab_ = faces.size();

    if (parallel_)
    {
        Foam::reduceOffsets
        (
            UPstream::worldComm,
            pointSlab_,
            cellSlab_
        );
    }


    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        format().tag
        (
            vtk::fileTag::PIECE,
            vtk::fileAttr::NUMBER_OF_POINTS, nTotalPoints(),
            vtk::fileAttr::NUMBER_OF_POLYS,  nTotalCells()
            // AND: vtk::fileAttr::NUMBER_OF_LINES,  0
        );
    }
}


void Foam::vtk::polyWriter::beginPiece
(
    const pointField& points,
    const bool useVerts
)
{
    // Basic sizes
    pointSlab_ = points.size();
    cellSlab_ = points.size();  // One vertex per point

    if (parallel_)
    {
        Foam::reduceOffsets
        (
            UPstream::worldComm,
            pointSlab_,
            cellSlab_
        );
    }

    if (!useVerts)
    {
        cellSlab_ = 0;
    }

    // Nothing else to do for legacy
    if (legacy()) return;

    if (format_)
    {
        if (useVerts)
        {
            format().tag
            (
                vtk::fileTag::PIECE,
                vtk::fileAttr::NUMBER_OF_POINTS, nTotalPoints(),
                vtk::fileAttr::NUMBER_OF_VERTS,  nTotalCells()
            );
        }
        else
        {
            format().tag
            (
                vtk::fileTag::PIECE,
                vtk::fileAttr::NUMBER_OF_POINTS, nTotalPoints()
            );
        }
    }
}


void Foam::vtk::polyWriter::writePoints
(
    const pointField& points
)
{
    this->beginPoints(nTotalPoints());

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), points);
    }
    else
    {
        vtk::writeList(format(), points);

    }

    this->endPoints();
}


void Foam::vtk::polyWriter::writeLines_legacy
(
    const UList<edge>& edges
)
{
    // The processor-local point offset
    const label pointOffset = pointSlab_.start();

    // Connectivity count without additional storage (done internally)
    const label nLocalLines = edges.size();
    const label nLocalConns = countConnectivity(edges);

    label nLines = nLocalLines;
    label nConns = nLocalConns;

    if (parallel_)
    {
        reduce(nLines, sumOp<label>());
        reduce(nConns, sumOp<label>());
    }

    if (nLines != nTotalCells())
    {
        FatalErrorInFunction
            << "Expecting " << nTotalCells()
            << " edges, but found " << nLines
            << exit(FatalError);
    }

    legacy::beginLines(os_, nLines, nConns);

    labelList vertLabels(nLocalLines + nLocalConns);

    {
        // Legacy: size + connectivity together
        // [nPts, id1, id2, ..., nPts, id1, id2, ...]

        auto iter = vertLabels.begin();

        const label off = pointOffset;

        for (const auto& e : edges)
        {
            *iter = e.size();   // The size prefix (always 2 for an edge)
            ++iter;

            *iter = off + e.first();    // Vertex labels
            ++iter;

            *iter = off + e.second();
            ++iter;
        }
    }


    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), vertLabels);
    }
    else
    {
        vtk::writeList(format(), vertLabels);
    }

    if (format_)
    {
        format().flush();
    }
}


void Foam::vtk::polyWriter::writeLines
(
    const UList<edge>& edges
)
{
    // The processor-local point offset
    const label pointOffset = pointSlab_.start();

    // Connectivity count without additional storage (done internally)
    const label nLocalLines = edges.size();
    const label nLocalConns = countConnectivity(edges);

    if (format_)
    {
        format().tag(vtk::fileTag::LINES);
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalConns);

        label nConns = nLocalConns;

        if (parallel_)
        {
            reduce(nConns, sumOp<label>());
        }

        if (format_)
        {
            const auto payLoad = vtk::sizeofData<label>(nConns);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);
        }

        {
            // XML: connectivity only
            // [id1, id2, ..., id1, id2, ...]

            auto iter = vertLabels.begin();

            const label off = pointOffset;

            for (const auto& e : edges)
            {
                // Edge vertex labels
                *iter = off + e.first();
                ++iter;

                *iter = off + e.second();
                ++iter;
            }
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertLabels);
        }
        else
        {
            vtk::writeList(format(), vertLabels);
        }

        this->endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        labelList vertOffsets(nLocalLines);
        label nOffs = vertOffsets.size();

        if (parallel_)
        {
            reduce(nOffs, sumOp<label>());
        }

        if (format_)
        {
            const auto payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }


        // processor-local connectivity offsets
        label off =
        (
            parallel_ ? UPstream::find_offset(nLocalConns) : 0
        );


        auto iter = vertOffsets.begin();

        for (const auto& e : edges)
        {
            off += e.size();   // End offset
            *iter = off;
            ++iter;
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertOffsets);
        }
        else
        {
            vtk::writeList(format_.ref(), vertOffsets);
        }


        this->endDataArray();
    }

    if (format_)
    {
        format().endTag(vtk::fileTag::LINES);
    }
}


template<class FaceType>
void Foam::vtk::polyWriter::writePolys_legacy_impl
(
    const UList<FaceType>& faces
)
{
    // The processor-local point offset
    const label pointOffset = pointSlab_.start();

    // Connectivity count without additional storage (done internally)
    const label nLocalPolys = faces.size();
    const label nLocalConns = countConnectivity(faces);

    label nPolys = nLocalPolys;
    label nConns = nLocalConns;

    if (parallel_)
    {
        reduce(nPolys, sumOp<label>());
        reduce(nConns, sumOp<label>());
    }

    if (nPolys != nTotalCells())
    {
        FatalErrorInFunction
            << "Expecting " << nTotalCells()
            << " faces, but found " << nPolys
            << exit(FatalError);
    }

    legacy::beginPolys(os_, nPolys, nConns);

    labelList vertLabels(nLocalPolys + nLocalConns);

    {
        // Legacy: size + connectivity together
        // [nPts, id1, id2, ..., nPts, id1, id2, ...]

        auto iter = vertLabels.begin();

        const label off = pointOffset;

        for (const auto& f : faces)
        {
            *iter = f.size();       // The size prefix
            ++iter;

            for (const label id : f)
            {
                *iter = id + off;   // Vertex label
                ++iter;
            }
        }
    }


    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), vertLabels);
    }
    else
    {
        vtk::writeList(format(), vertLabels);
    }

    if (format_)
    {
        format().flush();
    }
}


template<class FaceType>
void Foam::vtk::polyWriter::writePolys_impl
(
    const UList<FaceType>& faces
)
{
    // The processor-local point offset
    const label pointOffset = pointSlab_.start();

    // Connectivity count without additional storage (done internally)
    const label nLocalPolys = faces.size();
    const label nLocalConns = countConnectivity(faces);

    if (format_)
    {
        format().tag(vtk::fileTag::POLYS);
    }

    //
    // 'connectivity'
    //
    {
        labelList vertLabels(nLocalConns);

        label nConns = nLocalConns;

        if (parallel_)
        {
            reduce(nConns, sumOp<label>());
        }

        if (format_)
        {
            const auto payLoad = vtk::sizeofData<label>(nConns);

            format().beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
            format().writeSize(payLoad);
        }

        {
            // XML: connectivity only
            // [id1, id2, ..., id1, id2, ...]

            auto iter = vertLabels.begin();

            label off = pointOffset;

            for (const auto& f : faces)
            {
                for (const label id : f)
                {
                    *iter = id + off;  // Face vertex label
                    ++iter;
                }
            }
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertLabels);
        }
        else
        {
            vtk::writeList(format(), vertLabels);
        }

        this->endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        labelList vertOffsets(nLocalPolys);
        label nOffs = vertOffsets.size();

        if (parallel_)
        {
            reduce(nOffs, sumOp<label>());
        }

        if (format_)
        {
            const auto payLoad = vtk::sizeofData<label>(nOffs);

            format().beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
            format().writeSize(payLoad);
        }


        // processor-local connectivity offsets
        label off =
        (
            parallel_ ? UPstream::find_offset(nLocalConns) : 0
        );


        auto iter = vertOffsets.begin();

        for (const auto& f : faces)
        {
            off += f.size();  // End offset
            *iter = off;
            ++iter;
        }


        if (parallel_)
        {
            vtk::writeListParallel(format_.ref(), vertOffsets);
        }
        else
        {
            vtk::writeList(format_.ref(), vertOffsets);
        }


        this->endDataArray();
    }

    if (format_)
    {
        format().endTag(vtk::fileTag::POLYS);
    }
}


void Foam::vtk::polyWriter::writeVerts
(
    const label nTotalVerts
)
{
    // Note: nTotalVerts is usually identical to cellSlab_.total()
    // and should be identical to pointSlab_.total() as well

    if (legacy())
    {
        legacy::beginVerts(os_, nTotalVerts);

        // Legacy: size + connectivity together
        // [1, id1, 1, id2, ..., ...]

        // Have enough information to write on master only
        if (format_)
        {
            auto& fmt = format();

            // connectivity = 1 per vertex
            const label connect(1);

            for (label verti = 0; verti < nTotalVerts; ++verti)
            {
                vtk::write(fmt, connect);  // The size prefix
                vtk::write(fmt, verti);    // Vertex label
            }

            fmt.flush();
        }
    }
    else
    {
        // Same payload for connectivity and offsets
        const auto payLoad = vtk::sizeofData<label>(nTotalVerts);

        // Have enough information to write on master only
        if (format_)
        {
            auto& fmt = format();

            fmt.tag(vtk::fileTag::VERTS);

            // 'connectivity' = linear mapping onto points
            {
                fmt.beginDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY);
                fmt.writeSize(payLoad);

                vtk::writeIdentity(fmt, nTotalVerts);

                fmt.flush();
                fmt.endDataArray();
            }

            // 'offsets' (connectivity end offsets)
            // = linear mapping onto points (with 1 offset)
            {
                fmt.beginDataArray<label>(vtk::dataArrayAttr::OFFSETS);
                fmt.writeSize(payLoad);

                vtk::writeIdentity(fmt, nTotalVerts, 1);

                fmt.flush();
                fmt.endDataArray();
            }

            fmt.endTag(vtk::fileTag::VERTS);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::polyWriter::polyWriter
(
    const vtk::outputOptions opts
)
:
    vtk::fileWriter(vtk::fileTag::POLY_DATA, opts)
{
    // We do not currently support append mode
    opts_.append(false);
}


Foam::vtk::polyWriter::polyWriter
(
    const fileName& file,
    bool parallel
)
:
    // Default parameter fails for gcc-4.8.5, thus specify format here
    polyWriter(vtk::formatType::INLINE_BASE64)
{
    open(file, parallel);
}


Foam::vtk::polyWriter::polyWriter
(
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
:
    polyWriter(opts)
{
    open(file, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtk::polyWriter::writeGeometry()
{
    FatalErrorInFunction
        << "Method was not overloaded, called without a geometry!!" << nl
        << "    Indicates a programming error" << nl << endl
        << abort(FatalError);

    return false;
}


bool Foam::vtk::polyWriter::writeLineGeometry
(
    const pointField& points,
    const UList<edge>& edges
)
{
    enter_Piece();

    beginPiece(points, edges);

    writePoints(points);

    if (legacy())
    {
        writeLines_legacy(edges);
    }
    else
    {
        writeLines(edges);
    }

    return true;
}


template<class FaceType>
bool Foam::vtk::polyWriter::writePolyGeometry_impl
(
    const pointField& points,
    const UList<FaceType>& faces
)
{
    enter_Piece();

    beginPiece_impl(points, faces);

    writePoints(points);

    if (legacy())
    {
        writePolys_legacy_impl(faces);
    }
    else
    {
        writePolys_impl(faces);
    }

    return true;
}


bool Foam::vtk::polyWriter::writePolyGeometry
(
    const pointField& points,
    const UList<face>& faces
)
{
    return writePolyGeometry_impl(points, faces);
}


bool Foam::vtk::polyWriter::writePolyGeometry
(
    const pointField& points,
    const UList<triFace>& faces
)
{
    return writePolyGeometry_impl(points, faces);
}


bool Foam::vtk::polyWriter::writePolyGeometry
(
    const pointField& points,
    const UList<labelledTri>& faces
)
{
    return writePolyGeometry_impl(points, faces);
}


bool Foam::vtk::polyWriter::writeVertGeometry
(
    const pointField& points
)
{
    enter_Piece();

    beginPiece(points, true);  //< useVerts = true

    writePoints(points);

    writeVerts(nTotalCells());

    return true;
}


bool Foam::vtk::polyWriter::beginCellData(label nFields)
{
    return enter_CellData(nTotalCells(), nFields);
}


bool Foam::vtk::polyWriter::beginPointData(label nFields)
{
    return enter_PointData(nTotalPoints(), nFields);
}


bool Foam::vtk::polyWriter::writeProcIDs()
{
    return vtk::fileWriter::writeProcIDs
    (
        // Appropriate rank-local size:
        (this->isPointData() ? pointSlab_.size() : cellSlab_.size())
    );
}


void Foam::vtk::polyWriter::writeLocalIDs(const word& fieldName)
{
    vtk::fileWriter::writeLocalIDs
    (
        fieldName,
        // Appropriate rank-local size:
        (this->isPointData() ? pointSlab_.size() : cellSlab_.size())
    );
}


void Foam::vtk::polyWriter::writeGlobalIDs(const word& fieldName)
{
    vtk::fileWriter::writeGlobalIDs
    (
        fieldName,
        // Appropriate rank-local size:
        (this->isPointData() ? pointSlab_.size() : cellSlab_.size())
    );
}


// ************************************************************************* //
