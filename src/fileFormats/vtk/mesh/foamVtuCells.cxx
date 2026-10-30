/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "globalMeshData.H"
#include "foamVtuCells.H"
#include "foamVtkOutputOptions.H"

#define Foam_vtk_vtuCells_demandDriven_POLY_FACEIDS

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Need special care when concatenating begin/end offsets in parallel.
//
// From the first viable rank: use the begin/end offsets
// Other ranks: use their end offsets only
//
// Returns a labelList since the caller will invariably need a deep copy
// of values anyhow
static Foam::labelList parCoordinatedOffsets
(
    const Foam::labelUList& offsets,
    const int communicator
)
{
    using namespace Foam;

    if (UPstream::is_parallel(communicator))
    {
        // Length of addessed items (as per globalIndex)
        const auto len = (offsets.size() - 1);

        bool isLeader =
        (
            UPstream::find_first(len > 0, communicator)
         == UPstream::myProcNo(communicator)
        );

        if (isLeader)
        {
            // The full list
            return labelList(offsets);
        }
        else if (len > 0)
        {
            // Offsets without the first value
            return labelList(offsets.slice(1));
        }
        else
        {
            // No valid offsets
            return labelList();
        }
    }
    else
    {
        // Return a copy, but unlikely to be called in non-parallel
        return labelList(offsets);
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::vtuCells::vtuCells
(
    const contentType output,
    const bool decompose
)
:
    vtk::vtuSizing(),
    output_(output),
    decomposeRequest_(decompose),
    merged_(vtk::vtuSizing::mergeType::NONE)
{}


Foam::vtk::vtuCells::vtuCells
(
    const polyMesh& mesh,
    const contentType output,
    const bool decompose
)
:
    vtuCells(output, decompose)
{
    reset(mesh);
}


Foam::vtk::vtuCells::vtuCells
(
    const vtk::outputOptions opts,
    const bool decompose
)
:
    vtuCells
    (
        (
            opts.legacy()
          ? contentType::LEGACY
          : opts.is_hdf()
          ? contentType::HDF
          : contentType::XML
        ),
        decompose
    )
{}


Foam::vtk::vtuCells::vtuCells
(
    const polyMesh& mesh,
    const vtk::outputOptions opts,
    const bool decompose
)
:
    vtuCells(opts, decompose)
{
    reset(mesh);
}


Foam::vtk::vtuCells::vtuCells
(
    const UList<cellShape>& shapes,
    const contentType output
)
:
    vtuCells(output, false)  // decompose = false
{
    resetShapes(shapes);
}


Foam::vtk::vtuCells::vtuCells
(
    const UList<cellShape>& shapes,
    const vtk::outputOptions opts
)
:
    vtuCells(opts, false)  // decompose = false
{
    resetShapes(shapes);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::vtuCells::resize_all()
{
    cellTypes_.resize_nocopy(nFieldCells());
    vertLabels_.resize_nocopy(sizeOf<slotType::CELLS>(output_));
    vertOffset_.resize_nocopy(sizeOf<slotType::CELLS_OFFSETS>(output_));
    faceLabels_.resize_nocopy(sizeOf<slotType::FACES>(output_));
    faceOffset_.resize_nocopy(sizeOf<slotType::FACES_OFFSETS>(output_));

    #ifdef Foam_vtk_vtuCells_demandDriven_POLY_FACEIDS
    polyFaceIds_.clear();
    #else
    polyFaceIds_.resize_nocopy(sizeOf<slotType::POLY_FACEIDS>(output_));
    #endif
    polyFaceOffset_.resize_nocopy
    (
        sizeOf<slotType::POLY_FACEIDS_OFFSETS>(output_)
    );
}


void Foam::vtk::vtuCells::populateOutput(const polyMesh& mesh)
{
    // Already called
    // - vtuSizing::reset
    // - resize_all();

    merged_ = vtk::vtuSizing::mergeType::NONE;

    switch (output_)
    {
        case contentType::LEGACY :
        {
            populateLegacy
            (
                mesh,
                cellTypes_,
                vertLabels_,
                maps_
            );
            break;
        }

        case contentType::XML :
        {
            populateXml
            (
                mesh,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                maps_
            );
            break;
        }

        case contentType::INTERNAL2 :
        {
            populateInternal
            (
                mesh,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                maps_,
                output_
            );
            break;
        }

        case contentType::HDF :
        {
            // Note: polyFaceIds_ may be zero-sized (if demand-driven)
            populateHdf
            (
                mesh,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                polyFaceIds_,
                polyFaceOffset_,
                maps_
            );
            break;
        }
    }
}


void Foam::vtk::vtuCells::populateOutput(const UList<cellShape>& shapes)
{
    if (output_ != contentType::LEGACY && output_ != contentType::XML)
    {
        WarningInFunction
            << "Internal formats not supported for shape cells - using XML"
            << nl << nl;

        output_ = contentType::XML;
    }

    vtuSizing::resetShapes(shapes);
    merged_ = vtk::vtuSizing::mergeType::NONE;
    maps_.clear();
    resize_all();
    // Done in populate routine:
    // maps_.cellMap() = identity(vtuSizing::nCells());

    switch (output_)
    {
        case contentType::LEGACY:
        {
            populateShapesLegacy
            (
                shapes,
                cellTypes_,
                vertLabels_,
                maps_
            );
            break;
        }

        case contentType::XML:
        {
            populateShapesXml
            (
                shapes,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                maps_
            );
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unhandled VTK format " << int(output_) << nl
                << exit(FatalError);
            break;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::vtuCells::clear()
{
    vtuSizing::clear();
    merged_ = vtk::vtuSizing::mergeType::NONE;
    cellTypes_.clear();
    vertLabels_.clear();
    vertOffset_.clear();
    faceLabels_.clear();
    faceOffset_.clear();
    polyFaceIds_.clear();
    polyFaceOffset_.clear();

    maps_.clear();
}


void Foam::vtk::vtuCells::reset(const polyMesh& mesh)
{
    vtuSizing::reset(mesh, decomposeRequest_);
    merged_ = vtk::vtuSizing::mergeType::NONE;
    resize_all();

    populateOutput(mesh);
}


void Foam::vtk::vtuCells::reset
(
    const polyMesh& mesh,
    const labelUList& subsetCellsIds
)
{
    vtuSizing::reset(mesh, subsetCellsIds, decomposeRequest_);
    merged_ = vtk::vtuSizing::mergeType::NONE;
    resize_all();

    if (selectionMode() == selectionModeType::SUBSET_MESH)
    {
        maps_.cellMap() = subsetCellsIds;
    }

    populateOutput(mesh);
}


void Foam::vtk::vtuCells::reset
(
    const polyMesh& mesh,
    const enum contentType output,
    const bool decompose
)
{
    output_ = output;
    decomposeRequest_ = decompose;

    reset(mesh);
}


void Foam::vtk::vtuCells::resetShapes
(
    const UList<cellShape>& shapes
)
{
    if (output_ != contentType::LEGACY && output_ != contentType::XML)
    {
        WarningInFunction
            << "VTK internal format is not supported for shape cells"
            << " switching to xml" << nl << nl;

        output_ = contentType::XML;
    }

    decomposeRequest_ = false;
    merged_ = vtk::vtuSizing::mergeType::NONE;

    vtuSizing::resetShapes(shapes);

    maps_.clear();
    resize_all();

    // Create an identity map
    // The size == number of shapes, which is also vtuSizing::nCells()
    maps_.cellMap().resize_nocopy(vtuSizing::nCells());
    Foam::identity(maps_.cellMap());

    switch (output_)
    {
        case contentType::LEGACY:
        {
            populateShapesLegacy
            (
                shapes,
                cellTypes_,
                vertLabels_,
                maps_
            );
            break;
        }

        case contentType::XML:
        {
            populateShapesXml
            (
                shapes,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                maps_
            );
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unhandled VTK format " << int(output_) << nl
                << exit(FatalError);
            break;
        }
    }
}


void Foam::vtk::vtuCells::addPointCellLabels(const labelUList& cellIds)
{
    maps_.additionalIds() = cellIds;
    setNumAddPoints(maps_.additionalIds().size());
}


void Foam::vtk::vtuCells::renumberCells(const labelUList& mapping)
{
    maps_.renumberCells(mapping);
}


void Foam::vtk::vtuCells::renumberPoints(const labelUList& mapping)
{
    maps_.renumberPoints(mapping);
}


const Foam::labelList& Foam::vtk::vtuCells::polyFaceIds() const
{
    #ifdef Foam_vtk_vtuCells_demandDriven_POLY_FACEIDS
    if
    (
        const label numFaces =
        (
            (output_ == contentType::HDF) ? vtuSizing::nFacesPoly() : 0
        );
        (polyFaceIds_.size() < numFaces)
    )
    {
        polyFaceIds_.resize_nocopy(numFaces);
        std::iota(polyFaceIds_.begin(), polyFaceIds_.end(), 0);
    }
    #endif
    return polyFaceIds_;
}


Foam::refPtr<Foam::labelList>
Foam::vtk::vtuCells::vertLabels(label pointOffset) const
{
    if (pointOffset <= 0)
    {
        return refPtr<labelList>(vertLabels_);
    }

    auto tresult = refPtr<labelList>::New(vertLabels_);
    vtuSizing::renumberVertLabels(tresult.ref(), pointOffset, output_);
    return tresult;
}


Foam::refPtr<Foam::labelList>
Foam::vtk::vtuCells::vertOffsets
(
    label beginOffset,
    bool syncPar
) const
{
    // If the format uses begin/end offsets, need special care
    // when concatenating.

    const bool uses_beginEnd_style =
    (
        (output_ == contentType::INTERNAL2)
     || (output_ == contentType::HDF)
    );

    const auto& offsets = vertOffset_;

    const int communicator = UPstream::worldComm;

    if (syncPar && UPstream::parRun() && uses_beginEnd_style)
    {
        auto tresult = refPtr<labelList>::New();
        auto& result = tresult.ref();

        result = parCoordinatedOffsets(offsets, communicator);

        // Determine a consistent beginOffset
        if (beginOffset < 0)
        {
            beginOffset = UPstream::find_offset<label>
            (
                (offsets.size() > 1 ? offsets.back() : 0),
                communicator
            );
        }

        if (beginOffset > 0)
        {
            for (auto& val : result)
            {
                val += beginOffset;
            }
        }

        return tresult;
    }
    else
    {
        if (beginOffset <= 0)
        {
            return refPtr<labelList>(offsets);
        }

        auto tresult = refPtr<labelList>::New(offsets);
        vtuSizing::renumberFaceOffsets(tresult.ref(), beginOffset, output_);
        return tresult;
    }
}


Foam::refPtr<Foam::labelList>
Foam::vtk::vtuCells::faceLabels(label pointOffset) const
{
    if (pointOffset <= 0)
    {
        return refPtr<labelList>(faceLabels_);
    }

    auto tresult = refPtr<labelList>::New(faceLabels_);
    vtuSizing::renumberFaceLabels(tresult.ref(), pointOffset, output_);
    return tresult;
}


Foam::refPtr<Foam::labelList>
Foam::vtk::vtuCells::faceOffsets
(
    label beginOffset,
    bool syncPar
) const
{
    // If the format uses begin/end offsets, need special care
    // when concatenating.

    const bool uses_beginEnd_style =
    (
        (output_ == contentType::HDF)
    );

    const auto& offsets = faceOffset_;

    const int communicator = UPstream::worldComm;

    if (syncPar && UPstream::parRun() && uses_beginEnd_style)
    {
        auto tresult = refPtr<labelList>::New();
        auto& result = tresult.ref();

        result = parCoordinatedOffsets(offsets, communicator);

        // Determine a consistent beginOffset
        if (beginOffset < 0)
        {
            beginOffset = UPstream::find_offset<label>
            (
                (offsets.size() > 1 ? offsets.back() : 0),
                communicator
            );
        }

        if (beginOffset > 0)
        {
            for (auto& val : result)
            {
                val += beginOffset;
            }
        }

        return tresult;
    }
    else if (!offsets.empty())
    {
        if (beginOffset <= 0)
        {
            return refPtr<labelList>(offsets);
        }

        auto tresult = refPtr<labelList>::New(offsets);
        vtuSizing::renumberFaceOffsets(tresult.ref(), beginOffset, output_);
        return tresult;
    }
    else
    {
        auto tresult = refPtr<labelList>::New();
        auto& result = tresult.ref();

        // With VTKHDF, do not need dummy offsets
        if (contentType::HDF == output_)
        {
            return tresult;
        }

        const label numCells = cellTypes_.size();

        result = vtuSizing::dummyFaceOffsets(numCells, output_, beginOffset);

        return tresult;
    }
}


Foam::refPtr<Foam::labelList>
Foam::vtk::vtuCells::polyFaceIds(label beginOffset) const
{
    #ifdef Foam_vtk_vtuCells_demandDriven_POLY_FACEIDS
    if
    (
        const label numFaces =
        (
            (output_ == contentType::HDF) ? vtuSizing::nFacesPoly() : 0
        );
        (polyFaceIds_.size() < numFaces)
    )
    {
        auto tresult = refPtr<labelList>::New(numFaces);
        auto& result = tresult.ref();
        std::iota(result.begin(), result.end(), beginOffset);
        return tresult;
    }
    #endif

    if (beginOffset <= 0)
    {
        return refPtr<labelList>(polyFaceIds_);
    }

    auto tresult = refPtr<labelList>::New(polyFaceIds_);
    auto& result = tresult.ref();

    if (beginOffset > 0)
    {
        for (auto& val : result)
        {
            val += beginOffset;
        }
    }

    return tresult;
}


Foam::refPtr<Foam::labelList>
Foam::vtk::vtuCells::polyFaceOffsets
(
    label beginOffset,
    bool syncPar
) const
{
    // If the format uses begin/end offsets, need special care
    // when concatenating.

    const bool uses_beginEnd_style =
    (
        (output_ == contentType::HDF)
    );

    //const auto& offsets = polyFaceOffset_;

    const int communicator = UPstream::worldComm;

    if (syncPar && UPstream::parRun() && uses_beginEnd_style)
    {
        auto tresult = refPtr<labelList>::New();
        auto& result = tresult.ref();

        // Can assume that we enter this branch when the outside caller
        // has determined that polyhedrons are involved. So we can avoid
        // any additional reduction here.
        //
        // For VTKHDF format, need to have PolyhedronOffsets defined
        // everywhere as soon as *any* partition has polyhedrons.

        // Slightly wasteful with the extra allocations here,
        // but needed for the logic when concatenating the offsets.

        labelList dummyOffsets;

        if
        (
            (output_ == contentType::HDF)
         && (polyFaceOffset_.empty() && nCells() > 0)
        )
        {
            dummyOffsets.resize(nCells()+1, 0);
        }

        const auto& offsets =
        (
            dummyOffsets.empty() ? polyFaceOffset_ : dummyOffsets
        );

        result = parCoordinatedOffsets(offsets, communicator);

        // Determine a consistent beginOffset
        if (beginOffset < 0)
        {
            beginOffset = UPstream::find_offset<label>
            (
                (offsets.size() > 1 ? offsets.back() : 0),
                communicator
            );
        }

        if (beginOffset > 0)
        {
            for (auto& val : result)
            {
                val += beginOffset;
            }
        }

        return tresult;
    }
    else
    {
        const auto& offsets = polyFaceOffset_;

        if (beginOffset <= 0)
        {
            return refPtr<labelList>(offsets);
        }

        auto tresult = refPtr<labelList>::New(offsets);
        auto& result = tresult.ref();

        if (beginOffset > 0)
        {
            for (auto& val : result)
            {
                val += beginOffset;
            }
        }

        return tresult;
    }
}


void Foam::vtk::vtuCells::mergePoints(const polyMesh& mesh)
{
    // NB: caller is responsible for ensuring consistent merge states
    // on all ranks!

    if (UPstream::parRun() && (merged_ == vtk::vtuSizing::mergeType::NONE))
    {
        labelList pointToGlobal;
        labelList uniqueMeshPoints;

        const autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints(pointToGlobal, uniqueMeshPoints);

        if (returnReduceOr(this->nAddPoints() > 0))
        {
            // Needs additional effort with decomposed polyhedrals...
            const globalIndex& gi = globalPointsPtr();

            // Bookkeeping for the additional cell-centre points
            const globalIndex giAdded
            (
                this->nAddPoints(),
                UPstream::worldComm,
                UPstream::parRun()
            );

            const auto myProci = UPstream::myProcNo(UPstream::worldComm);

            // Re-index to account for the change in offsets due to
            // the additional cell-centre points
            for (auto& globalPointi : pointToGlobal)
            {
                if
                (
                    auto proci = gi.findProc(myProci, globalPointi);
                    proci >= 0
                )
                {
                    globalPointi += giAdded.localStart(proci);
                }
            }

            // Extend to include the additional points.
            // These are contiguous addressing starting after the number of
            // regular (unique, merged) mesh points.

            const label nOrigPoints = pointToGlobal.size();
            pointToGlobal.resize(nOrigPoints + this->nAddPoints());

            auto endSlice = pointToGlobal.slice(nOrigPoints);
            Foam::identity
            (
                endSlice,
                (giAdded.localStart(myProci) + gi.localEnd(myProci))
            );
        }

        vtuSizing::renumberVertLabels(vertLabels_, pointToGlobal, output_);
        vtuSizing::renumberFaceLabels(faceLabels_, pointToGlobal, output_);

        // The unique mesh points define the point map
        maps_.pointMap() = std::move(uniqueMeshPoints);
        this->setNumPoints(maps_.pointMap().size());

        merged_ = vtk::vtuSizing::mergeType::MERGED;
    }
}


// ************************************************************************* //
