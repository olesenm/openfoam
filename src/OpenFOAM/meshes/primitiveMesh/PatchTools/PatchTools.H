/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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

Class
    Foam::PatchTools

Description
    A collection of tools for searching, sorting PrimitivePatch information.

    The class could also be extended to include more than just static methods.

SourceFiles
    PatchTools.C
    PatchToolsCheck.C
    PatchToolsEdgeOwner.C
    PatchToolsGatherAndMerge.C
    PatchToolsMatch.C
    PatchToolsNormals.C
    PatchToolsSearch.C
    PatchToolsSortEdges.C
    PatchToolsSortPoints.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PatchTools_H
#define Foam_PatchTools_H

#include "autoPtr.H"
#include "globalIndex.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class bitSet;
class boundBox;

/*---------------------------------------------------------------------------*\
                         Class PatchTools Declaration
\*---------------------------------------------------------------------------*/

class PatchTools
{
    //- Check for orientation issues.
    //  Returns true if problems were found.
    //  If a normal flips across an edge, places it in the \tparam Inserter
    template<class FaceList, class PointField, class Inserter>
    static bool checkOrientation_impl
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const bool report,
        //! [out] edges that flipped
        Inserter inserter
    );

    //- Fill faceZone with \p currentZone for every face reachable
    //- from \p facei without crossing edge marked in \b borderEdge.
    //  \note \p faceZone has to be sized nFaces before calling.
    template<class BoolListType, class FaceList, class PointField>
    static void markZone_impl
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const BoolListType& borderEdge,
        const label facei,
        const label currentZone,
        labelUList& faceZone,
        DynamicList<label>& work0,
        DynamicList<label>& work1
    );


public:

    //- Check for orientation issues.
    //  Returns true if problems were found.
    //  If a normal flips across an edge, places it in the HashSet
    template<class FaceList, class PointField>
    static bool checkOrientation
    (
        const PrimitivePatch<FaceList, PointField>& p,
        bool report = false,
        labelHashSet* marked = nullptr
    );

    //- Check for orientation issues.
    //  Returns true if problems were found.
    //  If a normal flips across an edge, places it in the HashSet
    template<class FaceList, class PointField>
    static bool checkOrientation
    (
        const PrimitivePatch<FaceList, PointField>& p,
        labelHashSet& marked,
        bool report = false
    );

    //- Check for orientation issues.
    //  Returns true if problems were found.
    //  If a normal flips across an edge, places it in the HashSet
    template<class FaceList, class PointField>
    static bool checkOrientation
    (
        const PrimitivePatch<FaceList, PointField>& p,
        bitSet& marked,
        const bool report = false
    );

    //- Check that faces are in compact face walk order, no gaps etc.
    //
    //  Returns true if the face points are in strictly ascending order
    template<class FaceList>
    static bool checkFaceWalkOrder
    (
        const FaceList& faces,
        //! The total of support points for the faces.
        //! Use -1 to suppress this final check
        const label nPoints = -1
    );

    //- Check if the patch support points have a 1-to-1 correspondence
    //- to the local points. This implies face walk order, no gaps etc.
    //
    //  Returns true if the face points are in strictly ascending order
    template<class FaceList, class PointField>
    static bool checkCompactPoints
    (
        const PrimitivePatch<FaceList, PointField>& p
    );

    //- Fill faceZone with \p currentZone for every face reachable
    //- from \p facei without crossing edge marked in \b borderEdge.
    //  \note \p faceZone has to be sized nFaces before calling.
    template<class BoolListType, class FaceList, class PointField>
    static void markZone
    (
        const PrimitivePatch<FaceList, PointField>&,
        const BoolListType& borderEdge,
        const label facei,
        const label currentZone,
        labelUList& faceZone
    );

    //- Size and fills faceZone with zone of face.
    //  Zone is area reachable by edge crossing without crossing borderEdge.
    //  Returns number of zones.
    template<class BoolListType, class FaceList, class PointField>
    static label markZones
    (
        const PrimitivePatch<FaceList, PointField>&,
        const BoolListType& borderEdge,
        labelList& faceZone
    );

    //- Determine the mapping for a sub-patch.
    //- Only include faces for which bool-list entry is true.
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  faceSelection faces to include
    //  \param[out] faceMap  mapping new to old faces
    //
    //  The \p faceSelection face selection type is anything that
    //  can be tested like a bitSet/boolList/labelHashSet
    template<class BoolListType, class FaceList, class PointField>
    static void subsetMapFaces
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const BoolListType& faceSelection,
        labelList& faceMap
    );

    //- Determine the mapping for a sub-patch.
    //- Only include faces for which bool-list entry is true.
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  faceSelection faces to include
    //  \param[out] pointMap mapping new to old localPoints
    //  \param[out] faceMap  mapping new to old faces
    //
    //  The \p includeFaces face selection type is anything that
    //  can be tested like a bitSet/boolList/labelHashSet
    template<class BoolListType, class FaceList, class PointField>
    static void subsetMap
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const BoolListType& faceSelection,
        labelList& pointMap,
        labelList& faceMap
    );

    //- Determine the mapping for a sub-patch based on its colouring.
    //- Only include faces with a matching colour.
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  colour select within \p colouring
    //  \param[in]  colouring the face colours (eg from a decomposition method)
    //  \param[out] faceMap  mapping new to old faces
    //
    //  \note if the colouring is shorter than the number of faces,
    //     the out-of-bounds faces are simply ignored (not selected).
    //     This is consistent with the boolean selector versions.
    template<class FaceList, class PointField>
    static void subsetMapFaces
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const label colour,
        const labelUList& colouring,
        labelList& faceMap
    );

    //- Determine the mapping for a sub-patch based on its colouring.
    //- Only include faces with a matching colour.
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  colour select within \p colouring
    //  \param[in]  colouring the face colours (eg from a decomposition method)
    //  \param[out] pointMap mapping new to old localPoints
    //  \param[out] faceMap  mapping new to old faces
    //
    //  \note if the colouring is shorter than the number of faces,
    //     the out-of-bounds faces are simply ignored (not selected).
    //     This is consistent with the boolean selector versions.
    template<class FaceList, class PointField>
    static void subsetMap
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const label colour,
        const labelUList& colouring,
        labelList& pointMap,
        labelList& faceMap
    );

    //- Determine the mapping for a sub-patch.
    //- Only include faces for which the predicate tests true (on face index)
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  pred acceptance predicate testing on the face \em index
    //  \param[out] faceMap  mapping new to old faces
    template<class FaceList, class PointField, class Predicate>
    static void subsetMapFaces_if
    (
        const PrimitivePatch<FaceList, PointField>& p,
        Predicate pred,
        labelList& faceMap
    );

    //- Determine the mapping for a sub-patch.
    //- Only include faces for which the predicate tests true (on face index)
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  pred acceptance predicate testing on the face \em index
    //  \param[out] pointMap mapping new to old localPoints
    //  \param[out] faceMap  mapping new to old faces
    template<class FaceList, class PointField, class Predicate>
    static void subsetMap_if
    (
        const PrimitivePatch<FaceList, PointField>& p,
        Predicate pred,
        labelList& pointMap,
        labelList& faceMap
    );

    //- Determine the mapping for a sampled sub-patch.
    //- This differs from subsetMap in that the \p faceMap is used
    //- to determine the locations to select.
    //
    //  On input the \p faceMap specifies the old faces to be selected.
    //  This map is sorted, de-duplicated and bounds checked.
    //  After this sanitizing step, it then becomes the output face map.
    //
    //  \param[in]  p patch to be searched on
    //  \param[in,out] faceMap  mapping new to old faces
    template<class FaceList, class PointField>
    static void sampleMapFaces
    (
        const PrimitivePatch<FaceList, PointField>& p,
        labelList& faceMap
    );

    //- Determine the mapping for a sampled sub-patch.
    //- This differs from subsetMap in that the \p faceMap is used
    //- to determine the locations to select.
    //
    //  On input the \p faceMap specifies the old faces to be selected.
    //  This map is sorted, de-duplicated and bounds checked.
    //  After this sanitizing step, it then becomes the output face map.
    //
    //  \param[in]  p patch to be searched on
    //  \param[out] pointMap mapping new to old localPoints
    //  \param[in,out] faceMap  mapping new to old faces
    template<class FaceList, class PointField>
    static void sampleMap
    (
        const PrimitivePatch<FaceList, PointField>& p,
        labelList& pointMap,
        labelList& faceMap
    );

    //- Determine the local point mapping for a sub-patch.
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  faceMap the mapping new to old faces of the faces to use
    //  \param[out] pointMap mapping new to old localPoints
    template<class FaceList, class PointField>
    static void subsetMapPoints
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const labelUList& faceMap,
        labelList& pointMap
    );

    //- Determine the (non-local) point mapping for a sub-patch.
    //- This can be used to extract faces/points without
    //- an intermediate meshPoints() or localFaces().
    //
    //  \param[in]  p patch to be searched on
    //  \param[in]  faceMap the mapping new to old faces of the faces to use
    //  \param[out] pointMap mapping new to old (non-local) patch points
    //  \param[out] reversePointMap mapping old (non-local) patch points to new,
    //      which is normally used as a work array.
    template<class FaceList, class PointField>
    static void subsetMapPoints
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const labelUList& faceMap,
        labelList& pointMap,
        DynamicList<label>& reversePointMap
    );

    //- Copy a subset of faces (regular or local) and renumber their
    //- point labels using the \p oldToNewPoints mapping.
    //  The subset is given by the \p faceMap.
    //
    //  \param[in] inputFaces the regular/local faces to copy
    //  \param[out] newFaces the subsetted, renumbered faces as
    //      List/DynamicList etc of the face type
    //  \param[in] oldToNewPoints mapping of the old point labels
    //      (consistent with the \p inputFaces numbering) to new points.
    //      Must be valid for all points used by the selected faces.
    //  \param[in] faceMap  mapping new to old faces
    template<class OutFaceListType, class FaceList>
    static void copyMapFaces
    (
        const FaceList& inputFaces,
        OutFaceListType& newFaces,
        const labelUList& oldToNewPoints,
        const labelUList& faceMap
    );

    //- Copy a subset of faces from the patch.
    //- The subset if given by the \p faceMap and the new faces are
    //- renumbered according to the \p pointMap compact numbering
    //
    //  \param[in] p patch to copy
    //  \param[out] newFaces the subsetted, renumbered faces as
    //      List/DynamicList etc of the PrimitivePatch::face_type.
    //  \param[in] pointMap mapping new to old points
    //  \param[in] faceMap  mapping new to old faces
    template<class OutFaceListType, class FaceList, class PointField>
    static void copyFaces
    (
        const PrimitivePatch<FaceList, PointField>& p,
        OutFaceListType& newFaces,
        const labelUList& pointMap,
        const labelUList& faceMap
    );

    //- Copy a subset of localFaces from the patch.
    //- The subset if given by the \p faceMap and the new faces are
    //- renumbered according to the \p pointMap compact numbering
    //
    //  \param[in] p patch to copy
    //  \param[out] newFaces the subsetted, renumbered faces as
    //      List/DynamicList etc of the PrimitivePatch::face_type.
    //  \param[in] pointMap mapping new to old localPoints
    //  \param[in] faceMap  mapping new to old faces
    template<class OutFaceListType, class FaceList, class PointField>
    static void copyLocalFaces
    (
        const PrimitivePatch<FaceList, PointField>& p,
        OutFaceListType& newFaces,
        const labelUList& pointMap,
        const labelUList& faceMap
    );

    //- Copy a subset of localPoints from the patch, using the
    //- \p pointMap to select the points.
    //- Avoids triggering demand-driven meshPoints()/localFaces()
    //- calculation when they do not already exist.
    //
    //  \param[in] p patch to copy
    //  \param[out] newPoints the subsetted points as
    //      Field/List/DynamicList etc of the PrimitivePatch::point_type.
    //  \param[in] pointMap mapping new to old localPoints
    template<class OutPointListType, class FaceList, class PointField>
    static void copyLocalPoints
    (
        const PrimitivePatch<FaceList, PointField>& p,
        OutPointListType& newPoints,
        const labelUList& pointMap
    );

    //- Copy a subset of localFaces and localPoints from the patch.
    //- The subset is given by the \p faceMap (unique and sorted)
    //- and \p pointMap.
    //- Combines copyLocalFaces and copyLocalPoints.
    //
    //  \param[in] p patch to copy
    //  \param[out] newFaces the subsetted, renumbered faces
    //  \param[out] newPoints the subsetted points
    //  \param[in] pointMap mapping new to old localPoints
    //  \param[in] faceMap  mapping new to old faces
    template
    <
        class OutFaceListType,
        class OutPointListType,
        class FaceList,
        class PointField
    >
    static void copyLocalSubset
    (
        const PrimitivePatch<FaceList, PointField>& p,
        OutFaceListType& newFaces,
        OutPointListType& newPoints,
        const labelUList& pointMap,
        const labelUList& faceMap
    );

    //- Copy a subset of (local) faces from the patch.
    //- The subset is given by the \p faceMap (unique and sorted).
    //
    //  \param[in] p patch to copy
    //  \param[out] newFaces the subsetted, renumbered faces
    //  \param[out] newPoints the subsetted points
    //  \param[in] faceMap  mapping new to old faces
    template
    <
        class OutFaceListType,
        class OutPointListType,
        class FaceList,
        class PointField
    >
    static void copySubset
    (
        const PrimitivePatch<FaceList, PointField>& p,
        OutFaceListType& newFaces,
        OutPointListType& newPoints,
        const labelUList& faceMap
    )
    {
        labelList pointMap;
        PatchTools::subsetMapPoints(p, faceMap, pointMap);
        PatchTools::copyLocalSubset(p, newFaces, newPoints, pointMap, faceMap);
    }


    //- Calculate the bounding box of the patch without meshPoints etc.
    template<class FaceList, class PointField>
    static void calcBounds
    (
        const PrimitivePatch<FaceList, PointField>& p,
        //! [out] the bounding box
        boundBox& bb,
        //! [out] the number of unique points used by the patch
        label& nPoints
    );

    //- Return edge-face addressing sorted by angle around the edge.
    //  Orientation is anticlockwise looking from edge.vec(localPoints())
    template<class FaceList, class PointField>
    static labelListList sortedEdgeFaces
    (
        const PrimitivePatch<FaceList, PointField>&
    );

    //- Return point-edge addressing sorted by order around the point.
    template<class FaceList, class PointField>
    static labelListList sortedPointEdges
    (
        const PrimitivePatch<FaceList, PointField>&
    );

    //- If 2 face neighbours: label of face where ordering of edge
    //  is consistent with righthand walk.
    //  If 1 neighbour: label of only face.
    //  If >2 neighbours: undetermined.
    template<class FaceList, class PointField>
    static labelList edgeOwner
    (
        const PrimitivePatch<FaceList, PointField>&
    );


    //- Find corresponding points on patches sharing the same points
    //  p1PointLabels : points on p1 that were matched
    //  p2PointLabels : corresponding points on p2
    template
    <
        class FaceList1, class PointField1,
        class FaceList2, class PointField2
    >
    static void matchPoints
    (
        const PrimitivePatch<FaceList1, PointField1>& p1,
        const PrimitivePatch<FaceList2, PointField2>& p2,

        labelList& p1PointLabels,
        labelList& p2PointLabels
    );

    //- Find corresponding edges on patches sharing the same points
    //  p1EdgeLabels    : edges on p1 that were matched
    //  p2EdgeLabels    : corresponding edges on p2
    //  sameOrientation : same orientation?
    template
    <
        class FaceList1, class PointField1,
        class FaceList2, class PointField2
    >
    static void matchEdges
    (
        const PrimitivePatch<FaceList1, PointField1>& p1,
        const PrimitivePatch<FaceList2, PointField2>& p2,

        labelList& p1EdgeLabels,
        labelList& p2EdgeLabels,
        bitSet& sameOrientation
    );


    //- Return parallel consistent point normals for patches using mesh points.
    template<class FaceList, class PointField>
    static tmp<pointField> pointNormals
    (
        const polyMesh&,
        const PrimitivePatch<FaceList, PointField>&,
        const bitSet& flipMap = bitSet::null()
    );

    //- Return parallel consistent point normals for patches using mesh points.
    //  Get all geometry from provided patch local points. (assumed to be
    //  consistent across coupled patches)
    template<class FaceList, class PointField>
    static tmp<pointField> pointNormals
    (
        const polyMesh&,
        const PrimitivePatch<FaceList, PointField>&,
        const pointField& localPoints,
        const bitSet& flipMap = bitSet::null()
    );

    //- Return parallel consistent edge normals for patches using mesh points.
    //  Supply with patch matching info from matchEdges.
    template<class FaceList, class PointField>
    static tmp<pointField> edgeNormals
    (
        const polyMesh&,
        const PrimitivePatch<FaceList, PointField>&,
        const labelList& patchEdges,
        const labelList& coupledEdges,
        const bitSet& flipMap = bitSet::null()
    );

    //- Return parallel consistent edge normals for patches using mesh points.
    //  Get all geometry from provided patch local points. (assumed to be
    //  consistent across coupled patches)
    //  Supply with patch matching info from matchEdges.
    template<class FaceList, class PointField>
    static tmp<pointField> edgeNormals
    (
        const polyMesh&,
        const PrimitivePatch<FaceList, PointField>&,
        const pointField& localPoints,
        const labelList& patchEdges,
        const labelList& coupledEdges,
        const bitSet& flipMap = bitSet::null()
    );


    //- Gather points and faces onto master and merge into single patch.
    //  Note: Normally uses faces/points (not localFaces/localPoints)
    //
    //  \param[in] mergeDist Geometric merge tolerance for Foam::mergePoints
    //  \param[in] pp The patch to merge
    //  \param[out] mergedPoints merged points (master only, empty elsewhere)
    //  \param[out] mergedFaces merged faces (master only, empty elsewhere)
    //  \param[out] pointAddr Points globalIndex gather addressing
    //      (master only, empty elsewhere)
    //  \param[out] faceAddr Faces globalIndex gather addressing
    //      (master only, empty elsewhere)
    //  \param[out] pointMergeMap An old-to-new mapping from original
    //      point index to the index into merged points.
    //  \param[in]  useLocal  gather/merge patch localFaces/localPoints
    //      instead of faces/points
    //
    //  \note
    //  - OpenFOAM-v2112 and earlier: geometric merge on all patch points.
    //  - OpenFOAM-v2206 and later: geometric merge on patch boundary points.
    template<class FaceList, class PointField>
    static void gatherAndMerge
    (
        const scalar mergeDist,
        const PrimitivePatch<FaceList, PointField>& pp,
        Field
        <
            typename PrimitivePatch<FaceList, PointField>::point_type
        >& mergedPoints,
        List
        <
            typename PrimitivePatch<FaceList, PointField>::face_type
        >& mergedFaces,
        globalIndex& pointAddr,
        globalIndex& faceAddr,
        labelList& pointMergeMap = const_cast<labelList&>(labelList::null()),
        const bool useLocal = false
    );

    //- Gather points and faces onto master and merge into single patch.
    //  Note: Normally uses faces/points (not localFaces/localPoints)
    //
    //  \param[in] mergeDist Geometric merge tolerance for Foam::mergePoints
    //  \param[in] pp The patch to merge
    //  \param[out] mergedPoints merged points (master only, empty elsewhere)
    //  \param[out] mergedFaces merged faces (master only, empty elsewhere)
    //  \param[out] pointMergeMap An old-to-new mapping from original
    //      point index to the index into merged points.
    //  \param[in] useLocal gather/merge patch localFaces/localPoints
    //      instead of faces/points
    //
    //  \note
    //  - OpenFOAM-v2112 and earlier: geometric merge on all patch points.
    //  - OpenFOAM-v2206 and later: geometric merge on patch boundary points.
    template<class FaceList, class PointField>
    static void gatherAndMerge
    (
        const scalar mergeDist,
        const PrimitivePatch<FaceList, PointField>& pp,
        Field
        <
            typename PrimitivePatch<FaceList, PointField>::point_type
        >& mergedPoints,
        List
        <
            typename PrimitivePatch<FaceList, PointField>::face_type
        >& mergedFaces,
        labelList& pointMergeMap = const_cast<labelList&>(labelList::null()),
        const bool useLocal = false
    );

    //- Gather (mesh!) points and faces onto master and merge collocated
    //  points into a single patch. Uses coupled point mesh
    //  structure so does not need tolerances.
    //  On master and sub-ranks returns:
    //  - pointToGlobal : for every local point index the global point index
    //  - uniqueMeshPointLabels : my local mesh points
    //  - globalPoints : global numbering for the global points
    //  - globalFaces : global numbering for the faces
    //  .
    //  On master only returns:
    //  - mergedFaces : the merged faces
    //  - mergedPoints : the merged points
    template<class FaceList>
    static void gatherAndMerge
    (
        const polyMesh& mesh,
        const FaceList& faces,
        const labelList& meshPoints,
        const Map<label>& meshPointMap,

        labelList& pointToGlobal,
        labelList& uniqueMeshPointLabels,
        autoPtr<globalIndex>& globalPoints,
        autoPtr<globalIndex>& globalFaces,
        List<typename FaceList::value_type>& mergedFaces,
        pointField& mergedPoints
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchTools.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
