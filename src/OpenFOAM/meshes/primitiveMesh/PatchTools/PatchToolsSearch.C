/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

Description
    Searching and marking zones of the patch.

\*---------------------------------------------------------------------------*/

#include "PatchTools.H"
#include "bitSet.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class BoolListType, class FaceList, class PointField>
void Foam::PatchTools::markZone_impl
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& borderEdge,
    const label facei,
    const label currentZone,
    labelUList& faceZone,

    DynamicList<label>& work0,
    DynamicList<label>& work1
)
{
    const labelListList& faceEdges = p.faceEdges();
    const labelListList& edgeFaces = p.edgeFaces();

    // Faces to be visited vs next faces to be visited (wavefront/queue).
    //
    // The two buffers are ping-ponged via swap() to reuse their allocated
    // capacity each pass.

    auto& pendingFaces = work0;
    auto& nextFaces = work1;

    pendingFaces.clear();
    nextFaces.clear();

    // Seed with a single face
    pendingFaces.push_back(facei);

    while (!pendingFaces.empty())
    {
        nextFaces.clear();

        // Pick up neighbours of pendingFaces
        for (const label facei : pendingFaces)
        {
            for (const label edgei : faceEdges[facei])
            {
                if (!borderEdge.test(edgei))
                {
                    for (const label nbrFacei : edgeFaces[edgei])
                    {
                        if (faceZone[nbrFacei] == -1)
                        {
                            faceZone[nbrFacei] = currentZone;
                            nextFaces.push_back(nbrFacei);
                        }
                        else if (faceZone[nbrFacei] != currentZone)
                        {
                            FatalErrorInFunction
                                << "Zones " << faceZone[nbrFacei]
                                << " at face " << nbrFacei
                                << " connects to zone " << currentZone
                                << " at face " << facei
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        // Swap the queue buffers (O(1) - no reallocation)
        pendingFaces.swap(nextFaces);
    }
}


template<class BoolListType, class FaceList, class PointField>
void Foam::PatchTools::markZone
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& borderEdge,
    const label facei,
    const label currentZone,
    labelUList& faceZone
)
{
    DynamicList<label> work0, work1;

    PatchTools::markZone_impl
    (
        p,
        borderEdge,
        facei,
        currentZone,
        faceZone,
        work0,
        work1
    );
}


template<class BoolListType, class FaceList, class PointField>
Foam::label Foam::PatchTools::markZones
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& borderEdge,
    labelList& faceZone
)
{
    // NB: even if borderEdge is empty or only false, the patch may still
    // be split into multiple disconnected zones, so always need to walk.

    DynamicList<label> work0, work1;  //<- work arrays

    faceZone.resize_fill(p.size(), -1);
    label zoneIndex = 0;

    // Find unvisited faces, mark the zone, find the next unvisited face.
    // Needs to check for (seedFacei >= 0) since find() returns -1 when
    // nothing more is found.
    for
    (
        label seedFacei = 0;
        (seedFacei >= 0 && seedFacei < faceZone.size());
        ++zoneIndex
    )
    {
        faceZone[seedFacei] = zoneIndex;

        PatchTools::markZone_impl
        (
            p,
            borderEdge,
            seedFacei,
            zoneIndex,
            faceZone,
            work0, work1
        );

        // Continue from the next unvisited face
        ++seedFacei;
        seedFacei = faceZone.find(-1, seedFacei);
    }

    return zoneIndex;
}


// Generate faceMap based on boolean-like selection
template<class BoolListType, class FaceList, class PointField>
void Foam::PatchTools::subsetMapFaces
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& faceSelection,
    labelList& faceMap
)
{
    faceMap.resize_nocopy(p.size());

    label newFacei = 0;

    forAll(faceMap, oldFacei)
    {
        if (faceSelection.test(oldFacei))
        {
            // Compact storage for new faces
            faceMap[newFacei++] = oldFacei;
        }
    }

    faceMap.resize(newFacei);
}


// Generate faceMap based on predicate selector on face index
template<class FaceList, class PointField, class Predicate>
void Foam::PatchTools::subsetMapFaces_if
(
    const PrimitivePatch<FaceList, PointField>& p,
    Predicate pred,
    labelList& faceMap
)
{
    faceMap.resize_nocopy(p.size());

    label newFacei = 0;

    forAll(faceMap, oldFacei)
    {
        if (pred(oldFacei))
        {
            // Compact storage for new faces
            faceMap[newFacei++] = oldFacei;
        }
    }

    faceMap.resize(newFacei);
}


// Generate faceMap based on colouring
template<class FaceList, class PointField>
void Foam::PatchTools::subsetMapFaces
(
    const PrimitivePatch<FaceList, PointField>& p,
    const label colour,
    const labelUList& colouring,
    labelList& faceMap
)
{
    // If the colouring is shorter than the number of faces,
    // the out-of-bounds faces are simply ignored (not selected).
    // This is consistent with the boolean selector versions.

    const label loopLen = std::min(p.size(), colouring.size());

    faceMap.resize_nocopy(loopLen);

    label newFacei = 0;

    for
    (
        label oldFacei = colouring.find(colour);
        (oldFacei >= 0 && oldFacei < loopLen);
        /*nil*/
    )
    {
        // Compact storage for new faces
        faceMap[newFacei++] = oldFacei;

        // Find the next with this colour
        ++oldFacei;
        oldFacei = colouring.find(colour, oldFacei);
    }

    faceMap.resize(newFacei);
}


// Generate/re-generate the faceMap into sanitized canonical form
template<class FaceList, class PointField>
void Foam::PatchTools::sampleMapFaces
(
    const PrimitivePatch<FaceList, PointField>& p,
    labelList& faceMap
)
{
    const label maxLen = p.size();

    Foam::sort(faceMap);

    label newFacei = 0;
    label prevFacei = -1;

    for (label oldFacei : faceMap)
    {
        if (oldFacei >= maxLen)
        {
            break;
        }

        // bounds checking and de-duplication
        if ((oldFacei >= 0) && (oldFacei != prevFacei))
        {
            // Compact storage for new faces
            faceMap[newFacei++] = oldFacei;

            prevFacei = oldFacei;
        }
    }

    faceMap.resize(newFacei);
}


template<class FaceList, class PointField>
void Foam::PatchTools::subsetMapPoints
(
    const PrimitivePatch<FaceList, PointField>& p,
    const labelUList& faceMap,
    labelList& pointMap
)
{
    pointMap.clear();

    if (faceMap.empty())
    {
        return;
    }

    bitSet pointUsed;

    if (p.hasLocalFaces())
    {
        const auto& localFaces = p.localFaces();

        pointUsed.resize_fill(p.meshPoints().size(), false);

        for (label oldFacei : faceMap)
        {
            // Local points used by face
            pointUsed.set(localFaces[oldFacei]);
        }
    }
    else if (PatchTools::checkCompactPoints(p))
    {
        // Compact points: local == global numbering, no mapping needed.

        pointUsed.resize_fill(p.points().size(), false);

        for (label oldFacei : faceMap)
        {
            pointUsed.set(p[oldFacei]);
        }
    }
    else
    {
        // No localFaces - two stages:
        // 1. fake meshPoints() as a dense map (non-local -> local)
        //    using the face walk order of the entire patch
        // 2. subset using the local point numbering

        labelList oldToLocal(p.points().size(), -1);

        label nLocalPoints = 0;

        for (const auto& f : p)
        {
            for (label oldPointi : f)
            {
                if (oldToLocal[oldPointi] == -1)
                {
                    oldToLocal[oldPointi] = nLocalPoints;
                    ++nLocalPoints;
                }
            }
        }

        pointUsed.resize_fill(nLocalPoints, false);

        for (label oldFacei : faceMap)
        {
            for (label oldPointi : p[oldFacei])
            {
                pointUsed.set(oldToLocal[oldPointi]);
            }
        }
    }

    // The newToOld mappings
    pointMap = pointUsed.sortedToc();
}


template<class FaceList, class PointField>
void Foam::PatchTools::subsetMapPoints
(
    const PrimitivePatch<FaceList, PointField>& p,
    const labelUList& faceMap,
    labelList& pointMap,
    DynamicList<label>& reversePointMap
)
{
    if (faceMap.empty())
    {
        pointMap.clear();
        reversePointMap.clear();
        return;
    }

    pointMap.resize_nocopy(p.points().size());
    reversePointMap.resize_fill(p.points().size(), -1);

    label pointi = 0;
    label maxPointi = -1;

    for (label oldFacei : faceMap)
    {
        // Face label renumbering (walk order)
        for (label oldPointi : p[oldFacei])
        {
            if (reversePointMap[oldPointi] == -1)
            {
                maxPointi = std::max(maxPointi, oldPointi);
                reversePointMap[oldPointi] = pointi;
                pointMap[pointi] = oldPointi;
                ++pointi;
            }
        }
    }

    reversePointMap.resize(maxPointi+1);
    pointMap.resize(pointi);
}


template<class BoolListType, class FaceList, class PointField>
void Foam::PatchTools::subsetMap
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& faceSelection,
    labelList& pointMap,
    labelList& faceMap
)
{
    pointMap.clear();

    // Get the corresponding face indices selected
    PatchTools::subsetMapFaces(p, faceSelection, faceMap);
    PatchTools::subsetMapPoints(p, faceMap, pointMap);
}


template<class FaceList, class PointField, class Predicate>
void Foam::PatchTools::subsetMap_if
(
    const PrimitivePatch<FaceList, PointField>& p,
    Predicate pred,
    labelList& pointMap,
    labelList& faceMap
)
{
    pointMap.clear();

    // Get the corresponding face indices selected
    PatchTools::subsetMapFaces_if(p, pred, faceMap);
    PatchTools::subsetMapPoints(p, faceMap, pointMap);
}


template<class FaceList, class PointField>
void Foam::PatchTools::subsetMap
(
    const PrimitivePatch<FaceList, PointField>& p,
    const label colour,
    const labelUList& colouring,
    labelList& pointMap,
    labelList& faceMap
)
{
    pointMap.clear();

    // Get the corresponding face indices selected
    PatchTools::subsetMapFaces(p, colour, colouring, faceMap);
    PatchTools::subsetMapPoints(p, faceMap, pointMap);
}


template<class FaceList, class PointField>
void Foam::PatchTools::sampleMap
(
    const PrimitivePatch<FaceList, PointField>& p,
    labelList& pointMap,
    labelList& faceMap
)
{
    // Get the corresponding face indices selected
    PatchTools::sampleMapFaces(p, faceMap);
    PatchTools::subsetMapPoints(p, faceMap, pointMap);
}


template<class OutFaceListType, class FaceList>
void Foam::PatchTools::copyMapFaces
(
    const FaceList& inputFaces,
    OutFaceListType& newFaces,
    const labelUList& oldToNewPoints,
    const labelUList& faceMap
)
{
    // Subset of faces.
    // Deep copy the local faces and apply face node renumbering
    newFaces.resize_nocopy(faceMap.size());

    forAll(faceMap, facei)
    {
        auto& f = newFaces[facei];

        // Deep copy the old face
        f = inputFaces[faceMap[facei]];

        // Renumber face node labels
        for (label& pointi : f)
        {
            pointi = oldToNewPoints[pointi];
        }
    }
}


template<class OutFaceListType, class FaceList, class PointField>
void Foam::PatchTools::copyFaces
(
    const PrimitivePatch<FaceList, PointField>& p,
    OutFaceListType& newFaces,
    const labelUList& pointMap,
    const labelUList& faceMap
)
{
    static_assert
    (
        // could relax this to std::is_base_of_v<output, input>
        std::is_same_v
        <
            typename OutFaceListType::value_type,
            typename PrimitivePatch<FaceList, PointField>::face_type
        >,
        "output list value_type must match the patch face_type"
    );

    // Face node label renumbering
    // Inverse point mapping - same as ListOps invert() without checks

    labelList oldToNewPoints(p.points().size(), -1);
    forAll(pointMap, pointi)
    {
        oldToNewPoints[pointMap[pointi]] = pointi;
    }

    PatchTools::copyMapFaces
    (
        p,  // real faces
        newFaces,
        oldToNewPoints,
        faceMap
    );
}


template<class OutFaceListType, class FaceList, class PointField>
void Foam::PatchTools::copyLocalFaces
(
    const PrimitivePatch<FaceList, PointField>& p,
    OutFaceListType& newFaces,
    const labelUList& pointMap,
    const labelUList& faceMap
)
{
    static_assert
    (
        // could relax this to std::is_base_of_v<output, input>
        std::is_same_v
        <
            typename OutFaceListType::value_type,
            typename PrimitivePatch<FaceList, PointField>::face_type
        >,
        "output list value_type must match the patch face_type"
    );

    // Face node label renumbering
    // Inverse point mapping - same as ListOps invert() without checks

    if (p.hasLocalFaces())
    {
        // Local faces exist: pointMap is in local point numbering
        labelList oldToNewPoints(p.meshPoints().size(), -1);
        forAll(pointMap, pointi)
        {
            oldToNewPoints[pointMap[pointi]] = pointi;
        }

        PatchTools::copyMapFaces
        (
            p.localFaces(),
            newFaces,
            oldToNewPoints,
            faceMap
        );
    }
    else if (PatchTools::checkCompactPoints(p))
    {
        // Compact points: local == global numbering
        labelList oldToNewPoints(p.points().size(), -1);
        forAll(pointMap, pointi)
        {
            oldToNewPoints[pointMap[pointi]] = pointi;
        }

        // copy the real faces
        PatchTools::copyMapFaces(p, newFaces, oldToNewPoints, faceMap);
    }
    else
    {
        // No localFaces and non-compact points.
        labelList oldToNewPoints(p.points().size(), -1);

        // pointMap refers to local (face walk order) numbering, but we will
        // use the regular (non-local) faces which have global points.
        // Thus compose global -> local -> new.

        // Mark the unused local mapping with -2 so that unused entries will
        // track as 'seen, but not unused' for the oldToNewPoints mapping

        labelList oldToNewLocal(p.points().size(), -2);  // -2, not -1 here.
        forAll(pointMap, pointi)
        {
            oldToNewLocal[pointMap[pointi]] = pointi;
        }

        label nLocalPoints = 0;

        for (const auto& f : p)
        {
            for (label oldPointi : f)
            {
                if (oldToNewPoints[oldPointi] == -1)
                {
                    // Local number is nLocalPoints (always < points size).
                    // Unselected points carry -2 (seen, but not used).
                    oldToNewPoints[oldPointi] = oldToNewLocal[nLocalPoints];
                    ++nLocalPoints;
                }
            }
        }

        // copy the real faces
        PatchTools::copyMapFaces(p, newFaces, oldToNewPoints, faceMap);
    }
}


template<class OutPointListType, class FaceList, class PointField>
void Foam::PatchTools::copyLocalPoints
(
    const PrimitivePatch<FaceList, PointField>& p,
    OutPointListType& newPoints,
    const labelUList& pointMap
)
{
    static_assert
    (
        std::is_same_v
        <
            typename OutPointListType::value_type,
            typename PrimitivePatch<FaceList, PointField>::point_type
        >,
        "output list value_type must match the patch point_type"
    );

    const auto& pts = p.points();

    if (pointMap.empty())
    {
        newPoints.clear();
        return;
    }

    newPoints.resize_nocopy(pointMap.size());

    if (p.hasMeshPoints())
    {
        // Already have local addressing, so no extra cost
        p.localPoints(pointMap, newPoints);
    }
    else if (PatchTools::checkCompactPoints(p))
    {
        // Compact points: local == global numbering
        forAll(pointMap, pointi)
        {
            newPoints[pointi] = pts[pointMap[pointi]];
        }
    }
    else
    {
        // No local addressing and non-compact points.
        // Walk the faces to establish the local (face walk order)
        // numbering without constructing meshPoints().

        // local -> new (-1: not selected)
        labelList localToNew(pts.size(), -1);
        forAll(pointMap, pointi)
        {
            localToNew[pointMap[pointi]] = pointi;
        }

        bitSet pointUsed(pts.size());
        label nLocalPoints = 0;

        for (const auto& f : p)
        {
            for (label oldPointi : f)
            {
                if (pointUsed.set(oldPointi))
                {
                    label newPointi = localToNew[nLocalPoints];
                    ++nLocalPoints;

                    if (newPointi >= 0)
                    {
                        newPoints[newPointi] = pts[oldPointi];
                    }
                }
            }
        }
    }
}


template
<
    class OutFaceListType,
    class OutPointListType,
    class FaceList,
    class PointField
>
void Foam::PatchTools::copyLocalSubset
(
    const PrimitivePatch<FaceList, PointField>& p,
    OutFaceListType& newFaces,
    OutPointListType& newPoints,
    const labelUList& pointMap,
    const labelUList& faceMap
)
{
    if
    (
        pointMap.empty()
     || p.hasMeshPoints() || PatchTools::checkCompactPoints(p)
    )
    {
        // Cheap cases (cached addressing, or local == global numbering)
        PatchTools::copyLocalPoints(p, newPoints, pointMap);
        PatchTools::copyLocalFaces(p, newFaces, pointMap, faceMap);
        return;
    }

    // No local addressing and non-compact points.
    // Single walk of the faces to compose global -> local -> new,
    // which is used for both the points and the faces.

    const auto& pts = p.points();

    // local -> new (-2: not selected, later as 'seen, but not used')
    labelList localToNew(pts.size(), -2);
    forAll(pointMap, pointi)
    {
        localToNew[pointMap[pointi]] = pointi;
    }

    // global -> new (-1: not yet visited)
    labelList oldToNewPoints(pts.size(), -1);

    newPoints.resize_nocopy(pointMap.size());

    label nLocalPoints = 0;

    for (const auto& f : p)
    {
        for (label oldPointi : f)
        {
            if (oldToNewPoints[oldPointi] == -1)
            {
                // First visit defines the local point number
                label newPointi = localToNew[nLocalPoints];
                ++nLocalPoints;

                oldToNewPoints[oldPointi] = newPointi;

                if (newPointi >= 0)
                {
                    newPoints[newPointi] = pts[oldPointi];
                }
            }
        }
    }

    PatchTools::copyMapFaces(p, newFaces, oldToNewPoints, faceMap);
}


template<class FaceList, class PointField>
void Foam::PatchTools::calcBounds
(
    const PrimitivePatch<FaceList, PointField>& p,
    boundBox& bb,
    label& nPoints
)
{
    // Unfortunately nPoints constructs meshPoints() so do compact version
    // ourselves

    const auto& pts = p.points();

    bitSet pointUsed(pts.size());

    nPoints = 0;
    bb.reset();

    for (const auto& f : p)
    {
        for (const label pointi : f)
        {
            if (pointUsed.set(pointi))
            {
                bb.add(pts[pointi]);
                ++nPoints;
            }
        }
    }
}


// ************************************************************************* //
