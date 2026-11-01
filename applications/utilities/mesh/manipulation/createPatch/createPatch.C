/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    createPatch

Group
    grpMeshManipulationUtilities

Description
    Create patches out of selected boundary faces, which are either
    from existing patches or from a faceSet.

    More specifically it:
    - creates new patches (from selected boundary faces).
      Synchronise faces on coupled patches.
    - synchronises points on coupled boundaries
    - remove patches with 0 faces in them

    Mapping fields:
    - always : map within a patch (default)
    - if possible : map from old patch to new patch (if we're just renaming
      patches)
    - if combining patches : apply the same mapping as the patch faces
    - if splitting patches : reverse the mapping used to combine the patches

    Combining/splitting patches:
    If used to merge patches into a single patch ('constructFrom patches'), it:
    - creates faceZones from patches-to-merge
    - takes over patchGroup info into faceZoneGroup
    - creates faceZoneGroup 'createPatch:<patchName>' to hold all original
      patch/faceZones
    - maps (in faceZone order) the original patch contents into the new
      patch.

    It can be run in reverse by supplying the -undo option which will
    split the patch back into original components (by using the faceZone)
    information.

    OK - have option to avoid creating faceZones
       'allowUndo true;'
    OK - if constructFromSet : clone patch from min patch in set. Map fields
      as if constructFromPatches.
    - if constructFromSet : allow undoing. Any set can only contain boundary
      faces so we can always construct a equivalent set of (partial) patches
      to map from.

        createPatchDict:
            name            combined_patch;
            allowUndo       true;
            constructFrom   set;
            name            some_boundary_faces;

        with faceSet 'some_boundary_faces' containing faces from patches A, B, C

      So:
      - createPatch will create a new patch 'combined_patch' with the selected
        faces from A, B, C
      - if allowUndo:
        - createPatch will create a faceZone for each of A, B, C
        - createPatch will create a faceZoneGroup 'createPatch:combined_patch'
          containing the faceZones for A, B, C
      - if the A,B,C become empty then they will be removed from the mesh
      - if actuallly -undo
        - use the 'combined_patch' faceZoneGroup to find out the original patch
          names A,B,C
        - recreate the patches A, B, C (if not present)
        - move the faces from the faceZones given by 'combined_patch'
          back into patches A, B, C
    - if 'patchFields' present warn that no undoing is possible. Or only warn
      if -undo actually used.
    OK - the number - index of patch seems missing for patches which don't have
      faces on processor0 (eg on motorBike it says
      'Mapping fields on patch motorBike_rr-wh-rim:005%5 from patches ()')
    OK - the -undo operation will need -lib libturbulenceModels.so
      since wallFunctions are not known to createPatch.
    - undo still check against allowUndo in createPatchDict.
      Not undo if allowUndo false.
    - patchInfo seems to trigger getting calculated bcs which then conflict
      with existing bcs.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"
#include "IOPtrList.H"
#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "wordRes.H"
#include "processorMeshes.H"
#include "IOdictionary.H"
#include "regionProperties.H"
#include "faceAreaWeightAMI2D.H"
#include "fvMeshTools.H"
#include "ReadFields.H"
#include "fvFieldReconstructor.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
}


word patchName(const word& name, const fvMesh& mesh0, const fvMesh& mesh1)
{
    word pName(name != "none" ? name : word::null);
    pName += mesh0.name();
    pName += "_to_";
    pName += mesh1.name();
    return pName;
}


void matchPatchFaces
(
    const word& entryName,
    const word& AMIMethod,
    const dictionary& AMIDict,
    const PtrList<fvMesh>& meshes,

    const label meshi,
    const label nSourcei,
    const label sourcei,
    const labelList& patchesi,

    const label meshj,
    const label nSourcej,
    const label sourcej,
    const labelList& patchesj,

    DynamicList<labelList>& interfaceMesh0,
    DynamicList<label>& interfaceSource0,
    DynamicList<labelList>& interfacePatch0,
    DynamicList<wordList>& interfaceNames0,
    DynamicList<List<DynamicList<label>>>& interfaceFaces0,

    DynamicList<labelList>& interfaceMesh1,
    DynamicList<label>& interfaceSource1,
    DynamicList<labelList>& interfacePatch1,
    DynamicList<wordList>& interfaceNames1,
    DynamicList<List<DynamicList<label>>>& interfaceFaces1
)
{
    // Now we have:
    // - meshi, sourcei, patchesi
    // - meshj, sourcej, patchesj


    // Attempt to match patches
    forAll(patchesi, i)
    {
        const auto& ppi = meshes[meshi].boundaryMesh()[patchesi[i]];

        forAll(patchesj, j)
        {
            const auto& ppj = meshes[meshj].boundaryMesh()[patchesj[j]];

            // Use AMI to try and find matches
            auto AMPtr(AMIInterpolation::New(AMIMethod, AMIDict));

            AMPtr->calculate(ppi, ppj, nullptr);
            if
            (
                gAverage(AMPtr->tgtWeightsSum()) > SMALL
             || gAverage(AMPtr->srcWeightsSum()) > SMALL
            )
            {
                const label inti = interfaceMesh0.size();

                Info<< "Introducing interface " << inti << " between"
                    << " mesh " << meshes[meshi].name()
                    //<< " source:" << sourcei
                    << " patch " << ppi.name()
                    << " and mesh " << meshes[meshj].name()
                    //<< " source:" << sourcej
                    << " patch " << ppj.name()
                    << endl;

                // Mesh 0
                //~~~~~~~

                auto& intMesh0 = interfaceMesh0.emplace_back(nSourcei, -1);
                intMesh0[sourcei] = meshi;

                interfaceSource0.push_back(sourcei);

                auto& intPatch0 = interfacePatch0.emplace_back(nSourcei, -1);
                intPatch0[sourcei] = ppi.index();

                auto& intNames0 = interfaceNames0.emplace_back(nSourcei);
                intNames0[sourcei] =
                    patchName(entryName, meshes[meshi], meshes[meshj]);


                // Mesh 1
                //~~~~~~~

                auto& intMesh1 = interfaceMesh1.emplace_back(nSourcej, -1);
                intMesh1[sourcej] = meshj;

                interfaceSource1.push_back(sourcej);

                auto& intPatch1 = interfacePatch1.emplace_back(nSourcej, -1);
                intPatch1[sourcej] = ppj.index();

                auto& intNames1 = interfaceNames1.emplace_back(nSourcej);
                intNames1[sourcej] =
                    patchName(entryName, meshes[meshj], meshes[meshi]);

                auto& intFaces0 = interfaceFaces0.emplace_back(nSourcei);
                DynamicList<label>& faces0 = intFaces0[sourcei];
                faces0.setCapacity(ppi.size());

                auto& intFaces1 = interfaceFaces1.emplace_back(nSourcej);
                DynamicList<label>& faces1 = intFaces1[sourcej];
                faces1.setCapacity(ppj.size());


                // Mark any mesh faces:
                // - that have a weight > 0.5.
                // - contribute to these faces (even if < 0.5)

                faces0.clear();
                faces1.clear();

                // Marked as donor of face with high weight
                scalarField targetMask;
                {
                    const scalarField& weights = AMPtr->srcWeightsSum();
                    scalarField mask(weights.size(), Zero);
                    forAll(weights, facei)
                    {
                        const scalar sum = weights[facei];
                        if (sum > 0.5)
                        {
                            mask[facei] = sum;
                        }
                    }

                    // Push source field mask to target
                    targetMask = AMPtr->interpolateToTarget(mask);
                }

                scalarField sourceMask;
                {
                    const scalarField& weights = AMPtr->tgtWeightsSum();
                    scalarField mask(weights.size(), Zero);
                    forAll(weights, facei)
                    {
                        const scalar sum = weights[facei];
                        if (sum > 0.5)
                        {
                            mask[facei] = sum;
                        }
                    }

                    // Push target mask back to source
                    sourceMask = AMPtr->interpolateToSource(mask);
                }

                {
                    const scalarField& weights = AMPtr->srcWeightsSum();
                    forAll(weights, facei)
                    {
                        if (weights[facei] > 0.5 || sourceMask[facei] > SMALL)
                        {
                            faces0.push_back(ppi.start()+facei);
                        }
                    }
                }
                {
                    const scalarField& weights = AMPtr->tgtWeightsSum();
                    forAll(weights, facei)
                    {
                        if (weights[facei] > 0.5 || targetMask[facei] > SMALL)
                        {
                            faces1.push_back(ppj.start()+facei);
                        }
                    }
                }

                faces0.shrink();
                faces1.shrink();
            }
        }
    }
}


void matchPatchFaces
(
    const bool includeOwn,
    const PtrList<fvMesh>& meshes,
    const List<DynamicList<label>>& interRegionSources,
    const List<wordList>& patchNames,
    const labelListListList& matchPatchIDs,
    const List<wordList>& AMIMethods,
    List<PtrList<dictionary>> patchInfoDicts,

    DynamicList<labelList>& interfaceMesh0,
    DynamicList<label>& interfaceSource0,
    DynamicList<labelList>& interfacePatch0,
    DynamicList<List<DynamicList<label>>>& interfaceFaces0,
    DynamicList<wordList>& interfaceNames0,

    DynamicList<labelList>& interfaceMesh1,
    DynamicList<label>& interfaceSource1,
    DynamicList<labelList>& interfacePatch1,
    DynamicList<List<DynamicList<label>>>& interfaceFaces1,
    DynamicList<wordList>& interfaceNames1
)
{
    // Add approximate matches
    forAll(meshes, meshi)
    {
        const labelListList& allPatchesi = matchPatchIDs[meshi];

        const label meshjStart(includeOwn ? meshi : meshi+1);

        for (label meshj = meshjStart; meshj < meshes.size(); meshj++)
        {
            const labelListList& allPatchesj = matchPatchIDs[meshj];

            for (const label sourcei : interRegionSources[meshi])
            {
                const labelList& patchesi = allPatchesi[sourcei];

                for (const label sourcej : interRegionSources[meshj])
                {
                    const labelList& patchesj = allPatchesj[sourcej];

                    // Now we have:
                    // - meshi, sourcei, patchesi
                    // - meshj, sourcej, patchesj

                    matchPatchFaces
                    (
                        patchNames[meshi][sourcei],
                        (
                            AMIMethods[meshi][sourcei].size()
                          ? AMIMethods[meshi][sourcei]
                          : faceAreaWeightAMI2D::typeName
                        ),
                        patchInfoDicts[meshi][sourcei],
                        meshes,

                        meshi,
                        allPatchesi.size(), // nSourcei,
                        sourcei,
                        patchesi,

                        meshj,
                        allPatchesj.size(), // nSourcej,
                        sourcej,
                        patchesj,

                        interfaceMesh0,
                        interfaceSource0,
                        interfacePatch0,
                        interfaceNames0,
                        interfaceFaces0,

                        interfaceMesh1,
                        interfaceSource1,
                        interfacePatch1,
                        interfaceNames1,
                        interfaceFaces1
                    );
                }
            }
        }
    }
}


void changePatchID
(
    const bool modify,
    const fvMesh& mesh,
    const label faceID,
    const label patchID,
    const label zoneID,     // Optional additional zone
    const bool addZone,
    polyTopoChange& meshMod
)
{
    const auto& faceZones = mesh.faceZones();

    DynamicList<label> zoneIDs;
    DynamicList<bool> zoneFlips;

    faceZones.whichZones(faceID, zoneIDs);
    for (const label zonei : zoneIDs)
    {
        const label index = faceZones[zonei].localID(faceID);
        zoneFlips.append(faceZones[zonei].flipMap()[index]);
    }
    // Optional additional zone
    if (zoneID != -1)
    {
        if (addZone && !zoneIDs.found(zoneID))
        {
            zoneIDs.append(zoneID);
            zoneFlips.append(false);
        }
        else if (!addZone && zoneIDs.found(zoneID))
        {
            const label index = zoneIDs.find(zoneID);
            zoneIDs.remove(index);
            zoneFlips.remove(index);
        }
    }

    if (modify)
    {
        meshMod.modifyFace
        (
            mesh.faces()[faceID],               // face
            faceID,                             // face ID
            mesh.faceOwner()[faceID],           // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            zoneIDs,                            // zone ID
            zoneFlips                           // zone flip
        );
    }
    else
    {
        meshMod.addFace
        (
            mesh.faces()[faceID],       // modified face
            mesh.faceOwner()[faceID],   // owner
            -1,                         // neighbour
            -1,                         // master point
            -1,                         // master edge
            faceID,                     // master face
            false,                      // face flip
            patchID,                    // patch for face
            zoneIDs,                    // zone for face
            zoneFlips                   // face flip in zone
        );
    }
}


void changePatchID
(
    const fvMesh& mesh,
    const labelUList& faceLabels,
    const label patchID,
    const label zoneID,
    const bool addZone,
    bitSet& isRepatchedBoundary,
    polyTopoChange& meshMod
)
{
    for (const label facei : faceLabels)
    {
        if (mesh.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Face " << facei
                << " is not an external face of the mesh." << endl
                << "This application can only repatch"
                << " existing boundary faces." << exit(FatalError);
        }

        const bool isFirst =
            isRepatchedBoundary.set(facei-mesh.nInternalFaces());
        if (!isFirst)
        {
            static label nWarnings = 0;
            if (nWarnings == 0)
            {
                const label newPatchi = meshMod.region()[facei];
                //FatalErrorInFunction
                WarningInFunction
                    << "Face " << facei
                    << " at " << mesh.faceCentres()[facei]
                    << " marked for patch " << patchID
                    << " name " << mesh.boundaryMesh()[patchID].name()
                    << " is already marked for patch " << newPatchi
                    << " name " << mesh.boundaryMesh()[newPatchi].name()
                    << ". Creating duplicate face. Suppressing further warnings"
                    //<< exit(FatalError);
                    << endl;
            }
            nWarnings++;
        }

        changePatchID(isFirst, mesh, facei, patchID, zoneID, addZone, meshMod);
    }
}


void changePatchID
(
    const bool allowUndo,
    const fvMesh& mesh,
    const labelList& sortedFaceLabels,
    const label destPatchi,
    bitSet& isRepatchedBoundary,
    polyTopoChange& meshMod
)
{
    const auto& patches = mesh.boundaryMesh();
    const auto& fzs = mesh.faceZones();

    if (sortedFaceLabels.empty())
    {
        return;
    }

    label patchStart = 0;
    const label facei = sortedFaceLabels[0];
    label curPatchi = patches.whichPatch(facei);
    forAll(sortedFaceLabels, i)
    {
        const label facei = sortedFaceLabels[i];
        const label patchi = patches.whichPatch(facei);
        if
        (
            patchi != curPatchi
         || i == sortedFaceLabels.size()-1
        )
        {
            // Change of patch or last face. Change all faces from patchStart
            // to i-1

            // Add to optional backup faceZone
            const label zonei =
            (
                allowUndo
              ? fzs.findZoneID(patches[curPatchi].name())
              : -1
            );
            const label nFaces =
            (
                i == sortedFaceLabels.size()-1
              ? i-patchStart+1
              : i-patchStart
            );

            const SubList<label> slice(sortedFaceLabels, nFaces, patchStart);

            changePatchID
            (
                mesh,
                slice,
                destPatchi,
                zonei,  // backup faceZone
                true,   // preserve faceZones
                isRepatchedBoundary,
                meshMod
            );
            curPatchi = patchi;
            patchStart = i;
        }
    }
}


// Dump for all patches the current match
void dumpCyclicMatch(const fileName& prefix, const polyMesh& mesh)
{
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            // Dump patches
            {
                OFstream str(prefix+cycPatch.name()+".obj");
                Pout<< "Dumping " << cycPatch.name()
                    << " faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    cycPatch,
                    cycPatch.points()
                );
            }

            {
                OFstream str(prefix+nbrPatch.name()+".obj");
                Pout<< "Dumping " << nbrPatch.name()
                    << " faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    nbrPatch,
                    nbrPatch.points()
                );
            }


            // Lines between corresponding face centres
            OFstream str(prefix+cycPatch.name()+nbrPatch.name()+"_match.obj");
            label vertI = 0;

            Pout<< "Dumping cyclic match as lines between face centres to "
                << str.name() << endl;

            forAll(cycPatch, facei)
            {
                const point& fc0 = mesh.faceCentres()[cycPatch.start()+facei];
                meshTools::writeOBJ(str, fc0);
                vertI++;
                const point& fc1 = mesh.faceCentres()[nbrPatch.start()+facei];
                meshTools::writeOBJ(str, fc1);
                vertI++;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
}


void separateList
(
    const vectorField& separation,
    UList<vector>& field
)
{
    if (separation.size() == 1)
    {
        // Single value for all.

        forAll(field, i)
        {
            field[i] += separation[0];
        }
    }
    else if (separation.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] += separation[i];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}


// Synchronise points on both sides of coupled boundaries.
template<class CombineOp>
void syncPoints
(
    const polyMesh& mesh,
    pointField& points,
    const CombineOp& cop,
    const point& nullValue
)
{
    if (points.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << points.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (UPstream::parRun())
    {
        const labelList& procPatches = mesh.globalData().processorPatches();

        // Send
        for (const label patchi : procPatches)
        {
            const polyPatch& pp = patches[patchi];
            const auto& procPatch = refCast<const processorPolyPatch>(pp);

            if (pp.nPoints() && procPatch.owner())
            {
                // Get data per patchPoint in neighbouring point numbers.
                pointField patchInfo(procPatch.nPoints(), nullValue);

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                forAll(nbrPts, pointi)
                {
                    label nbrPointi = nbrPts[pointi];
                    if (nbrPointi >= 0 && nbrPointi < patchInfo.size())
                    {
                        patchInfo[nbrPointi] = points[meshPts[pointi]];
                    }
                }

                OPstream toNbr
                (
                    UPstream::commsTypes::buffered,
                    procPatch.neighbProcNo()
                );
                toNbr << patchInfo;
            }
        }


        // Receive and set.

        for (const label patchi : procPatches)
        {
            const polyPatch& pp = patches[patchi];
            const auto& procPatch = refCast<const processorPolyPatch>(pp);

            if (pp.nPoints() && !procPatch.owner())
            {
                // We do not know the number of points on the other side
                // so cannot use UIPstream::read

                pointField nbrPatchInfo;
                IPstream::recv(nbrPatchInfo, procPatch.neighbProcNo());

                // Null any value which is not on neighbouring processor
                nbrPatchInfo.setSize(procPatch.nPoints(), nullValue);

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                const labelList& meshPts = procPatch.meshPoints();

                forAll(meshPts, pointi)
                {
                    label meshPointi = meshPts[pointi];
                    points[meshPointi] = nbrPatchInfo[pointi];
                }
            }
        }
    }

    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();
            const labelList& nbrMeshPts = nbrPatch.meshPoints();

            pointField half0Values(coupledPoints.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point0 = meshPts[e[0]];
                half0Values[i] = points[point0];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.reverseT(), half0Values);
            }
            else if (cycPatch.separated())
            {
                hasTransformation = true;
                separateList(cycPatch.separation(), half0Values);
            }

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point1 = nbrMeshPts[e[1]];
                points[point1] = half0Values[i];
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly necessary.
    //Pstream::reduceOr(hasTransformation);

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningInFunction
                << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }


        // Values on shared points.
        pointField sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = points[meshPointi];
        }

        // Combine - globally consistent
        Pstream::listCombineReduce(sharedPts, cop);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            points[meshPointi] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}


Foam::label createPatch
(
    fvMesh& mesh,
    const word& patchName,
    const label sourcePatchi,    // -1 or index of patch to clone
    const dictionary& patchDict, // If not cloning, dict with patch definition
    const word& groupName        // Optional group for new patch (for undoing)
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label newPatchi = patches.findPatchID(patchName);
    if (newPatchi == -1)
    {
        autoPtr<polyPatch> ppPtr;
        if (sourcePatchi != -1)
        {
            Info<< "Creating new patch " << patchName
                << " by cloning "
                << patches[sourcePatchi].name() << endl;
            ppPtr = patches[sourcePatchi].clone(patches);
            ppPtr->name() = patchName;
        }
        else
        {
            Info<< "Creating new patch " << patchName
                << " from " << patchDict
                << endl;

            dictionary allDict(patchDict);
            allDict.set("nFaces", 0);
            allDict.set("startFace", 0);
            ppPtr = polyPatch::New
            (
                patchName,
                allDict,
                0,          // index, overwritten
                patches
            );
        }
        // Put new patch in combined group for undoing.
        if (!groupName.empty() && !ppPtr().inGroups().found(groupName))
        {
            ppPtr().inGroups().append(groupName);
        }
        // Add patch to mesh & fields
        newPatchi = fvMeshTools::addPatch
        (
            mesh,
            ppPtr(),
            patchDict.subOrEmptyDict("patchFields"),
            fvPatchFieldBase::calculatedType(),
            true,
            sourcePatchi    // source patch to clone fields from
        );
    }
    return newPatchi;
}


template<class GeoField>
void storeBFields
(
    const PtrList<fvMesh>& meshes,
    const List<PtrList<GeoField>>& flds,
    List<PtrList<typename GeoField::Boundary>>& savedFlds
)
{
    savedFlds.setSize(flds.size());
    forAll(meshes, meshi)
    {
        auto& savedBflds = savedFlds[meshi];
        savedBflds.setSize(flds[meshi].size());
        forAll(flds[meshi], fldi)
        {
            const auto& fld = flds[meshi][fldi];
            Info<< "Saving boundary field of " << fld.name()
                << endl;
            const auto& ifld = fld.internalField();
            const auto& bfld = fld.boundaryField();
            savedBflds.set(fldi, new typename GeoField::Boundary(ifld, bfld));
        }
    }
}


template<class GeoField>
void mapBFields
(
    const label patchi,
    const DynamicList<label>& oldPatchIDs,
    const labelListList& newPatchFaces,
    const PtrList<typename GeoField::Boundary>& oldBflds,
    PtrList<GeoField>& flds
)
{
    // Combine old patch fields. Like reconstructPar / fvFieldReconstructor.
    // newPatchFaces : for each old patch, each old patch face, the face in
    // the new patch. (so usable in rmap)

    if (oldPatchIDs.empty())
    {
        return;
    }

    forAll(oldBflds, fldi)
    {
        const auto& oldBfld = oldBflds[fldi];
        typename GeoField::Boundary& bfld = flds[fldi].boundaryFieldRef();

        if (!bfld.set(patchi))
        {
            // Construct new patchField. Type and size only. See
            // fvFieldReconstructor.txx
            bfld.set
            (
                patchi,
                GeoField::Patch::New
                (
                    oldBfld[oldPatchIDs[0]],
                    bfld[patchi].patch(),
                    bfld[patchi].internalField(),
                    fvFieldReconstructor::fvPatchFieldReconstructor
                    (
                        bfld[patchi].patch().size()
                    )
                )
            );
        }

        // Now we either have a new patchField or an existing one.
        auto& newPfld = bfld[patchi];

        for (const label oldPatchi : oldPatchIDs)
        {
            const auto& oldPfld = oldBfld[oldPatchi];

            // Pout<< "Mapping patch field " << flds[fldi].name()
            //     << " for patch " << newPfld.patch().name()
            //     << " on combined patch:" << patchi
            //     << " from old patch " << oldPatchi
            //     << " old type:" << oldBfld[oldPatchi].type()
            //     << " old values:" << oldPfld
            //     << nl
            //     << " with map:" << flatOutput(newPatchFaces[oldPatchi]) << endl;

            if (oldPfld.type() != newPfld.type())
            {
                WarningInFunction
                    << "For field " << flds[fldi].name()
                    << " when trying to combine patches " << oldPatchIDs
                    << " into single patch " << newPfld.patch().name()
                    << " the original field type " << oldPfld.type()
                    << " differs from the new type " << newPfld.type()
                    << " so we cannot map the patch fields."
                    << endl;
            }
            else
            {
                // Use rmap for mapping. This requires the map to be the
                // size as the source field and have -1 for any unused source
                // value.
                const auto& map = newPatchFaces[oldPatchi];
                labelList expandedMap(oldPfld.size(), -1);
                SubList<label>(expandedMap, map.size()) = map;
                newPfld.rmap(oldPfld, expandedMap);
            }
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create patches out of selected boundary faces, which are either"
        " from existing patches or from a faceSet"
    );

    #include "addOverwriteOption.H"
    #include "addAllRegionOptions.H"

    argList::addOption("dict", "file", "Alternative createPatchDict");
    argList::addBoolOption
    (
        "undo",
        "Undo the last createPatch operation"
        " (split patches back into original components)"
    );
    argList::addBoolOption
    (
        "writeObj",
        "Write obj files showing the cyclic matching process"
    );
    argList::addBoolOption
    (
        "no-fields",
        "Suppress conversion of fields (volume, finite-area, lagrangian)"
    );


    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "getAllRegionOptions.H"

    const bool overwrite = args.found("overwrite");

    #include "createNamedMeshes.H"

    // Disable reading/writing fields. Not compatible with patchFields entries
    const bool noFields = args.found("no-fields");
    const bool undo = args.found("undo");
    const bool writeObj = args.found("writeObj");


    // Read dictionaries and extract various
    wordList oldInstances(meshes.size());
    PtrList<dictionary> dicts(meshes.size());
    List<wordList> patchNames(meshes.size());
    List<labelListList> matchPatchIDs(meshes.size());
    List<wordList> matchMethods(meshes.size());
    List<PtrList<dictionary>> patchInfoDicts(meshes.size());

    List<DynamicList<label>> interRegionSources(meshes.size());
    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // If running parallel check same patches everywhere
        patches.checkParallelSync(true);

        oldInstances[meshi] = mesh.pointsInstance();

        const word dictName("createPatchDict");
        #include "setSystemMeshDictionaryIO.H"
        Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;

        dicts.set(meshi, new IOdictionary(dictIO));

        const auto& dict = dicts[meshi];
        PtrList<dictionary> patchSources(dict.lookup("patches"));
        patchNames[meshi].setSize(patchSources.size());
        matchPatchIDs[meshi].setSize(patchSources.size());
        matchMethods[meshi].setSize(patchSources.size());
        patchInfoDicts[meshi].setSize(patchSources.size());

        // Read patch construct info from dictionary
        forAll(patchSources, sourcei)
        {
            const auto& pDict = patchSources[sourcei];
            patchNames[meshi][sourcei] = pDict.getOrDefault<word>
            (
                "name",
                word::null,
                keyType::LITERAL
            );

            patchInfoDicts[meshi].set
            (
                sourcei,
                new dictionary(pDict.subOrEmptyDict("patchInfo"))
            );

            dictionary& patchDict = patchInfoDicts[meshi][sourcei];
            if (patchDict.found("AMIMethod"))
            {
                matchMethods[meshi][sourcei] = patchDict.get<word>("AMIMethod");
                // Disable full matching since we're trying to use AMIMethod to
                // find out actual overlap
                patchDict.add("requireMatch", false);
            }

            wordRes matchNms;
            if (!undo)
            {
                if (pDict.readIfPresent("patches", matchNms, keyType::LITERAL))
                {
                    matchPatchIDs[meshi][sourcei] =
                        patches.patchSet(matchNms).sortedToc();
                }
                else if (pDict.get<word>("constructFrom") == "set")
                {
                    // If constructing from a set use any patch in the set
                    // to clone the patch from.
                    const word setName = pDict.get<word>("set");
                    const faceSet set(mesh, setName);

                    Info<< "Read " << returnReduce(set.size(), sumOp<label>())
                        << " faces from faceSet " << set.name() << endl;
                    // label minPatchi = labelMax;
                    for (const label facei : set)
                    {
                        const label patchi = patches.whichPatch(facei);
                        if (patchi == -1)
                        {
                            FatalErrorInFunction
                                << "Face " << facei
                                << " is not a boundary face."
                                << exit(FatalError);
                        }
                        // minPatchi = min(minPatchi, patchi);
                        matchPatchIDs[meshi][sourcei].append(patchi);
                    }
                    // minPatchi = returnReduce(minPatchi, minOp<label>());
                    // matchPatchIDs[meshi][sourcei].setSize(1, minPatchi);
                }
            }

            if (pDict.get<word>("constructFrom") == "autoPatch")
            {
                interRegionSources[meshi].push_back(sourcei);
            }
        }
    }



    // Determine matching faces. This is an interface between two
    // sets of faces:
    //  - originating mesh
    //  - originating patch
    //  - originating (mesh)facelabels
    // It matches all mesh against each other. Lower numbered mesh gets
    // postfix 0, higher numbered mesh postfix 1.

    // Per interface, per patchSource:
    //      1. the lower numbered mesh
    DynamicList<labelList> interfaceMesh0;
    //      1b. the source index (i.e. the patch dictionary)
    DynamicList<label> interfaceSource0;
    //      2. the patch on the interfaceMesh0
    DynamicList<labelList> interfacePatch0;
    //      3. the facelabels on the interfaceMesh0
    DynamicList<List<DynamicList<label>>> interfaceFaces0;
    //      4. generated interface name
    DynamicList<wordList> interfaceNames0;

    // Same for the higher numbered mesh
    DynamicList<labelList> interfaceMesh1;
    DynamicList<label> interfaceSource1;
    DynamicList<labelList> interfacePatch1;
    DynamicList<List<DynamicList<label>>> interfaceFaces1;
    DynamicList<wordList> interfaceNames1;
    if (!undo)
    {
        // Whether to match to patches in own mesh
        const bool includeOwn = (meshes.size() == 1);

        matchPatchFaces
        (
            includeOwn,
            meshes,
            interRegionSources,
            patchNames,
            matchPatchIDs,
            matchMethods,
            patchInfoDicts,

            interfaceMesh0,
            interfaceSource0,
            interfacePatch0,
            interfaceFaces0,
            interfaceNames0,

            interfaceMesh1,
            interfaceSource1,
            interfacePatch1,
            interfaceFaces1,
            interfaceNames1
        );
    }



    // Read fields
    List<PtrList<volScalarField>> vsFlds(meshes.size());
    List<PtrList<volVectorField>> vvFlds(meshes.size());
    List<PtrList<volSphericalTensorField>> vstFlds(meshes.size());
    List<PtrList<volSymmTensorField>> vsymtFlds(meshes.size());
    List<PtrList<volTensorField>> vtFlds(meshes.size());

    List<PtrList<surfaceScalarField>> ssFlds(meshes.size());
    List<PtrList<surfaceVectorField>> svFlds(meshes.size());
    List<PtrList<surfaceSphericalTensorField>> sstFlds(meshes.size());
    List<PtrList<surfaceSymmTensorField>> ssymtFlds(meshes.size());
    List<PtrList<surfaceTensorField>> stFlds(meshes.size());

    forAll(meshes, meshi)
    {
        const fvMesh& mesh = meshes[meshi];

        if (noFields)
        {
            for (const auto& d : patchInfoDicts[meshi])
            {
                if (d.found("patchFields"))
                {
                    FatalIOErrorInFunction(dicts[meshi])
                        << "Cannot read patchFields when no-fields option"
                        << " is specified."
                        << exit(FatalIOError);
                }
            }
        }
        if (!noFields)
        {
            // Read objects in time directory
            IOobjectList objects(mesh, runTime.timeName());

            // Read volume fields.
            ReadFields(mesh, objects, vsFlds[meshi]);
            ReadFields(mesh, objects, vvFlds[meshi]);
            ReadFields(mesh, objects, vstFlds[meshi]);
            ReadFields(mesh, objects, vsymtFlds[meshi]);
            ReadFields(mesh, objects, vtFlds[meshi]);

            // Read surface fields.
            ReadFields(mesh, objects, ssFlds[meshi]);
            ReadFields(mesh, objects, svFlds[meshi]);
            ReadFields(mesh, objects, sstFlds[meshi]);
            ReadFields(mesh, objects, ssymtFlds[meshi]);
            ReadFields(mesh, objects, stFlds[meshi]);
        }
    }


    const label nMesh = meshes.size();

    // Store old patch fields (used when mapping after topo changes)
    List<PtrList<volScalarField::Boundary>> oldVss(nMesh);
    List<PtrList<volVectorField::Boundary>> oldVvs(nMesh);
    List<PtrList<volSphericalTensorField::Boundary>> oldVsts(nMesh);
    List<PtrList<volSymmTensorField::Boundary>> oldVsyms(nMesh);
    List<PtrList<volTensorField::Boundary>> oldVts(nMesh);

    storeBFields(meshes, vsFlds, oldVss);
    storeBFields(meshes, vvFlds, oldVvs);
    storeBFields(meshes, vstFlds, oldVsts);
    storeBFields(meshes, vsymtFlds, oldVsyms);
    storeBFields(meshes, vtFlds, oldVts);

    List<PtrList<surfaceScalarField::Boundary>> oldSsBflds(nMesh);
    List<PtrList<surfaceVectorField::Boundary>> oldSvBflds(nMesh);
    List<PtrList<surfaceSphericalTensorField::Boundary>> oldSstBflds(nMesh);
    List<PtrList<surfaceSymmTensorField::Boundary>> oldSsymtBflds(nMesh);
    List<PtrList<surfaceTensorField::Boundary>> oldStBflds(nMesh);

    storeBFields(meshes, ssFlds, oldSsBflds);
    storeBFields(meshes, svFlds, oldSvBflds);
    storeBFields(meshes, sstFlds, oldSstBflds);
    storeBFields(meshes, ssymtFlds, oldSsymtBflds);
    storeBFields(meshes, stFlds, oldStBflds);


    // Maintain list of added patches so we exclude them from filtering
    // later on
    List<DynamicList<word>> allAddedPatches(meshes.size());

    // Loop over all regions

    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];
        auto& fzs = mesh.faceZones();

        Info<< "\n\nAdding patches to mesh " << mesh.name() << nl << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const dictionary& dict = dicts[meshi];

        if (writeObj)
        {
            dumpCyclicMatch("initial_", mesh);
        }

        // Read patch construct info from dictionary
        PtrList<dictionary> patchSources(dict.lookup("patches"));


        // 1. Add all new patches
        // ~~~~~~~~~~~~~~~~~~~~~~

        forAll(patchSources, sourcei)
        {
            const dictionary& dict = patchSources[sourcei];
            const word sourceType(dict.get<word>("constructFrom"));
            const word& patchName = patchNames[meshi][sourcei];
            const word zoneGroupName("createPatch:" + patchName);
            const bool allowUndo = dict.getOrDefault("allowUndo", false);

            const dictionary& patchDict = patchInfoDicts[meshi][sourcei];


            if (sourceType == "autoPatch")
            {
                // Special : automatically create the necessary inter-region
                // coupling patches

                if (undo)
                {
                    FatalErrorInFunction
                        << "Cannot undo autoPatch construction. Aborting."
                        << exit(FatalError);
                }

                forAll(interfaceMesh0, inti)
                {
                    const labelList& allMeshes0 = interfaceMesh0[inti];
                    const wordList& allNames0 = interfaceNames0[inti];
                    const labelList& allMeshes1 = interfaceMesh1[inti];
                    const wordList& allNames1 = interfaceNames1[inti];

                    if
                    (
                        interfaceSource0[inti] == sourcei
                     && allMeshes0[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh0. mesh1 is the remote mesh.
                        const label sourcej = interfaceSource1[inti];
                        const word& patchName = allNames0[sourcei];
                        if (patches.findPatchID(patchName) == -1)
                        {
                            dictionary allDict(patchDict);
                            const auto& mesh1 = meshes[allMeshes1[sourcej]];
                            allDict.set("sampleRegion", mesh1.name());
                            const auto& destPatch = allNames1[sourcej];
                            allDict.set("samplePatch", destPatch);
                            allDict.set("neighbourPatch", destPatch);

                            Info<< "Adding new patch " << patchName
                                << " from " << allDict << endl;
                            createPatch
                            (
                                mesh,
                                patchName,
                                -1,
                                allDict,
                                word::null
                            );
                            allAddedPatches[meshi].append(patchName);
                        }
                    }
                }

                forAll(interfaceMesh1, inti)
                {
                    const labelList& allMeshes0 = interfaceMesh0[inti];
                    const wordList& allNames0 = interfaceNames0[inti];
                    const labelList& allMeshes1 = interfaceMesh1[inti];
                    const wordList& allNames1 = interfaceNames1[inti];

                    if
                    (
                        interfaceSource1[inti] == sourcei
                     && allMeshes1[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh1. mesh0 is the remote mesh.

                        const label sourcej = interfaceSource0[inti];
                        const word& patchName = allNames1[sourcei];
                        if (patches.findPatchID(patchName) == -1)
                        {
                            dictionary allDict(patchDict);
                            const auto& destPatch = allNames0[sourcej];
                            const auto& mesh0 = meshes[allMeshes0[sourcej]];
                            allDict.set("sampleRegion", mesh0.name());
                            allDict.set("samplePatch", destPatch);
                            allDict.set("neighbourPatch", destPatch);

                            Info<< "Adding new patch " << patchName
                                << " from " << allDict << endl;
                            createPatch
                            (
                                mesh,
                                patchName,
                                -1,
                                allDict,
                                word::null
                            );
                            allAddedPatches[meshi].append(patchName);
                        }
                    }
                }
            }
            else
            {
                if (undo)
                {
                    // Note: fvMeshTools::addPatch clears out zoneGroup info
                    //       so take copy instead of using info directly.
                    const auto patchSources = fzs.groupZoneIDs()[zoneGroupName];

                    // Ignore any patchInfo overrides. TBD.
                    Info<< "Splitting patch " << patchName << endl;
                    const label patchi = patches.findPatchID(patchName, false);
                    for (const label zonei : patchSources)
                    {
                        const auto& fz = fzs[zonei];
                        const label newPatchi = createPatch
                        (
                            mesh,
                            fz.name(),
                            patchi,             // patch to clone
                            patchDict,          // optional patch fields
                            word::null
                        );

                        const polyPatch& pp = patches[newPatchi];
                        wordList& gNames = const_cast<wordList&>(pp.inGroups());

                        // Append faceZone groups to patch groups
                        for (const word& groupName : fz.inGroups())
                        {
                            if
                            (
                                groupName != zoneGroupName
                            && !gNames.found(groupName)
                            )
                            {
                                gNames.append(groupName);
                            }
                        }
                        allAddedPatches[meshi].append(pp.name());
                    }
                }
                else
                {
                    const labelList& patchSources =
                        matchPatchIDs[meshi][sourcei];

                    // If in patch-combine mode add old patches as faceZones
                    if (allowUndo)
                    {
                        fzs.clearAddressing();

                        for (const label patchi : patchSources)
                        {
                            const polyPatch& pp = patches[patchi];
                            if (fzs.findZoneID(pp.name()) == -1)
                            {
                                const label index = fzs.size();
                                Info<< "Adding faceZone " << pp.name()
                                    << " to store original faces of patch "
                                    << pp.name() << endl;
                                fzs.append(new faceZone(pp.name(), index, fzs));

                                wordList& gNames = fzs.last().inGroups();
                                // Append pp groups to faceZone groups
                                fzs.clearAddressing();
                                forAll(pp.inGroups(), i)
                                {
                                    const word& groupName = pp.inGroups()[i];
                                    if (!gNames.found(groupName))
                                    {
                                        gNames.append(groupName);
                                    }
                                }
                                // Put new faceZone in combined group for
                                // undoing.
                                if (!gNames.found(zoneGroupName))
                                {
                                    gNames.append(zoneGroupName);
                                }
                            }
                        }
                    }

                    if (patches.findPatchID(patchName) == -1)
                    {
                        createPatch
                        (
                            mesh,
                            patchName,
                            (
                                (
                                    patchInfoDicts[meshi][sourcei].empty()
                                 && patchSources.size()
                                )
                              ? patchSources[0]
                              : -1
                            ),
                            patchDict,
                            word::null      // zoneGroupName
                        );
                        allAddedPatches[meshi].append(patchName);
                    }
                }
            }
        }

        Info<< endl;
    }

    // Patch indices have changed so clear the matchPatchIDs.
    matchPatchIDs.clear();


    // 2. Repatch faces
    // ~~~~~~~~~~~~~~~~

    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];
        auto& fzs = mesh.faceZones();

        Info<< "\n\nRepatching mesh " << mesh.name() << nl << endl;


        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const dictionary& dict = dicts[meshi];

        // Whether to synchronise points
        const bool pointSync(dict.getOrDefault<bool>("pointSync", false));


        // Read patch construct info from dictionary
        PtrList<dictionary> patchSources(dict.lookup("patches"));


        polyTopoChange meshMod(mesh);


        // Mark all repatched faces. This makes sure that the faces to repatch
        // do not overlap
        bitSet isRepatchedBoundary(mesh.nBoundaryFaces());

        forAll(patchSources, sourcei)
        {
            const dictionary& dict = patchSources[sourcei];
            const word& patchName = patchNames[meshi][sourcei];
            const word zoneGroupName("createPatch:" + patchName);
            const word sourceType(dict.get<word>("constructFrom"));
            const bool allowUndo = dict.getOrDefault("allowUndo", false);

            if (sourceType == "autoPatch")
            {
                forAll(interfaceMesh0, inti)
                {
                    const labelList& allMeshes0 = interfaceMesh0[inti];
                    const wordList& allNames0 = interfaceNames0[inti];
                    const auto& allFaces0 = interfaceFaces0[inti];

                    const labelList& allMeshes1 = interfaceMesh1[inti];
                    const wordList& allNames1 = interfaceNames1[inti];
                    const auto& allFaces1 = interfaceFaces1[inti];

                    if
                    (
                        interfaceSource0[inti] == sourcei
                     && allMeshes0[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh0. mesh1 is the remote mesh.

                        const label destPatchi =
                            patches.findPatchID(allNames0[sourcei], false);

                        //const auto& mesh1 =
                        //    meshes[allMeshes1[interfaceSource1[inti]]];
                        //Pout<< "Matched mesh:" << mesh.name()
                        //    << " to mesh:" << mesh1.name()
                        //    << " through:" << allNames0[sourcei] << endl;

                        changePatchID
                        (
                            mesh,
                            allFaces0[sourcei],
                            destPatchi,
                            -1,
                            true,   // optionally add to backup faceZone
                            isRepatchedBoundary,
                            meshMod
                        );
                    }
                    if
                    (
                        interfaceSource1[inti] == sourcei
                     && allMeshes1[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh1. mesh0 is the remote mesh.

                        const label destPatchi =
                            patches.findPatchID(allNames1[sourcei], false);

                        //const auto& mesh0 =
                        //    meshes[allMeshes0[interfaceSource0[inti]]];
                        //Pout<< "Matched mesh:" << mesh.name()
                        //    << " to mesh:" << mesh0.name()
                        //    << " through:" << allNames1[sourcei] << endl;

                        changePatchID
                        (
                            mesh,
                            allFaces1[sourcei],
                            destPatchi,
                            -1,
                            true,   // optionally add to backup faceZone
                            isRepatchedBoundary,
                            meshMod
                        );
                    }
                }
            }
            else if (undo)
            {
                // Undoing : only undoes patching. Does not restore
                // faceZones or faceSets.
                const auto patchSources = fzs.groupZoneIDs()[zoneGroupName];
                wordList sourcePatchNames(patchSources.size());
                forAll(patchSources, i)
                {
                    sourcePatchNames[i] = fzs[patchSources[i]].name();
                }

                Info<< "Splitting patch " << patchName
                    << " into " << flatOutput(sourcePatchNames) << endl;
                for (const label zonei : patchSources)
                {
                    const auto& fz = fzs[zonei];
                    const label patchi = patches.findPatchID(fz.name(), false);
                    changePatchID
                    (
                        mesh,
                        fz.addressing(),
                        patchi,
                        zonei,   // backup faceZone
                        false,   // remove from backup faceZone
                        isRepatchedBoundary,
                        meshMod
                    );
                }
            }
            else if (sourceType == "patches")
            {
                const label destPatchi = patches.findPatchID(patchName, false);

                // Could use the faceZoneGroup to find which patches to
                // combine to have more consistent ordering (if allowUndo)
                // const auto& patchSources =
                //     fzs.groupZoneIDs()[zoneGroupName];
                const labelList patchSources
                (
                    patches.patchSet
                    (
                        dict.get<wordRes>("patches")
                    ).sortedToc()
                );
                wordList sourcePatchNames(patchSources.size());
                forAll(patchSources, i)
                {
                    sourcePatchNames[i] = patches[patchSources[i]].name();
                }
                Info<< "Moving faces into " << patchName
                    << " from patches " << flatOutput(sourcePatchNames)
                    << endl;

                // Repatch faces of the patches.
                for (const label patchi : patchSources)
                {
                    const polyPatch& pp = patches[patchi];

                    // Add to optional backup faceZone
                    const label zonei =
                    (
                        allowUndo
                        ? fzs.findZoneID(pp.name())
                        : -1
                    );
                    changePatchID
                    (
                        mesh,
                        pp.start()+identity(pp.size()),
                        destPatchi,
                        zonei,  // backup faceZone
                        true,   // add to backup faceZone
                        isRepatchedBoundary,
                        meshMod
                    );
                }
            }
            else if (sourceType == "set")
            {
                const label destPatchi = patches.findPatchID(patchName, false);
                const word setName(dict.get<word>("set"));

                const faceSet set(mesh, setName);

                Info<< "Read " << returnReduce(set.size(), sumOp<label>())
                    << " faces from faceSet " << set.name() << endl;

                // Sort (since faceSet contains faces in arbitrary order)
                const labelList faceLabels(set.sortedToc());

                changePatchID
                (
                    allowUndo,
                    mesh,
                    faceLabels,
                    destPatchi,
                    isRepatchedBoundary,
                    meshMod
                );
            }
            else if (sourceType == "zone")
            {
                const label destPatchi = patches.findPatchID(patchName, false);
                const word zoneName(dict.get<word>("zone"));

                const auto& set = fzs[zoneName];

                Info<< "Read " << returnReduce(set.size(), sumOp<label>())
                    << " faces from faceZone " << set.name() << endl;

                labelList faceLabels(set.addressing());
                sort(faceLabels);

                changePatchID
                (
                    allowUndo,
                    mesh,
                    faceLabels,
                    destPatchi,
                    isRepatchedBoundary,
                    meshMod
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Invalid source type " << sourceType << endl
                    << "Valid source types are 'patches' 'set'"
                    << exit(FatalError);
            }
        }


        // Change mesh, use inflation to reforce calculation of transformation
        // tensors.
        Info<< "Doing topology modification to order faces." << nl << endl;
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);

        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
           // Force calculation of transformation tensors
            mesh.movePoints(pointField(mesh.points()));
        }

        // Update fields (only map within patch)
        mesh.updateMesh(map());


        // Update faceZones if undoing
        if (undo)
        {
            forAll(patchSources, sourcei)
            {
                const word& patchName = patchNames[meshi][sourcei];
                const word zoneGroupName("createPatch:" + patchName);

                // Get the patch ordering from the faceZoneGroup
                const auto patchSources = fzs.groupZoneIDs()[zoneGroupName];

                DynamicList<word> removeFaceZones;
                for (const label zonei : patchSources)
                {
                    const auto& fz = fzs[zonei];

                    // Remove combined group for undoing.
                    wordList& gNames = const_cast<wordList&>(fz.inGroups());
                    fzs.clearAddressing();
                    const label gIndex = gNames.find(zoneGroupName);
                    if (gIndex != -1)
                    {
                        std::swap(gNames[gIndex], gNames.last());
                        gNames.pop_back();
                    }
                    // Clear out any empty faceZones.
                    if (returnReduce(fz.size(), sumOp<label>()) == 0)
                    {
                        removeFaceZones.append(fz.name());
                    }
                }
                for (const word& fzName : removeFaceZones)
                {
                    const label zonei = fzs.findZoneID(fzName);
                    if (zonei != -1)
                    {
                        Info<< "Removing empty faceZone " << fzName << endl;
                        labelList oldToNew(identity(fzs.size()));
                        std::swap(oldToNew[zonei], oldToNew.last());
                        fzs.reorder(oldToNew);
                        fzs.clearAddressing();
                        fzs.setSize(fzs.size()-1);
                    }
                }
            }
        }

        // E.g. new patch 6 gets constructed (through constructFrom 'set') from
        // face 12 in patch 2, face 5 in patch 0, face 3 in patch 2 etc.
        // We now want to construct:
        //   - oldPatchIDs: (2 0)
        //   - patchFaceMap[0] = (5)
        //   - patchFaceMap[2] = (12 3)
        // and then use that to
        //   - construct patch field on patch 6 from patch 1 for faces (5)
        //     and from patch 0 for faces (12 3)
        {
            const auto& faceMap = map().faceMap();
            const auto& oldPatchStarts = map().oldPatchStarts();
            const auto& oldPatchSizes = map().oldPatchSizes();

            // Detect which patches completely originate from a subset of other
            // patches
            for (const auto& ppName : allAddedPatches[meshi])
            {
                const label patchi = patches.findIndex(ppName);
                const auto& pp = patches[patchi];

                labelListList patchFaceMap(oldPatchStarts.size());

                label maxOldPatchi = -1;
                forAll(pp, i)
                {
                    const label facei = pp.start()+i;
                    const label oldFacei = faceMap[facei];
                    const label oldPatchi =
                        findLower(oldPatchStarts, oldFacei+1);
                    maxOldPatchi = Foam::max(maxOldPatchi, oldPatchi);
                    // Make sure enough space. Nop for further invocations.
                    patchFaceMap[oldPatchi].setSize
                    (
                        oldPatchSizes[oldPatchi],
                        -1
                    );
                    const label localFacei = oldFacei-oldPatchStarts[oldPatchi];
                    patchFaceMap[oldPatchi][localFacei] = i;
                }
                reduce(maxOldPatchi, maxOp<label>());


                // Collect any non-zero sized patches
                DynamicList<label> oldPatchIDs(patches.size());
                for (label i = 0; i <= maxOldPatchi; i++)
                {
                    if (returnReduceOr(patchFaceMap[i].size()))
                    {
                        oldPatchIDs.append(i);
                    }
                }

                Info<< "Mapping fields on patch " << ppName
                    << " from patches " << flatOutput(oldPatchIDs) << endl;

                // Overwrite patch field with mapped one
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldVss[meshi],
                    vsFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldVvs[meshi],
                    vvFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldVsts[meshi],
                    vstFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldVsyms[meshi],
                    vsymtFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldVts[meshi],
                    vtFlds[meshi]
                );

                // Surface fields
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldSsBflds[meshi],
                    ssFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldSvBflds[meshi],
                    svFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldSstBflds[meshi],
                    sstFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldSsymtBflds[meshi],
                    ssymtFlds[meshi]
                );
                mapBFields
                (
                    patchi,
                    oldPatchIDs,
                    patchFaceMap,
                    oldStBflds[meshi],
                    stFlds[meshi]
                );
            }
        }


        // Update numbering pointing to meshi
        const auto& oldToNew = map().reverseFaceMap();
        forAll(interfaceMesh0, inti)
        {
            const labelList& allMeshes0 = interfaceMesh0[inti];
            forAll(allMeshes0, sourcei)
            {
                if (allMeshes0[sourcei] == meshi)
                {
                    inplaceRenumber(oldToNew, interfaceFaces0[inti][sourcei]);
                }
            }
        }
        forAll(interfaceMesh1, inti)
        {
            const labelList& allMeshes1 = interfaceMesh1[inti];
            forAll(allMeshes1, sourcei)
            {
                if (allMeshes1[sourcei] == meshi)
                {
                    inplaceRenumber(oldToNew, interfaceFaces1[inti][sourcei]);
                }
            }
        }


        if (writeObj)
        {
            dumpCyclicMatch("coupled_", mesh);
        }

        // Synchronise points.
        if (!pointSync)
        {
            Info<< "Not synchronising points." << nl << endl;
        }
        else
        {
            Info<< "Synchronising points." << nl << endl;

            // This is a bit tricky. Both normal and position might be out and
            // current separation also includes the normal
            // ( separation_ = (nf&(Cr - Cf))*nf ).

            // For cyclic patches:
            // - for separated ones use user specified offset vector

            forAll(mesh.boundaryMesh(), patchi)
            {
                const polyPatch& pp = mesh.boundaryMesh()[patchi];

                if (pp.size() && isA<coupledPolyPatch>(pp))
                {
                    const coupledPolyPatch& cpp =
                        refCast<const coupledPolyPatch>(pp);

                    if (cpp.separated())
                    {
                        Info<< "On coupled patch " << pp.name()
                            << " separation[0] was "
                            << cpp.separation()[0] << endl;

                        if (isA<cyclicPolyPatch>(pp) && pp.size())
                        {
                            const auto& cycpp =
                                refCast<const cyclicPolyPatch>(pp);

                            if
                            (
                                cycpp.transform()
                             == cyclicPolyPatch::TRANSLATIONAL
                            )
                            {
                                // Force to wanted separation
                                Info<< "On cyclic translation patch "
                                    << pp.name()
                                    << " forcing uniform separation of "
                                    << cycpp.separationVector() << endl;
                                const_cast<vectorField&>(cpp.separation()) =
                                    pointField(1, cycpp.separationVector());
                            }
                            else
                            {
                                const auto& nbr = cycpp.neighbPatch();
                                const_cast<vectorField&>(cpp.separation()) =
                                    pointField
                                    (
                                        1,
                                        nbr[0].centre(mesh.points())
                                      - cycpp[0].centre(mesh.points())
                                    );
                            }
                        }
                        Info<< "On coupled patch " << pp.name()
                            << " forcing uniform separation of "
                            << cpp.separation() << endl;
                    }
                    else if (!cpp.parallel())
                    {
                        Info<< "On coupled patch " << pp.name()
                            << " forcing uniform rotation of "
                            << cpp.forwardT()[0] << endl;

                        const_cast<tensorField&>
                        (
                            cpp.forwardT()
                        ).setSize(1);
                        const_cast<tensorField&>
                        (
                            cpp.reverseT()
                        ).setSize(1);

                        Info<< "On coupled patch " << pp.name()
                            << " forcing uniform rotation of "
                            << cpp.forwardT() << endl;
                    }
                }
            }

            Info<< "Synchronising points." << endl;

            pointField newPoints(mesh.points());

            syncPoints
            (
                mesh,
                newPoints,
                minMagSqrEqOp<vector>(),
                point(GREAT, GREAT, GREAT)
            );

            scalarField diff(mag(newPoints-mesh.points()));
            Info<< "Points changed by average:" << gAverage(diff)
                << " max:" << gMax(diff) << nl << endl;

            mesh.movePoints(newPoints);
        }


        // 3. Remove zeros-sized patches
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Info<< "Removing patches with no faces in them." << nl << endl;
        const wordList oldPatchNames(mesh.boundaryMesh().names());
        const wordList oldPatchTypes(mesh.boundaryMesh().types());
        fvMeshTools::removeEmptyPatches(mesh, allAddedPatches[meshi], true);
        forAll(oldPatchNames, patchi)
        {
            const word& pName = oldPatchNames[patchi];
            if (mesh.boundaryMesh().findPatchID(pName) == -1)
            {
                Info<< "Removed zero-sized patch " << pName
                    << " type " << oldPatchTypes[patchi]
                    << " at position " << patchi << endl;
            }
        }


        if (writeObj)
        {
            dumpCyclicMatch("final_", mesh);
        }
    }

    if (!overwrite)
    {
        ++runTime;
    }

    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];
        mesh.setInstance(overwrite ? oldInstances[meshi] : runTime.timeName());
    }

    // More precision (for points data)
    IOstream::minPrecision(10);

    // Write resulting mesh
    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];

        // Override bcs with explicitly provided info. Done late so there
        // are already patch faces
        forAll(patchInfoDicts[meshi], sourcei)
        {
            const dictionary& patchDict = patchInfoDicts[meshi][sourcei];
            const word& patchName = patchNames[meshi][sourcei];
            const label patchID = mesh.boundary().findPatchID(patchName);
            if (patchID != -1 && patchDict.found("patchFields"))
            {
                const dictionary& pfd = patchDict.subDict("patchFields");
                fvMeshTools::setPatchFields(mesh, patchID, pfd);
            }
        }

        const word& fieldInst = runTime.timeName();
        auto setInstance = [&fieldInst](auto& flds)
        {
            for (auto& fld : flds) { fld.instance() = fieldInst; }
        };
        // Volume fields.
        setInstance(vsFlds[meshi]);
        setInstance(vvFlds[meshi]);
        setInstance(vstFlds[meshi]);
        setInstance(vsymtFlds[meshi]);
        setInstance(vtFlds[meshi]);

        // Surface fields.
        setInstance(ssFlds[meshi]);
        setInstance(svFlds[meshi]);
        setInstance(sstFlds[meshi]);
        setInstance(ssymtFlds[meshi]);
        setInstance(stFlds[meshi]);

        Info<< "\n\nWriting repatched mesh " << mesh.name()
            << " to " << runTime.timeName() << nl << endl;
        mesh.clearOut();    // remove meshPhi
        mesh.write();
        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
