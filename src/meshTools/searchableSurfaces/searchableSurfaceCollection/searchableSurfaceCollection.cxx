/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "searchableSurfaceCollection.H"
#include "CompactListList.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceCollection, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableSurfaceCollection,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableSurfaceCollection,
        dict,
        collection
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceCollection::findNearest_impl
(
    const pointField& samples,
    scalarField& minDistSqr,
    List<pointIndexHit>& nearestInfo,
    labelList* nearestSurface
) const
{
    // Initialise
    nearestInfo.resize_fill(samples.size(), pointIndexHit());

    if (nearestSurface)
    {
        nearestSurface->resize_fill(samples.size(), -1);
    }

    List<pointIndexHit> hitInfo(samples.size());

    const scalarField localMinDistSqr(samples.size(), GREAT);

    forAll(subGeom_, surfI)
    {
        subGeom_[surfI].findNearest
        (
            cmptDivide  // Transform then divide
            (
                transform_[surfI].localPosition(samples),
                scale_[surfI]
            ),
            localMinDistSqr,
            hitInfo
        );

        forAll(hitInfo, pointi)
        {
            if (hitInfo[pointi].hit())
            {
                // Rework back into global coordinate sys. Multiply then
                // transform
                point globalPt = transform_[surfI].globalPosition
                (
                    cmptMultiply
                    (
                        hitInfo[pointi].point(),
                        scale_[surfI]
                    )
                );

                scalar distSqr = globalPt.distSqr(samples[pointi]);

                if (distSqr < minDistSqr[pointi])
                {
                    minDistSqr[pointi] = distSqr;
                    nearestInfo[pointi].hitPoint(globalPt);
                    nearestInfo[pointi].setIndex
                    (
                        hitInfo[pointi].index()
                      + indexOffset_.localStart(surfI)
                    );

                    if (nearestSurface)
                    {
                        (*nearestSurface)[pointi] = surfI;
                    }
                }
            }
        }
    }
}


// Sort hits into per-surface bins. Misses are rejected.
// Maintains map back to position
void Foam::searchableSurfaceCollection::sortHits
(
    const UList<pointIndexHit>& info,
    CompactListList<pointIndexHit>& surfInfo,
    CompactListList<label>& infoMap
) const
{
    // Count hits per surface.
    labelList nHits(subGeom_.size(), Zero);

    forAll(info, pointi)
    {
        if (info[pointi].hit())
        {
            label index = info[pointi].index();
            label surfi = indexOffset_.findRow(index);
            nHits[surfi]++;
        }
    }

    // Per surface the hit
    surfInfo.resize_nocopy(nHits);

    // Per surface the original position
    infoMap.resize_nocopy(nHits);

    nHits = 0;

    forAll(info, pointi)
    {
        if (info[pointi].hit())
        {
            const label index = info[pointi].index();
            const label surfi = indexOffset_.findRow(index);

            // current (local row) output location
            const label curr = nHits[surfi]++;

            // Store for correct surface and adapt indices back to local
            // ones

            auto hits = surfInfo[surfi];
            auto map = infoMap[surfi];

            hits[curr].hitPoint
            (
                info[pointi].point(),
                indexOffset_.toLocal(surfi, index)
            );
            map[curr] = pointi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceCollection::searchableSurfaceCollection
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    instance_(dict.size()),
    scale_(dict.size()),
    transform_(dict.size()),
    subGeom_(dict.size()),
    mergeSubRegions_(dict.get<bool>("mergeSubRegions"))
{
    Info<< "SearchableCollection : " << name() << endl;

    DynamicList<label> surfSizes(dict.size());
    label surfI = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            instance_[surfI] = dEntry.keyword();

            const dictionary& sDict = dEntry.dict();

            sDict.readEntry("scale", scale_[surfI]);

            // Mandatory 'transform' sub-dictionary
            transform_.set
            (
                surfI,
                new coordSystem::cartesian(sDict, "transform")
            );

            const word subGeomName(sDict.get<word>("surface"));
            //Pout<< "Trying to find " << subGeomName << endl;

            searchableSurface& s =
                io.db().lookupObjectRef<searchableSurface>(subGeomName);

            // I don't know yet how to handle the globalSize combined with
            // regionOffset. Would cause non-consecutive indices locally
            // if all indices offset by globalSize() of the local region...
            if (s.size() != s.globalSize())
            {
                FatalErrorInFunction
                    << "Cannot use a distributed surface in a collection."
                    << exit(FatalError);
            }

            subGeom_.set(surfI, &s);

            surfSizes.push_back(s.size());

            Info<< "    instance : " << instance_[surfI] << endl;
            Info<< "    surface  : " << s.name() << endl;
            Info<< "    scale    : " << scale_[surfI] << endl;
            Info<< "    transform: " << transform_[surfI] << endl;

            surfI++;
        }
    }

    instance_.setSize(surfI);
    scale_.setSize(surfI);
    transform_.setSize(surfI);
    subGeom_.setSize(surfI);

    indexOffset_.reset(surfSizes);

    // Bounds is the overall bounds - prepare for min/max ops
    bounds().reset();

    forAll(subGeom_, surfI)
    {
        const boundBox& surfBb = subGeom_[surfI].bounds();

        // Transform back to global coordinate sys.
        const point surfBbMin = transform_[surfI].globalPosition
        (
            cmptMultiply
            (
                surfBb.min(),
                scale_[surfI]
            )
        );
        const point surfBbMax = transform_[surfI].globalPosition
        (
            cmptMultiply
            (
                surfBb.max(),
                scale_[surfI]
            )
        );

        bounds().min() = min(bounds().min(), surfBbMin);
        bounds().max() = max(bounds().max(), surfBbMax);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceCollection::~searchableSurfaceCollection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableSurfaceCollection::regions() const
{
    if (regions_.empty() && !subGeom_.empty())
    {
        regionOffset_.resize_nocopy(subGeom_.size());

        DynamicList<word> allRegions;

        if (mergeSubRegions_)
        {
            // Single name per surface, regardless how many regions any
            // sub-surface has.

            // The region offset (when merged) is an identity map
            Foam::identity(regionOffset_);

            forAll(subGeom_, i)
            {
                allRegions.push_back(instance_[i] + "_" + Foam::name(i));
            }
        }
        else
        {
            // Retain all sub-surface region names
            forAll(subGeom_, i)
            {
                regionOffset_[i] = allRegions.size();

                for (const word& regionName : subGeom_[i].regions())
                {
                    allRegions.push_back(instance_[i] + "_" + regionName);
                }
            }
        }
        regions_ = std::move(allRegions);
    }
    return regions_;
}


Foam::label Foam::searchableSurfaceCollection::size() const
{
    // The total number of triangles etc in the collection
    return indexOffset_.totalSize();
}


Foam::label
Foam::searchableSurfaceCollection::nPoints() const
{
    label total = 0;
    for (const auto& geo : subGeom_)
    {
        total += geo.nPoints();
    }
    return total;
}


Foam::tmp<Foam::pointField>
Foam::searchableSurfaceCollection::coordinates() const
{
    auto tctrs = tmp<pointField>::New(size());
    auto& ctrs = tctrs.ref();

    // Append individual coordinates
    label coordI = 0;

    forAll(subGeom_, surfI)
    {
        const pointField subCoords(subGeom_[surfI].coordinates());

        forAll(subCoords, i)
        {
            ctrs[coordI++] = transform_[surfI].globalPosition
            (
                cmptMultiply
                (
                    subCoords[i],
                    scale_[surfI]
                )
            );
        }
    }

    return tctrs;
}


void Foam::searchableSurfaceCollection::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.resize_nocopy(size());
    radiusSqr.resize_nocopy(centres.size());

    // Append individual coordinates
    label coordI = 0;

    forAll(subGeom_, surfI)
    {
        const scalar maxScaleSqr = Foam::sqr(cmptMax(scale_[surfI]));

        pointField subCentres;
        scalarField subRadiusSqr;
        subGeom_[surfI].boundingSpheres(subCentres, subRadiusSqr);

        forAll(subCentres, i)
        {
            centres[coordI] = transform_[surfI].globalPosition
            (
                cmptMultiply
                (
                    subCentres[i],
                    scale_[surfI]
                )
            );
            radiusSqr[coordI] = maxScaleSqr*subRadiusSqr[i];
            coordI++;
        }
    }
}


Foam::tmp<Foam::pointField>
Foam::searchableSurfaceCollection::points() const
{
    // Get overall size
    label nPoints = 0;

    for (const auto& geo : subGeom_)
    {
        nPoints += geo.nPoints();
    }

    auto tpts = tmp<pointField>::New(nPoints);
    auto& pts = tpts.ref();

    // Append individual coordinates
    nPoints = 0;

    forAll(subGeom_, surfI)
    {
        const pointField subCoords(subGeom_[surfI].points());

        forAll(subCoords, i)
        {
            pts[nPoints++] = transform_[surfI].globalPosition
            (
                cmptMultiply
                (
                    subCoords[i],
                    scale_[surfI]
                )
            );
        }
    }

    return tpts;
}


void Foam::searchableSurfaceCollection::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& nearestInfo
) const
{
    // How to scale distance?
    scalarField minDistSqr(nearestDistSqr);
    findNearest_impl(samples, minDistSqr, nearestInfo);
}


void Foam::searchableSurfaceCollection::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize_fill(start.size(), pointIndexHit());

    // Current nearest (to start) intersection
    pointField nearest(end);

    List<pointIndexHit> hitInfo(start.size());

    forAll(subGeom_, surfI)
    {
        // Starting point
        tmp<pointField> e0 = cmptDivide
        (
            transform_[surfI].localPosition
            (
                start
            ),
            scale_[surfI]
        );

        // Current best end point
        tmp<pointField> e1 = cmptDivide
        (
            transform_[surfI].localPosition
            (
                nearest
            ),
            scale_[surfI]
        );

        subGeom_[surfI].findLine(e0, e1, hitInfo);

        forAll(hitInfo, pointi)
        {
            if (hitInfo[pointi].hit())
            {
                // Transform back to global coordinate sys.
                nearest[pointi] = transform_[surfI].globalPosition
                (
                    cmptMultiply
                    (
                        hitInfo[pointi].point(),
                        scale_[surfI]
                    )
                );
                info[pointi] = hitInfo[pointi];
                info[pointi].point() = nearest[pointi];
                info[pointi].setIndex
                (
                    hitInfo[pointi].index()
                  + indexOffset_.localStart(surfI)
                );
            }
        }
    }


    // Debug check
    if constexpr (false)
    {
        forAll(info, pointi)
        {
            if (info[pointi].hit())
            {
                vector n(end[pointi] - start[pointi]);
                scalar magN = mag(n);

                if (magN > SMALL)
                {
                    n /= mag(n);

                    scalar s = ((info[pointi].point()-start[pointi])&n);

                    if (s < 0 || s > 1)
                    {
                        FatalErrorInFunction
                            << "point:" << info[pointi]
                            << " s:" << s
                            << " outside vector "
                            << " start:" << start[pointi]
                            << " end:" << end[pointi]
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


void Foam::searchableSurfaceCollection::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    // To be done ...
    findLine(start, end, info);
}


void Foam::searchableSurfaceCollection::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    // To be done. Assume for now only one intersection.
    List<pointIndexHit> nearestInfo;
    findLine(start, end, nearestInfo);

    info.resize_nocopy(start.size());

    forAll(info, pointi)
    {
        info[pointi].clear();

        if (nearestInfo[pointi].hit())
        {
            info[pointi].resize_fill(1, nearestInfo[pointi]);
        }
    }
}


void Foam::searchableSurfaceCollection::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    if (subGeom_.size() == 0)
    {
        region.clear();
    }
    else if (subGeom_.size() == 1)
    {
        if (mergeSubRegions_)
        {
            // The region is identical to the surface index (0)
            region.resize_fill(info.size(), 0);
        }
        else
        {
            subGeom_[0].getRegion(info, region);
        }
    }
    else
    {
        // Multiple surfaces. Sort by surface.

        if (!mergeSubRegions_)
        {
            // The regionOffset_ must be populated (when not merged)
            (void) this->regions();
        }

        // Per surface the hit
        CompactListList<pointIndexHit> surfInfo;
        // Per surface the original index
        CompactListList<label> infoMap;
        sortHits(info, surfInfo, infoMap);

        region.resize_fill(info.size(), -1);

        // Do region tests

        if (mergeSubRegions_)
        {
            // Actually no need for surfInfo.
            // The region is identical to the surface index.
            forAll(infoMap, surfi)
            {
                // The region offset is fixed:
                const label regOffset = surfi;

                UIndirectList<label>(region, infoMap[surfi]) = regOffset;
            }
        }
        else
        {
            // The regionOffset_ must be populated (when not merged)
            (void) this->regions();

            forAll(infoMap, surfI)
            {
                // The region offset is variable:
                const label regOffset = regionOffset_[surfI];

                labelList surfRegion;
                subGeom_[surfI].getRegion
                (
                    // The API needs an update!
                    static_cast<const List<pointIndexHit>&>(surfInfo[surfI]),
                    surfRegion
                );

                UIndirectList<label> mapped(region, infoMap[surfI]);

                std::transform
                (
                    surfRegion.begin(), surfRegion.end(),
                    mapped.begin(),
                    [regOffset](label r){ return (r + regOffset); }
                );
            }
        }
    }
}


void Foam::searchableSurfaceCollection::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (subGeom_.size() == 0)
    {}
    else if (subGeom_.size() == 1)
    {
        subGeom_[0].getNormal(info, normal);
    }
    else
    {
        // Multiple surfaces. Sort by surface.

        // Per surface the hit
        CompactListList<pointIndexHit> surfInfo;
        // Per surface the original index
        CompactListList<label> infoMap;
        sortHits(info, surfInfo, infoMap);

        normal.resize_nocopy(info.size());

        // Do region tests
        forAll(surfInfo, surfI)
        {
            vectorField surfNormal;
            subGeom_[surfI].getNormal
            (
                // The API needs an update!
                static_cast<const List<pointIndexHit>&>(surfInfo[surfI]),
                surfNormal
            );

            // Transform back to global coordinate sys.
            surfNormal = transform_[surfI].globalVector(surfNormal);

            UIndirectList<vector>(normal, infoMap[surfI]) = surfNormal;
        }
    }
}


void Foam::searchableSurfaceCollection::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorInFunction
        << "Volume type not supported for collection."
        << exit(FatalError);
}


void Foam::searchableSurfaceCollection::distribute
(
    const List<treeBoundBox>& bbs,
    const bool keepNonLocal,
    autoPtr<mapDistribute>& faceMap,
    autoPtr<mapDistribute>& pointMap
)
{
    forAll(subGeom_, surfI)
    {
        // Note:Transform the bounding boxes? Something like
        // pointField bbPoints =
        // cmptDivide
        // (
        //     transform_[surfI].localPosition(bbs[i].points()),
        //     scale_[surfI]
        // );
        // treeBoundBox newBb(bbPoints);

        // Note: what to do with faceMap, pointMap from multiple surfaces?
        subGeom_[surfI].distribute
        (
            bbs,
            keepNonLocal,
            faceMap,
            pointMap
        );
    }
}


void Foam::searchableSurfaceCollection::setField(const labelList& values)
{
    forAll(subGeom_, surfI)
    {
        subGeom_[surfI].setField
        (
            // The API needs an update!
            static_cast<const labelList&>
            (
                SubList<label>
                (
                    values,
                    subGeom_[surfI].size(),
                    indexOffset_.localStart(surfI)
                )
            )
        );
    }
}


void Foam::searchableSurfaceCollection::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    if (subGeom_.size() == 0)
    {
        return;
    }
    if (subGeom_.size() == 1)
    {
        subGeom_[0].getField(info, values);
    }
    else
    {
        // Multiple surfaces. Sort by surface.

        // Per surface the hit
        CompactListList<pointIndexHit> surfInfo;
        // Per surface the original index
        CompactListList<label> infoMap;
        sortHits(info, surfInfo, infoMap);

        // Do surface tests
        forAll(surfInfo, surfI)
        {
            labelList surfValues;
            subGeom_[surfI].getField
            (
                // The API needs an update!
                static_cast<const List<pointIndexHit>&>(surfInfo[surfI]),
                surfValues
            );

            if (surfValues.size())
            {
                // Size values only when we have a surface that supports it.
                values.resize(info.size(), -1);

                // We do not trust that surfValues will always match in
                // length (smaller or larger), so reslice ourselves

                const label subSize
                (
                    std::min(infoMap.localSize(surfI), surfValues.size())
                );

                const auto subInfoMap =
                (
                    infoMap.values().slice(infoMap.localStart(surfI), subSize)
                );

                UIndirectList<label>(values, subInfoMap) =
                    surfValues.slice(0, subSize);
            }
        }
    }
}


Foam::label Foam::searchableSurfaceCollection::measure
(
    List<vector>* centres,
    List<scalar>* radiusSqr,
    List<label>* regions,
    List<vector>* normals
) const
{
    if (subGeom_.size() == 0)
    {
        // Same as searchableSurface::measure() but be explicit about it
        if (centres) centres->clear();
        if (radiusSqr) radiusSqr->clear();
        if (regions) regions->clear();
        if (normals) normals->clear();
        return 0;
    }
    else if (subGeom_.size() == 1)
    {
        // Note, when merged: every entry belongs to the single collection region,
        // so the sub-surface regions are not needed

        const label output_size = subGeom_[0].size();

        const label subCount = subGeom_[0].measure
        (
            centres,
            radiusSqr,
            (mergeSubRegions_ ? nullptr : regions),
            normals
        );

        bool usedParams = (subCount > 0);

        if (mergeSubRegions_ && regions)
        {
            // Single merged surface, so region == 0
            usedParams = true;
            regions->resize_fill(output_size, 0);
        }

        const vector& scale = scale_[0];
        const auto& transf = transform_[0];

        if (centres)
        {
            // inplace transform the centres
            std::transform
            (
                centres->begin(), centres->end(),
                centres->begin(),

                [&scale,&transf](const vector& cc)
                {
                    return transf.globalPosition(cmptMultiply(cc, scale));
                }
            );
        }

        if (radiusSqr)
        {
            const scalar maxScaleSqr = Foam::sqr(cmptMax(scale));

            // inplace rescale the radius squared
            std::transform
            (
                radiusSqr->begin(), radiusSqr->end(),
                radiusSqr->begin(),

                [maxScaleSqr](scalar d2) { return maxScaleSqr*d2; }
            );
        }

        if (normals)
        {
            // inplace transform the normals
            std::transform
            (
                normals->begin(), normals->end(),
                normals->begin(),

                [&scale,&transf](const vector& n)
                {
                    return normalised
                    (
                        transf.globalVector(cmptDivide(n, scale))
                    );
                }
            );
        }

        return (usedParams ? output_size : 0);
    }

    const label output_size = size();

    bool usedParams = false;

    std::unique_ptr<List<point>> subCentres_ptr;
    if (centres)
    {
        usedParams = true;
        centres->resize_nocopy(output_size);
        subCentres_ptr = std::make_unique<List<point>>();
    }

    std::unique_ptr<List<scalar>> subRadiusSqr_ptr;
    if (radiusSqr)
    {
        usedParams = true;
        radiusSqr->resize_nocopy(output_size);
        subRadiusSqr_ptr = std::make_unique<List<scalar>>();
    }

    std::unique_ptr<List<label>> subRegions_ptr;
    if (regions)
    {
        usedParams = true;
        regions->resize_nocopy(output_size);

        if (!mergeSubRegions_)
        {
            // The regionOffset_ must be populated (when not merged)
            (void) this->regions();

            subRegions_ptr = std::make_unique<List<label>>();
        }
    }

    std::unique_ptr<List<vector>> subNormals_ptr;
    if (normals)
    {
        usedParams = true;
        normals->resize_nocopy(output_size);
        subNormals_ptr = std::make_unique<List<vector>>();
    }

    if (!usedParams)
    {
        return 0;
    }

    // After here we have one or more parameters to measure()
    // and fit back into the overall list

    forAll(subGeom_, surfi)
    {
        // Offset into the combined output, local size
        const label localStart = indexOffset_.localStart(surfi);
        const label localSize = indexOffset_.localSize(surfi);

        const vector& scale = scale_[surfi];
        const scalar maxScaleSqr = Foam::sqr(cmptMax(scale));
        const auto& transf = transform_[surfi];

        // Region offset can be fixed (mergeSubRegions) or variable
        const label regOffset
        (
            (regions && !mergeSubRegions_)
          ? regionOffset_[surfi]   // variable offset
          : surfi                  // fixed offset
        );

        const label subCount = subGeom_[surfi].measure
        (
            subCentres_ptr.get(),
            subRadiusSqr_ptr.get(),
            subRegions_ptr.get(),
            subNormals_ptr.get()
        );

        if (FOAM_UNLIKELY(subCount != localSize))
        {
            #ifdef FULLDEBUG
            {
                FatalErrorInFunction
                    << "measure of sub-surface:" << surfi
                    << " returned " << subCount << ", expected " << localSize
                    << endl;
            }
            #endif

            if (debug)
            {
                WarningInFunction
                    << "measure of sub-surface:" << surfi
                    << " returned " << subCount << ", expected " << localSize
                    << endl;
            }
        }

        if (centres && subCentres_ptr)
        {
            const auto& subCentres = *subCentres_ptr;

            // Transform back to global coordinate sys.
            std::transform
            (
                subCentres.begin(), subCentres.end(),
                (centres->begin() + localStart),

                [&scale,&transf](const vector& cc)
                {
                    return transf.globalPosition(cmptMultiply(cc, scale));
                }
            );
        }

        if (radiusSqr && subRadiusSqr_ptr)
        {
            const auto& subRadiusSqr = *subRadiusSqr_ptr;

            std::transform
            (
                subRadiusSqr.begin(), subRadiusSqr.end(),
                (radiusSqr->begin() + localStart),

                [maxScaleSqr](scalar d2) { return maxScaleSqr*d2; }
            );

        }

        if (regions && subRegions_ptr && !mergeSubRegions_)
        {
            const auto& subRegions = *subRegions_ptr;

            std::transform
            (
                subRegions.begin(), subRegions.end(),
                (regions->begin() + localStart),

                [regOffset](label idx) { return (regOffset + idx); }
            );
        }
        else if (regions)
        {
            // Merged (per surface region id)
            std::fill_n(regions->begin() + localStart, localSize, regOffset);
        }

        if (normals && subNormals_ptr)
        {
            const auto& subNormals = *subNormals_ptr;

            // Transform back to global coordinate sys.
            std::transform
            (
                subNormals.begin(), subNormals.end(),
                (normals->begin() + localStart),

                [&scale,&transf](const vector& n)
                {
                    return normalised
                    (
                        transf.globalVector(cmptDivide(n, scale))
                    );
                }
            );
        }
    }

    return (usedParams ? output_size : 0);
}


// ************************************************************************* //
