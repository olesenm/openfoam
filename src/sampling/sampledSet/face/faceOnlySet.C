/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "faceOnlySet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceOnlySet, 0);
    addToRunTimeSelectionTable(sampledSet, faceOnlySet, word);
}

const Foam::scalar Foam::faceOnlySet::tol = ROOTSMALL;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::faceOnlySet::trackToBoundary
(
    passiveParticleCloud& particleCloud,
    passiveParticle& singleParticle,
    const scalar smallDist,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<scalar>& samplingCurveDist
) const
{
    const vector offset = (end_ - start_);

    particle::trackingData td(particleCloud);

    point trackPt = singleParticle.position();

    while (true)
    {
        point oldPoint = trackPt;

        singleParticle.trackToAndHitFace(end_ - start_, 0, particleCloud, td);

        trackPt = singleParticle.position();

        if (singleParticle.face() != -1 && mag(oldPoint - trackPt) > smallDist)
        {
            // Reached face. Sample.
            samplingPts.append(trackPt);
            samplingCells.append(singleParticle.cell());
            samplingFaces.append(singleParticle.face());
            samplingCurveDist.append(mag(trackPt - start_));
        }

        if (-smallDist < ((trackPt - end_) & offset))
        {
            // Projected onto sampling vector
            // - done when we are near or past the end of the sampling vector
            return false;
        }
        else if (singleParticle.onBoundaryFace())
        {
            // Boundary reached
            return true;
        }
    }
}


void Foam::faceOnlySet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // Distance vector between sampling points
    if (mag(end_ - start_) < SMALL)
    {
        FatalErrorInFunction
            << "Incorrect sample specification :"
            << " start equals end point." << endl
            << "  start:" << start_
            << "  end:" << end_
            << exit(FatalError);
    }

    const vector offset = (end_ - start_);
    const vector normOffset = offset/mag(offset);
    const vector smallVec = tol*offset;
    const scalar smallDist = mag(smallVec);

    // Force calculation of cloud addressing on all processors
    const bool oldMoving = const_cast<polyMesh&>(mesh()).moving(false);
    passiveParticleCloud particleCloud(mesh());

    // Get all boundary intersections (includes processor patches for entry)
    List<pointIndexHit> bHits = searchEngine().intersections
    (
        start_ - smallVec,
        end_ + smallVec
    );

    point bPoint(GREAT, GREAT, GREAT);
    label bFacei = -1;

    if (bHits.size())
    {
        bPoint = bHits[0].hitPoint();
        bFacei = bHits[0].index();
    }

    // Get first tracking point. Use bPoint, bFacei if provided.
    point trackPt;
    label trackCelli = -1;
    label trackFacei = -1;

    // Pout<< "before getTrackingPoint : bPoint:" << bPoint
    //     << " bFacei:" << bFacei << endl;

    getTrackingPoint
    (
        start_,
        bPoint,
        bFacei,
        smallDist,
        trackPt,
        trackCelli,
        trackFacei
    );

    // Initial-start logic failed to find a containing cell; refine search
    if (trackCelli == -1)
    {
        // The line may pass through this processor subdomain without the
        // boundary intersection tree reporting a hit, e.g. neighbouring
        // processor patch ordering. Search along the segment.
        constexpr int nDiv = 100;
        for (int i = 0; i <= nDiv; ++i)
        {
            const point p(start_ + (scalar(i)/scalar(nDiv))*(end_ - start_));

            const label celli = mesh().findCell(p, searchEngine().decompMode());

            if (celli != -1)
            {
                trackPt = p;
                trackCelli = celli;
                trackFacei = -1;
                break;
            }
        }
    }

    // Pout<< "after getTrackingPoint : "
    //     << " trackPt:" << trackPt
    //     << " trackCelli:" << trackCelli
    //     << " trackFacei:" << trackFacei
    //     << endl;

    if (trackCelli == -1)
    {
        // Line start_ - end_ does not intersect domain at all.
        // (or is along edge)
        // Set points and cell/face labels to empty lists
        //Info<< "calcSamples : Both start_ and end_ outside domain"
        //    << endl;

        const_cast<polyMesh&>(mesh()).moving(oldMoving);
        return;
    }

    if (trackFacei == -1)
    {
        // No boundary face. Check for nearish internal face
        trackFacei = findNearFace(trackCelli, trackPt, smallDist);
    }

    // Pout<< "calcSamples : got first point to track from :"
    //     << "  trackPt:" << trackPt
    //     << "  trackCell:" << trackCelli
    //     << "  trackFace:" << trackFacei
    //     << endl;

    //
    // Track until hit end of all boundary intersections
    //

    // current segment number
    label segmentI = 0;

    // starting index of current segment in samplePts
    label startSegmentI = 0;

    // In parallel or cyclic cases, coupled patches are not physical exits from
    // overall domain, so advance to next meaningful non-coupled boundary hit
    label bHitI = 0;
    while (bHitI < bHits.size())
    {
        const label patchi =
            mesh().boundaryMesh().whichPatch(bHits[bHitI].index());

        if (mesh().boundaryMesh()[patchi].coupled())
        {
            ++bHitI;
            continue;
        }

        const scalar dist = (bHits[bHitI].hitPoint() - trackPt) & normOffset;

        if (dist > smallDist)
        {
            break;
        }

        ++bHitI;
    }

    while (true)
    {
        if (trackFacei != -1)
        {
            const label patchi =
                mesh().boundaryMesh().whichPatch(trackFacei);

            if (mesh().boundaryMesh()[patchi].coupled())
            {
                // Entering through processor/cyclic patch. Push into the
                // domain without sampling the coupled face.
                const label entryFacei = trackFacei;
                trackPt = pushIn(trackPt, entryFacei);
                trackCelli = mesh().findCell
                (
                    trackPt,
                    searchEngine().decompMode()
                );

                if (trackCelli == -1)
                {
                    trackCelli = getBoundaryCell(entryFacei);
                }

                trackFacei = -1;
            }
            else
            {
                // Pout<< "trackPt:" << trackPt << " on face so use." << endl;
                samplingPts.append(trackPt);
                samplingCells.append(trackCelli);
                samplingFaces.append(trackFacei);
                samplingCurveDist.append(mag(trackPt - start_));
            }
        }

        // Initialize tracking starting from trackPt
        passiveParticle singleParticle
        (
            mesh(),
            trackPt,
            trackCelli
        );

        bool reachedBoundary = trackToBoundary
        (
            particleCloud,
            singleParticle,
            smallDist,
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingCurveDist
        );

        // Fill sampleSegments
        for (label i = samplingPts.size() - 1; i >= startSegmentI; --i)
        {
            samplingSegments.append(segmentI);
        }

        if (!reachedBoundary)
        {
            // Pout<< "calcSamples : Reached end of samples: "
            //     << "  samplePt now:" << singleParticle.position()
            //     << endl;
            break;
        }

        if (singleParticle.onBoundaryFace())
        {
            const polyPatch& pp = mesh().boundaryMesh()[singleParticle.patch()];

            // Processor can own multiple disjoint pieces of the sampling line;
            // try to continue with a later local piece instead of stopping
            // permanently at the first coupled interface
            if (pp.coupled() && samplingPts.size() > startSegmentI)
            {
                const vector nf(mesh().faceAreas()[singleParticle.face()]);

                // Only terminate local tracking when leaving the subdomain
                if ((nf & normOffset) > 0)
                {
                    bool foundAhead = false;
                    const point pos(singleParticle.position());

                    for (int i = 1; i <= 100; ++i)
                    {
                        const point p
                        (
                            pos + (scalar(i)/scalar(100))*(end_ - pos)
                        );

                        const label celli =
                            mesh().findCell(p, searchEngine().decompMode());

                        if (celli != -1)
                        {
                            trackPt = p;
                            trackCelli = celli;
                            trackFacei = -1;
                            foundAhead = true;
                            break;
                        }
                    }

                    if (!foundAhead)
                    {
                        break;
                    }

                    ++segmentI;
                    startSegmentI = samplingPts.size();
                    continue;
                }
            }
        }

        bool foundValidB = false;

        while (bHitI < bHits.size())
        {
            const label patchi =
                mesh().boundaryMesh().whichPatch(bHits[bHitI].index());

            if (mesh().boundaryMesh()[patchi].coupled())
            {
                ++bHitI;
                continue;
            }

            scalar dist =
                (bHits[bHitI].hitPoint() - singleParticle.position())
              & normOffset;

            // Pout<< "Finding next boundary : "
            //     << "bPoint:" << bHits[bHitI].point()
            //     << "  tracking:" << singleParticle.position()
            //     << "  dist:" << dist
            //     << endl;

            if (dist > smallDist)
            {
                // Hit-point is past tracking position
                foundValidB = true;
                break;
            }
            else
            {
                ++bHitI;
            }
        }

        if (!foundValidB)
        {
            // No valid boundary intersection found beyond tracking position
            break;
        }

        // Update starting point for tracking
        trackFacei = bHits[bHitI].index();
        trackPt = pushIn(bHits[bHitI].hitPoint(), trackFacei);
        trackCelli = getBoundaryCell(trackFacei);

        ++segmentI;

        startSegmentI = samplingPts.size();
    }

    const_cast<polyMesh&>(mesh()).moving(oldMoving);
}


void Foam::faceOnlySet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceOnlySet::faceOnlySet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const point& start,
    const point& end
)
:
    sampledSet(name, mesh, searchEngine, axis),
    start_(start),
    end_(end)
{
    genSamples();
}


Foam::faceOnlySet::faceOnlySet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    start_(dict.get<point>("start")),
    end_(dict.get<point>("end"))
{
    genSamples();
}


// ************************************************************************* //
