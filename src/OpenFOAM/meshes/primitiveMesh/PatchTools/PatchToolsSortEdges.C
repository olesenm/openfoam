/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "PatchTools.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::labelListList
Foam::PatchTools::sortedEdgeFaces
(
    const PrimitivePatch<FaceList, PointField>& p
)
{
    using point_type = typename PointField::value_type;

    const edgeList& edges = p.edges();
    const labelListList& edgeFaces = p.edgeFaces();
    const auto& localFaces = p.localFaces();

    // const auto& localPoints = p.localPoints();
    // The same without the expense of demand-driven localPoints()
    const UIndirectList<point_type> localPoints(p.points(), p.meshPoints());

    // create the lists for the various results. (resized on completion)
    labelListList sortedEdgeFaces(edgeFaces.size());

    // Work array
    DynamicList<scalar> faceAngles;

    forAll(edgeFaces, edgeI)
    {
        const labelList& faceNbs = edgeFaces[edgeI];

        if (faceNbs.size() > 2)
        {
            // Get point on edge and normalized direction of edge (= e2 base
            // of our coordinate system)
            const edge& e = edges[edgeI];

            const point& edgePt = localPoints[e.start()];

            const vector e2 = normalised(localPoints[e.end()] - edgePt);

            scalar maxAngle = GREAT;
            vector maxAngleEdgeDir(vector::max);

            {
                // Get the vertex on 0th face that forms a vector with the first
                // edge point that has the largest angle with the edge

                constexpr label nbI = 0;

                for (label fpi : localFaces[faceNbs[nbI]])
                {
                    if (fpi != e.start())
                    {
                        const vector faceEdgeDir =
                            normalised(localPoints[fpi] - edgePt);

                        const scalar angle = e2 & faceEdgeDir;

                        if (Foam::mag(angle) < maxAngle)
                        {
                            maxAngle = angle;
                            maxAngleEdgeDir = faceEdgeDir;
                        }
                    }
                }
            }

            // Get vector normal both to e2 and to edge from opposite vertex
            // to edge (will be x-axis of our coordinate system)
            const vector e0 = normalised(e2 ^ maxAngleEdgeDir);

            // Get y-axis of coordinate system
            const vector e1 = e2 ^ e0;

            // e0 is reference so angle is 0
            faceAngles.clear();
            faceAngles.push_back(0);

            for (label nbI = 1; nbI < faceNbs.size(); nbI++)
            {
                // Get the vertex on face that forms a vector with the first
                // edge point that has the largest angle with the edge

                maxAngle = GREAT;
                maxAngleEdgeDir = vector::max;

                for (label fpi : localFaces[faceNbs[nbI]])
                {
                    if (fpi != e.start())
                    {
                        const vector faceEdgeDir =
                            normalised(localPoints[fpi] - edgePt);

                        const scalar angle = e2 & faceEdgeDir;

                        if (Foam::mag(angle) < maxAngle)
                        {
                            maxAngle = angle;
                            maxAngleEdgeDir = faceEdgeDir;
                        }
                    }
                }

                const vector vec = normalised(e2 ^ maxAngleEdgeDir);

                faceAngles.push_back(pseudoAngle(e0, e1, vec));
            }

            // Determine the sort order of the angles directly into the
            // destination slot, then transcribe in-place from angle-sort
            // index to the corresponding face label
            {
                auto& slot = sortedEdgeFaces[edgeI];
                Foam::sortedOrder(faceAngles, slot);

                for (label& idx : slot)
                {
                    idx = faceNbs[idx];
                }
            }
        }
        else
        {
            // No need to sort. Just copy.
            sortedEdgeFaces[edgeI] = faceNbs;
        }
    }

    return sortedEdgeFaces;
}


// ************************************************************************* //
