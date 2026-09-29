/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "boolList.H"
#include "CircularBuffer.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Orders the local points on the patch for most efficient search

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcLocalPointOrder() const
{
    // Note: Cannot use bandCompression as point-point addressing does
    // not exist and is not considered generally useful.

    DebugInFunction << "Calculating local point order" << endl;

    if (localPointOrderPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "local point order already calculated"
            << abort(FatalError);
    }

    const auto& lf = localFaces();

    const labelListList& ff = faceFaces();

    localPointOrderPtr_.reset(new labelList(meshPoints().size(), -1));
    auto& pointOrder = *localPointOrderPtr_;

    boolList visitedFace(lf.size(), false);
    boolList visitedPoint(pointOrder.size(), false);

    label nPoints = 0;

    // FIFO buffer managing point/face insertion order
    CircularBuffer<label> faceOrder(32);

    forAll(lf, facei)
    {
        if (!visitedFace[facei])
        {
            faceOrder.push_back(facei);

            while (!faceOrder.empty())
            {
                // Process as FIFO
                const label curFace = faceOrder.front();
                faceOrder.pop_front();

                if (!visitedFace[curFace])
                {
                    visitedFace[curFace] = true;

                    // mark points
                    for (const label pointi : lf[curFace])
                    {
                        if (!visitedPoint[pointi])
                        {
                            visitedPoint[pointi] = true;

                            pointOrder[nPoints] = pointi;

                            ++nPoints;
                        }
                    }

                    // Add unvisited face neighbours to the list

                    for (const label nbrFacei : ff[curFace])
                    {
                        if (!visitedFace[nbrFacei])
                        {
                            faceOrder.push_back(nbrFacei);
                        }
                    }
                }
            }
        }
    }

    DebugInfo << "Calculated local point order" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
template<class OutputPointFieldType>
void Foam::PrimitivePatch<FaceList, PointField>::localPoints
(
    Foam::identity_construct_t,
    OutputPointFieldType& newPoints
) const
{
    static_assert
    (
        std::is_same_v
        <
            typename OutputPointFieldType::value_type,
            point_type
        >,
        "output list value_type must match the patch point_type"
    );

    // Copy out all local points, avoid demand-driven localPoints().
    // They are based on meshPoints() so use that as the least costly solution.

    if (points_.empty() || this->nFaces() == 0)
    {
        newPoints.clear();
        return;
    }

    // The local to global mapping
    const auto& meshPts = meshPoints();

    if (label nLocalPoints = meshPts.size(); nLocalPoints > 0)
    {
        newPoints.resize_nocopy(nLocalPoints);
    }
    else
    {
        newPoints.clear();
        return;
    }

    auto iter = newPoints.begin();

    if (hasLocalPoints())
    {
        const auto& locPoints = this->localPoints();

        std::copy(locPoints.begin(), locPoints.end(), iter);
    }
    else
    {
        // Copy out all of the local points
        for (label globalPointi : meshPts)
        {
            *iter = points_[globalPointi];
            ++iter;
        }
    }
}


template<class FaceList, class PointField>
template<class OutputPointFieldType>
void Foam::PrimitivePatch<FaceList, PointField>::localPoints
(
    const labelUList& pointMap,
    OutputPointFieldType& newPoints
) const
{
    static_assert
    (
        std::is_same_v
        <
            typename OutputPointFieldType::value_type,
            point_type
        >,
        "output list value_type must match the patch point_type"
    );

    // Copy out a sub-sample of local points, avoid demand-driven localPoints().
    // They are based on meshPoints() so use that as the least costly solution.

    if (pointMap.empty())
    {
        newPoints.clear();
        return;
    }

    newPoints.resize_nocopy(pointMap.size());

    auto iter = newPoints.begin();

    if (hasLocalPoints())
    {
        const auto& locPoints = this->localPoints();

        // Subset the local points by the pointMap
        for (label localPointi : pointMap)
        {
            *iter = locPoints[localPointi];
            ++iter;
        }
    }
    else
    {
        // The local to global mapping
        const auto& meshPts = meshPoints();

        // Subset the local points by the pointMap
        for (label localPointi : pointMap)
        {
            *iter = points_[meshPts[localPointi]];
            ++iter;
        }
    }
}


template<class FaceList, class PointField>
Foam::Field<typename Foam::PrimitivePatch<FaceList, PointField>::point_type>
Foam::PrimitivePatch<FaceList, PointField>::localPoints
(
    const labelUList& pointMap
) const
{
    Field<point_type> newPoints;

    this->localPoints(pointMap, newPoints);

    return newPoints;
}


// ************************************************************************* //
