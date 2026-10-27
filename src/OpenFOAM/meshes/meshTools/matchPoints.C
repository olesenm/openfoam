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

#include "matchPoints.H"
#include "SortableList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<scalar>& matchDistances,
    const bool verbose,
    labelList& from0To1,
    const point& origin
)
{
    from0To1.resize_fill(pts0.size(), -1);

    if (pts0.empty())
    {
        // An empty subset, so no mismatches are possible
        return true;
    }

    // Sanity checks
    if (matchDistances.empty())
    {
        FatalErrorInFunction
            << "No match distances provided for any of the "
            << pts0.size() << " points" << endl
            << exit(FatalError);
    }

    // The last addressable value in matchDistances
    const scalar endMatchDistance = matchDistances.back();

    bool fullMatch = true;

    point compareOrigin = origin;

    if (origin == point::max)
    {
        if (pts1.size())
        {
            compareOrigin = Foam::sum(pts1)/pts1.size();
        }
        else
        {
            // Unusable, but avoid comparison with VGREAT!
            compareOrigin = point::zero;
        }
    }

    // Fill directly to avoid a temporary (point - origin) field
    SortableList<scalar> pts0MagSqr(pts0.size());
    forAll(pts0, i)
    {
        pts0MagSqr[i] = pts0[i].distSqr(compareOrigin);
    }
    pts0MagSqr.sort();

    SortableList<scalar> pts1MagSqr(pts1.size());
    forAll(pts1, i)
    {
        pts1MagSqr[i] = pts1[i].distSqr(compareOrigin);
    }
    pts1MagSqr.sort();


    forAll(pts0MagSqr, i)
    {
        const auto dist0 = pts0MagSqr[i];
        const auto index0 = pts0MagSqr.indices()[i];

        const scalar matchDist =
        (
            (index0 < matchDistances.size())
          ? matchDistances[index0]
          : endMatchDistance
        );
        const scalar matchDistSqr = Foam::sqr(matchDist);

        label startI = findLower(pts1MagSqr, 0.99999*dist0 - 2*matchDist);

        if (startI == -1)
        {
            startI = 0;
        }

        // Go through range of equal mag and find nearest vector.
        scalar minDistSqr = VGREAT;
        label minIndex1 = -1;

        for
        (
            label j = startI;
            (
                (j < pts1MagSqr.size())
             && (pts1MagSqr[j] < 1.00001*dist0 + 2*matchDist)
            );
            j++
        )
        {
            label index1 = pts1MagSqr.indices()[j];

            // Compare actual vectors
            scalar distSqr = pts0[index0].distSqr(pts1[index1]);

            if (distSqr <= matchDistSqr && distSqr < minDistSqr)
            {
                minDistSqr = distSqr;
                minIndex1 = index1;
            }
        }

        if (minIndex1 == -1)
        {
            fullMatch = false;

            if (verbose)
            {
                Pout<< "Cannot find point in pts1 matching point " << index0
                    << " coord:" << pts0[index0]
                    << " in pts0 when using tolerance " << matchDist << endl;

                // Go through range of equal mag and find equal vector.
                Pout<< "Searching started from:" << startI << " in pts1" << nl;

                for
                (
                    label j = startI;
                    (
                        (j < pts1MagSqr.size())
                     && (pts1MagSqr[j] < 1.00001*dist0 + 2*matchDist)
                    );
                    j++
                )
                {
                    label index1 = pts1MagSqr.indices()[j];

                    Pout<< "    Compared coord: " << pts1[index1]
                        << " at index " << j
                        << " with difference to point "
                        << pts0[index0].dist(pts1[index1]) << endl;
                }
            }
        }

        from0To1[index0] = minIndex1;
    }

    return fullMatch;
}


bool Foam::matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<point>& pts0Dir,
    const UList<point>& pts1Dir,
    const UList<scalar>& matchDistances,
    const bool verbose,
    labelList& from0To1,
    const point& origin
)
{
    from0To1.resize_fill(pts0.size(), -1);

    if (pts0.empty())
    {
        // An empty subset, so no mismatches are possible
        return true;
    }

    // Sanity checks
    if (matchDistances.empty())
    {
        FatalErrorInFunction
            << "No match distances provided for any of the "
            << pts0.size() << " points" << endl
            << exit(FatalError);
    }
    if
    (
        (pts0Dir.size() < pts0.size())
     || (pts1Dir.size() < pts1.size())
    )
    {
        FatalErrorInFunction
            << "Points 0 with " << pts0Dir.size() << "/" << pts0.size()
            << " directions" << nl
            << "Points 1 with " << pts1Dir.size() << "/" << pts1.size()
            << " directions" << endl
            << exit(FatalError);
    }


    // The last addressable value in matchDistances
    const scalar endMatchDistance = matchDistances.back();

    bool fullMatch = true;

    point compareOrigin = origin;

    if (origin == point::max)
    {
        if (pts1.size())
        {
            compareOrigin = Foam::sum(pts1)/pts1.size();
        }
        else
        {
            // Unusable, but avoid comparison with VGREAT!
            compareOrigin = point::zero;
        }
    }


    // Fill directly to avoid a temporary (point - origin) field
    SortableList<scalar> pts0MagSqr(pts0.size());
    forAll(pts0, i)
    {
        pts0MagSqr[i] = pts0[i].distSqr(compareOrigin);
    }
    pts0MagSqr.sort();

    SortableList<scalar> pts1MagSqr(pts1.size());
    forAll(pts1, i)
    {
        pts1MagSqr[i] = pts1[i].distSqr(compareOrigin);
    }
    pts1MagSqr.sort();


    forAll(pts0MagSqr, i)
    {
        const auto dist0 = pts0MagSqr[i];
        const auto index0 = pts0MagSqr.indices()[i];

        const scalar matchDist =
        (
            (index0 < matchDistances.size())
          ? matchDistances[index0]
          : endMatchDistance
        );
        const scalar matchDistSqr = Foam::sqr(matchDist);

        label startI = findLower(pts1MagSqr, 0.99999*dist0 - 2*matchDist);

        if (startI == -1)
        {
            startI = 0;
        }

        // Go through range of equal mag and find nearest vector.
        scalar minDistSqr = VGREAT;
        label minIndex1 = -1;

        // Valid candidate points should have opposite normal
        const scalar minDistNorm = 0;

        for
        (
            label j = startI;
            (
                (j < pts1MagSqr.size())
             && (pts1MagSqr[j] < 1.00001*dist0 + 2*matchDist)
            );
            j++
        )
        {
            label index1 = pts1MagSqr.indices()[j];

            // Compare actual vectors
            scalar distSqr = pts0[index0].distSqr(pts1[index1]);
            scalar distNorm = (pts0Dir[index0] & pts1Dir[index1]);

            if
            (
                pts0Dir[index0].magSqr() < Foam::sqr(SMALL)
             && pts1Dir[index1].magSqr() < Foam::sqr(SMALL)
            )
            {
                distNorm = -1;
            }

            if (distSqr <= matchDistSqr && distSqr < minDistSqr)
            {
                // Check that the normals point in equal and opposite directions
                if (distNorm < minDistNorm)
                {
                    minDistSqr = distSqr;
                    minIndex1 = index1;
                }
            }
        }

        if (minIndex1 == -1)
        {
            fullMatch = false;

            if (verbose)
            {
                Pout<< "Cannot find point in pts1 matching point " << index0
                    << " coord:" << pts0[index0]
                    << " in pts0 when using tolerance " << matchDist << endl;

                // Go through range of equal mag and find equal vector.
                Pout<< "Searching started from:" << startI << " in pts1"
                    << endl;
                for
                (
                    label j = startI;
                    (
                        (j < pts1MagSqr.size())
                     && (pts1MagSqr[j] < 1.00001*dist0 + 2*matchDist)
                    );
                    j++
                )
                {
                    label index1 = pts1MagSqr.indices()[j];

                    Pout<< "    Compared coord: " << pts1[index1]
                        << " at index " << j
                        << " with difference to point "
                        << pts0[index0].dist(pts1[index1]) << endl;
                }
            }
        }

        from0To1[index0] = minIndex1;
    }

    return fullMatch;
}


bool Foam::matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    scalar matchDistance,
    const bool verbose,
    labelList& from0To1,
    const point& origin
)
{
    // Wrap match distance as single-value UList
    scalar distanceValue = matchDistance;
    const UList distances(&distanceValue, 1);

    return matchPoints
    (
        pts0,
        pts1,
        distances,
        verbose,
        from0To1,
        origin
    );
}


bool Foam::matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<point>& pts0Dir,
    const UList<point>& pts1Dir,
    scalar matchDistance,
    const bool verbose,
    labelList& from0To1,
    const point& origin
)
{
    // Wrap match distance as single-value UList
    scalar distanceValue = matchDistance;
    const UList distances(&distanceValue, 1);

    return matchPoints
    (
        pts0,
        pts1,
        pts0Dir,
        pts1Dir,
        distances,
        verbose,
        from0To1,
        origin
    );
}


// ************************************************************************* //
