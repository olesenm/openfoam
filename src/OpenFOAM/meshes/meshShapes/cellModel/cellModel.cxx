/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017,2025 OpenCFD Ltd.
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

#include "cellModel.H"
#include "pyramid.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Average of points
// Note: use double precision to avoid overflows when summing
static inline Foam::doubleVector pointsAverage
(
    const UList<point>& points,
    const labelUList& pointLabels
)
{
    doubleVector avg(Zero);

    if (const auto n = pointLabels.size(); n)
    {
        for (const auto pointi : pointLabels)
        {
            avg += points[pointi];
        }
        avg /= n;
    }

    return avg;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::cellModel::centre
(
    const labelUList& pointLabels,
    const UList<point>& points
) const
{
    // Note: use double precision to avoid overflows when summing

    // Estimated centre by averaging the cell points
    const point centrePoint(pointsAverage(points, pointLabels));


    // Calculate the centre by breaking the cell into pyramids and
    // volume-weighted averaging their centres

    doubleScalar sumV(0);
    doubleVector sumVc(Zero);

    forAll(faces_, facei)
    {
        const Foam::face f(pointLabels, faces_[facei]);

        const scalar pyrVol = pyramidPointFaceRef(f, centrePoint).mag(points);

        if (pyrVol > SMALL)
        {
            WarningInFunction
                << "zero or negative pyramid volume: " << -pyrVol
                << " for face " << facei
                << endl;
        }

        sumV -= pyrVol;
        sumVc -= pyrVol * pyramidPointFaceRef(f, centrePoint).centre(points);
    }

    return sumVc/(sumV + VSMALL);
}


Foam::scalar Foam::cellModel::mag
(
    const labelUList& pointLabels,
    const UList<point>& points
) const
{
    // Note: use double precision to avoid overflows when summing

    // Estimated centre by averaging the cell points
    const point centrePoint(pointsAverage(points, pointLabels));


    // Calculate the magnitude by summing the -mags of the pyramids
    // The sign change is because the faces point outwards
    // and a pyramid is constructed from an inward pointing face
    // and the base centre-apex vector

    scalar sumV(0);

    forAll(faces_, facei)
    {
        const Foam::face f(pointLabels, faces_[facei]);

        const scalar pyrVol = pyramidPointFaceRef(f, centrePoint).mag(points);

        if (pyrVol > SMALL)
        {
            WarningInFunction
                << "zero or negative pyramid volume: " << -pyrVol
                << " for face " << facei
                << endl;
        }

        sumV -= pyrVol;
    }

    return sumV;
}


// ************************************************************************* //
