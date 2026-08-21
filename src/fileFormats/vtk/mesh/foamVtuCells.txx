/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025-2026 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PointType>
Foam::tmp<Foam::Field<PointType>>
Foam::vtk::vtuCells::points(const polyMesh& mesh) const
{
    typedef Foam::Field<PointType> pointFieldType;

    // Combine mesh points and any additional cellCentre points
    // into a single field

    const auto& pts = mesh.points();
    const auto& cc = mesh.cellCentres();

    // The additional cellCentre points
    const labelUList& addPoints = addPointCellLabels();

    if constexpr (std::is_same_v<Foam::point, PointType>)
    {
        if (addPoints.empty() && !this->usePointMap())
        {
            // No decomposed cells etc
            return mesh.points();
        }
    }

    if (const auto& map = this->pointMap(); !map.empty())
    {
        // With pointMap
        auto tpoints = tmp<pointFieldType>::New(map.size() + addPoints.size());

        auto iter = tpoints.ref().begin();

        // Normal points (unique mesh points)
        for (const auto pointi : map)
        {
            *iter = pts[pointi];
            ++iter;
        }

        // Cell centres
        for (const label celli : addPoints)
        {
            *iter = cc[celli];
            ++iter;
        }

        return tpoints;
    }
    else
    {
        // All mesh points

        auto tpoints = tmp<pointFieldType>::New(pts.size() + addPoints.size());

        auto iter = tpoints.ref().begin();

        // Normal points
        iter = std::copy(pts.begin(), pts.end(), iter);

        // Cell centres
        for (const label celli : addPoints)
        {
            *iter = cc[celli];
            ++iter;
        }

        return tpoints;
    }
}


// ************************************************************************* //
