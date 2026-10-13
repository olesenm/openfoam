/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2020 OpenFOAM Foundation
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

#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "globalMeshData.H"
#include "pointMesh.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::updateMesh(const mapPolyMesh& mpm)
{
    DebugInFunction
        << "Updating addressing and (optional) pointMesh/pointFields" << endl;

    // Update boundaryMesh (note that patches themselves already ok)
//    boundary_.updateMesh(mpm);
    boundary_.updateMesh();

    // Update zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    // Remove the stored tet base points
    tetBasePtIsPtr_.reset(nullptr);

    // Remove the cell tree
    cellTreePtr_.reset(nullptr);

    // Update parallel data
    if (globalMeshDataPtr_)
    {
        globalMeshDataPtr_->updateMesh();
    }

    setInstance(time().timeName());

    // Map the old motion points if present
    if (oldPointsPtr_)
    {
        // Make a copy of the original points
        pointField oldMotionPoints(*oldPointsPtr_);

        pointField& newMotionPoints = *oldPointsPtr_;

        // Resize the list to new size
        newMotionPoints.resize(points_.size());

        // Map the list
        if (mpm.hasMotionPoints())
        {
            const auto& pointMap = mpm.pointMap();

            // Any points created out-of-nothing get set to the current
            // coordinate for lack of anything better.
            forAll(pointMap, newPointi)
            {
                if
                (
                    label oldPointi = pointMap[newPointi];
                    (oldPointi >= 0 && oldPointi < oldMotionPoints.size())
                )
                {
                    newMotionPoints[newPointi] = oldMotionPoints[oldPointi];
                }
                else
                {
                    newMotionPoints[newPointi] = points_[newPointi];
                }
            }
        }
        else
        {
            const labelList& pointMap = mpm.pointMap();
            const labelList& revPointMap = mpm.reversePointMap();

            forAll(pointMap, newPointi)
            {
                label oldPointi = pointMap[newPointi];
                if
                (
                    oldPointi >= 0
                 && oldPointi < revPointMap.size()
                 && oldPointi < oldMotionPoints.size()
                 && revPointMap[oldPointi] == newPointi
                ) // master point
                {
                    newMotionPoints[newPointi] = oldMotionPoints[oldPointi];
                }
                else
                {
                    newMotionPoints[newPointi] = points_[newPointi];
                }
            }
        }
    }

    // Map the old motion cell-centres if present
    if (oldCellCentresPtr_)
    {
        // Make a copy of the original cell-centres
        pointField oldMotionCellCentres(*oldCellCentresPtr_);

        pointField& newMotionCellCentres = *oldCellCentresPtr_;

        // Resize the list to new size
        newMotionCellCentres.resize(cellCentres().size());

        // Any points created out-of-nothing get set to the current coordinate
        // for lack of anything better.
        const auto& cellMap = mpm.cellMap();
        forAll(cellMap, celli)
        {
            if
            (
                label oldCelli = cellMap[celli];
                (oldCelli >= 0 && oldCelli < oldMotionCellCentres.size())
            )
            {
                newMotionCellCentres[celli] = oldMotionCellCentres[oldCelli];
            }
            else
            {
                newMotionCellCentres[celli] = cellCentres()[celli];
            }
        }
    }

    meshObject::updateMesh<polyMesh>(*this, mpm);
    meshObject::updateMesh<pointMesh>(*this, mpm);

    // Reset valid directions (could change by faces put into empty patches)
    geometricD_ = Zero;
    solutionD_ = Zero;

    const_cast<Time&>(time()).functionObjects().updateMesh(mpm);
}


// ************************************************************************* //
