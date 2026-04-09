/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Class
    Foam::zoneMotion

Description
    Intermediate class for handling "zonified" motion.
    Can select motion for the entire mesh, or restrict to a subset
    of points based on a specified cellSet or cellZone(s).

Usage
    \verbatim
    {
        cellSet     <name>;

        // and/or
        cellZone    <name or regex>;
    }
    \endverbatim

Note
    The cellZone can be specified as a wordRe and will select by name
    or group. The special name "none" is treated as ignored.

SourceFiles
    zoneMotion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zoneMotion_H
#define Foam_zoneMotion_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class polyMesh;
class wordRe;

/*---------------------------------------------------------------------------*\
                           Class zoneMotion Declaration
\*---------------------------------------------------------------------------*/

class zoneMotion
{
public:

    // Public Data

        //- The selection modes
        enum class selectionModes : char
        {
            smAll = 0,              //!< "all" cells
            smCellSet = 1,          //!< "cellSet"
            smCellZone = 2,         //!< "cellZone"
            smCellSetAndZone = 3    //!< both cellSet and cellZone specified
        };

private:

     // Private Data

        //- Points to move when cell zones or sets are supplied
        labelList pointIDs_;

        //- Selected cell zones (for 'cellZone' selection mode)
        labelList cellZoneIDs_;

        //- Cell selection mode
        selectionModes selectionMode_;


    // Private Member Functions

        //- Select pointIDs according to the current cellZoneIDs
        void selectZonePoints(const polyMesh& mesh);


public:

    // Constructors

        //- Default construct (all points)
        zoneMotion();

        //- Construct from dictionary (uses cellSet, cellZone entries)
        zoneMotion(const polyMesh& mesh, const dictionary& dict);

        //- Construct from dictionary (uses cellSet, cellZone entries)
        zoneMotion(const dictionary& dict, const polyMesh& mesh);

        //- Construct from cellZones specified by ids.
        //- An empty list corresponds to "all"
        zoneMotion(const polyMesh& mesh, const labelUList& cellZoneIds);

        //- Construct from cellZones specified by name.
        //- An empty name corresponds to "all"
        zoneMotion(const polyMesh& mesh, const wordRe& cellZoneSelection);


    //- Destructor
    virtual ~zoneMotion() = default;


    // Member Functions

        //- The point ids (for cell set/zone subset)
        const labelUList& pointIDs() const noexcept { return pointIDs_; }

        //- Selected cell zones (if any)
        const labelUList& cellZoneIDs() const noexcept { return cellZoneIDs_; }

        //- Move all cells?
        bool moveAllCells() const noexcept
        {
            return (selectionMode_ == selectionModes::smAll);
        }

        //- Reset to use all points
        void resetZone();

        //- Reset with cellZones specified by name.
        //- An empty name corresponds to "all"
        void resetZone(const polyMesh& mesh, const wordRe& cellZoneSelection);

        //- Reset with cellZones specified by IDs
        void resetZone(const polyMesh& mesh, const labelUList& cellZoneIds);

        //- Reset from dictionary (uses cellSet, cellZone entries)
        void resetZone(const polyMesh& mesh, const dictionary& dict);

        //- Rebuild pointIDs_ and selection state from the current mesh
        void rebuildPointIDs(const dictionary& dict, const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
