/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

Description
    Determine correspondence between points. See below.

SourceFiles
    matchPoints.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_matchPoints_H
#define Foam_matchPoints_H

#include "scalarField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Function matchPoints Declaration
\*---------------------------------------------------------------------------*/

//- Determine correspondence between pointFields.
//- Gets passed local matching distance per point.
//- Returns map from 0 to 1.
//- Returns true if every point in 0 has a corresponding point in 1;
//- false otherwise.
//- Compares mag(pt - origin) and does proper comparison for points
//- with (almost) identical mag. pts0 can be subset of pts1.
bool matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    //! The match distances. If the array is shorter than \p pts0, uses its
    //! last value. Can thus have an array of size 1 as \em uniform.
    const UList<scalar>& matchDistances,
    const bool verbose,
    //! [out] Map of 0 to 1
    labelList& from0To1,
    const point& origin = point::zero  //!< point::max = use average origin
);


//- Determine correspondence between pointFields using the supplied
//- \p pts0Dir and \p pts1Dir, which are the directions associated with
//- the points (eg, face normal associated with each face centre).
//- A match between similar points is only allowed if their directions are
//- equal and opposite
bool matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<point>& pts0Dir,
    const UList<point>& pts1Dir,
    //! The match distance. If the array is shorter than \p pts0, uses its
    //! last value. Can thus have an array of size 1 as \em uniform.
    const UList<scalar>& matchDistance,
    const bool verbose,
    //! [out] Map of 0 to 1
    labelList& from0To1,
    const point& origin = point::zero  //!< point::max = use average origin
);


//- Determine correspondence between pointFields.
//  The regular matchPoints algorithm but with an identical
//  matching distance for all points.
bool matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    //! The match distance (identical for all matches)
    scalar matchDistance,
    const bool verbose,
    //! [out] Map of 0 to 1
    labelList& from0To1,
    const point& origin = point::zero  //!< point::max = use average origin
);


//- Determine correspondence between pointFields using the supplied
//- \p pts0Dir and \p pts1Dir.
//  The regular matchPoints algorithm but with an identical
//  matching distance for all points.
bool matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<point>& pts0Dir,
    const UList<point>& pts1Dir,
    //! The match distance (identical for all matches)
    scalar matchDistance,
    const bool verbose,
    //! [out] Map of 0 to 1
    labelList& from0To1,
    const point& origin = point::zero  //!< point::max = use average origin
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
