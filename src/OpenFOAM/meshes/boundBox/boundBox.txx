/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "boundBox.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned N>
Foam::boundBox::boundBox
(
    const UList<point>& points,
    const FixedList<label, N>& indices,
    bool syncPar
)
:
    boundBox()
{
    add(points, indices);

    if (syncPar)
    {
        this->reduce();   // Uses UPstream::worldComm
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<unsigned N>
void Foam::boundBox::add
(
    const FixedList<point, N>& points
)
{
    for (const point& p : points)
    {
        add(p);
    }
}


template<class IndexContainer>
void Foam::boundBox::add
(
    const UList<point>& points,
    const IndexContainer& indices
)
{
    const label len = points.size();

    // Skip if points is empty
    if (!len)
    {
        return;
    }

    for (label i : indices)
    {
        // Ignore out-of-bound indices
        if (i >= 0 && i < len)
        {
            add(points[i]);
        }
    }
}


template<class IndexContainer>
inline bool Foam::boundBox::contains_all
(
    const UList<point>& points,
    const IndexContainer& indices
) const
{
    const label len = points.size();

    if (!len || indices.empty())
    {
        // Empty range == true (as per std::all_of).
        // ie, no test on individual elements can fail
        return true;
    }
    else if (!good())
    {
        // Inverted box - no match possible
        return false;
    }

    for (label i : indices)
    {
        // Treat out-of-bound indices as not part of the set
        if (i >= 0 && i < len)
        {
            if (!contains(points[i]))
            {
                return false;
            }
        }
    }

    return true;
}


template<class IndexContainer>
inline bool Foam::boundBox::contains_any
(
    const UList<point>& points,
    const IndexContainer& indices
) const
{
    const label len = points.size();

    if (!len || indices.empty())
    {
        // Trivial case - no match possible
        return false;
    }
    else if (!good())
    {
        // Inverted box - no match possible
        return false;
    }

    for (label i : indices)
    {
        // Treat out-of-bound indices as not part of the set
        if (i >= 0 && i < len)
        {
            if (contains(points[i]))
            {
                return true;
            }
        }
    }

    return false;
}


// ************************************************************************* //
