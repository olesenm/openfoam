/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "cellShape.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

int Foam::cellShape::compare(const cellShape& a, const cellShape& b)
{
    // Basic rule: we assume that the sequence of labels in each list
    // will be circular in the same order (but not necessarily in the
    // same direction). The limitation of this method is that with 3D
    // topologies I cannot guarantee that a congruent but not
    // identical cellShape (i.e. one sharing the same points but in a
    // different order) will necessarily be matched.

    const labelUList& labsA = a;
    const labelUList& labsB = b;

    const label sizeA = labsA.size();
    const label sizeB = labsB.size();

    if (sizeA != sizeB)
    {
        // Trivial reject: different sizes
        return 0;
    }
    else if (sizeA == 0)
    {
        // Both have zero vertices. Always identical
        return 1;
    }
    else if (sizeA == 1)
    {
        // Both have a single vertex? Simple check
        return (labsA[0] == labsB[0] ? 1 : 0);
    }


    // Search B for the starting point of A
    label Bptr = labsB.find(labsA[0]);

    if (Bptr < 0)
    {
        // No match found
        return 0;
    }

    // Now check B for the next label of A,
    // it could be in either direction
    label Aptr = 1;
    int dir = 0;

    ++Bptr;
    if (Bptr == labsB.size())
    {
        Bptr = 0;
    }

    if (labsA[Aptr] == labsB[Bptr])
    {
        // Yes - direction is 'up'
        dir = 1;
    }
    else
    {
        // No - so look downwards
        Bptr -= 2;
        if (Bptr < 0)
        {
            Bptr += labsB.size();
        }

        // Test whether downward label matches the second label in A
        if (labsA[Aptr] == labsB[Bptr])
        {
            // Yes - direction is 'down'
            dir = -1;
        }
    }

    // Check whether a match was made at all, and exit false if not
    if (dir == 0)
    {
        return 0;
    }

    // We now have both direction of search and next element
    // to search, so we can continue search until no more points.
    // Decrement size by 2 to account for first searches

    label remaining = (sizeA - 2);

    if (dir > 0)
    {
        while (remaining--)
        {
            ++Aptr;
            if (Aptr >= labsA.size())
            {
                Aptr = 0;
            }

            ++Bptr;
            if (Bptr >= labsB.size())
            {
                Bptr = 0;
            }

            if (labsA[Aptr] != labsB[Bptr])
            {
                return 0;
            }
        }
    }
    else
    {
        while (remaining--)
        {
            ++Aptr;
            if (Aptr >= labsA.size())
            {
                Aptr = 0;
            }

            --Bptr;
            if (Bptr < 0)
            {
                Bptr = labsB.size() - 1;
            }

            if (labsA[Aptr] != labsB[Bptr])
            {
                return 0;
            }
        }
    }

    // Equal values but perhaps different order
    return dir;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::cellShape::operator==(const cellShape& other) const
{
    return cellShape::compare(*this, other) != 0;
}

// bool Foam::cellShape::operator!=(const cellShape& other) const
// {
//     return cellShape::compare(*this, other) == 0;
// }


// ************************************************************************* //
