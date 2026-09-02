/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::AABBTreeBase::nLeafs() const noexcept
{
    const label len = (offsets_.size() - 1);
    return (len < 1) ? 0 : len;
}


inline bool Foam::AABBTreeBase::hasNodes() const noexcept
{
    return bool(nodes_);
}


inline const Foam::labelList::subList
Foam::AABBTreeBase::indices(label leaf) const
{
    if (leaf >= 0 && leaf < offsets_.size()-1)
    {
        return indices_.slice
        (
            offsets_[leaf],
            (offsets_[leaf+1]-offsets_[leaf])
        );
    }
    else
    {
        return labelList::subList();
    }
}


inline Foam::labelRange Foam::AABBTreeBase::range(label leaf) const
{
    if (leaf >= 0 && leaf < offsets_.size()-1)
    {
        return labelRange
        (
            offsets_[leaf],
            (offsets_[leaf+1]-offsets_[leaf])
        );
    }
    else
    {
        return labelRange();
    }
}


// ************************************************************************* //
