/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "ListOps.H"  // For uniqueOrder()

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SortList<T>::SortList(const UList<T>& values)
:
    parent_type(values, labelList())
{
    sort();
}


template<class T>
template<class Compare>
inline Foam::SortList<T>::SortList(const UList<T>& values, const Compare& comp)
:
    parent_type(values, labelList())
{
    sort<Compare>(comp);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::labelUList& Foam::SortList<T>::indices() const noexcept
{
    return this->addressing();
}


template<class T>
inline Foam::labelList& Foam::SortList<T>::indices() noexcept
{
    return this->addressing();
}


template<class T>
inline void Foam::SortList<T>::reverse()
{
    Foam::reverse(this->indices());
}


template<class T>
inline void Foam::SortList<T>::reset()
{
    this->indices().resize_nocopy(this->values().size());
    Foam::identity(this->indices());
}


template<class T>
template<class Compare>
inline void Foam::SortList<T>::sort(const Compare& comp)
{
    auto& vals = this->values();

    Foam::sortedOrder
    (
        vals,
        this->indices(),
        [&](label a, label b) -> bool { return comp(vals[a], vals[b]); }
    );
}


template<class T>
inline void Foam::SortList<T>::sort()
{
    Foam::sortedOrder(this->values(), this->indices());
}


template<class T>
inline void Foam::SortList<T>::uniqueSort()
{
    Foam::uniqueOrder(this->values(), this->indices());
}


template<class T>
inline void Foam::SortList<T>::reverseSort()
{
    Foam::sortedOrder
    (
        this->values(),
        this->indices(),
        ListOps::greater(this->values())
    );
}


// ************************************************************************* //
