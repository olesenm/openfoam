/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(label len)
:
    List<T>(len)
{}


template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(label len, Foam::zero)
:
    List<T>(len, Foam::zero{})
{}


template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(label len, const T& val)
:
    List<T>(len, val)
{}


template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(const SortableList<T, Forward>& lst)
:
    List<T>(lst),
    indices_(lst.indices())
{}


template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(SortableList<T, Forward>&& lst)
:
    List<T>(std::move(lst)),
    indices_(std::move(lst.indices_))
{}


template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(const UList<T>& values)
:
    List<T>(values)
{
    sort();
}


template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(List<T>&& values)
:
    List<T>(std::move(values))
{
    sort();
}


template<class T, bool Forward>
inline Foam::SortableList<T, Forward>::
SortableList(std::initializer_list<T> values)
:
    List<T>(values)
{
    sort();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::clear()
{
    List<T>::clear();
    indices_.clear();
}


template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::shrink_to_fit()
{
    indices_.clear();
}


template<class T, bool Forward>
inline Foam::List<T>&
Foam::SortableList<T, Forward>::shrink()
{
    indices_.clear();
    return static_cast<List<T>&>(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, bool Forward>
void Foam::SortableList<T, Forward>::sort()
{
    Foam::sortedOrder(*this, indices_, normalCompareOp());

    List<T> list(*this, indices_); // Copy with indices for mapping
    List<T>::transfer(list);
}


template<class T, bool Forward>
void Foam::SortableList<T, Forward>::reverseSort()
{
    Foam::sortedOrder(*this, indices_, inverseCompareOp());

    List<T> list(*this, indices_); // Copy with indices for mapping
    List<T>::transfer(list);
}


template<class T, bool Forward>
void Foam::SortableList<T, Forward>::partialSort(label n, label start)
{
    indices_.resize_nocopy(this->size());
    Foam::identity(indices_, 0);

    // Partial sort of the indices
    std::partial_sort
    (
        indices_.begin() + start,
        indices_.begin() + start + n,
        indices_.end(),
        normalCompareOp()
    );

    List<T> list(*this, indices_); // Copy with indices for mapping
    List<T>::transfer(list);
}


template<class T, bool Forward>
void Foam::SortableList<T, Forward>::swap(SortableList<T, Forward>& other)
{
    if (this == &other)
    {
        return;  // Self-swap is a no-op
    }

    List<T>::swap(other);
    indices_.swap(other.indices_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::
operator=(const T& val)
{
    indices_.clear();
    UList<T>::operator=(val);
}


template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::
operator=(const UList<T>& lst)
{
    indices_.clear();
    List<T>::operator=(lst);
}


template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::
operator=(const SortableList<T, Forward>& lst)
{
    if (this == &lst)
    {
        return;  // Self-assigment is a no-op
    }

    List<T>::operator=(lst);
    indices_ = lst.indices();
}


template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::
operator=(List<T>&& lst)
{
    indices_.clear();
    List<T>::transfer(lst);
}


template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::
operator=(SortableList<T, Forward>&& lst)
{
    if (this == &lst)
    {
        return;  // Self-assigment is a no-op
    }

    clear();
    this->swap(lst);
}


template<class T, bool Forward>
inline void Foam::SortableList<T, Forward>::
operator=(std::initializer_list<T> lst)
{
    List<T>::operator=(lst);
    sort();
}


// ************************************************************************* //
