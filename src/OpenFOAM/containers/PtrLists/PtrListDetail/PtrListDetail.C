/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "PtrListDetail.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::label Foam::Detail::PtrListDetail<T>::count_nonnull() const noexcept
{
    label n = 0;

    for (const T* ptr : *this)
    {
        if (ptr)
        {
            ++n;
        }
    }

    return n;
}


template<class T>
Foam::label Foam::Detail::PtrListDetail<T>::find_first(label pos)
const noexcept
{
    const label len = this->size();
    const auto ptrs = this->cdata();

    if (pos < 0) pos = 0;

    for (; pos < len; ++pos)
    {
        if (ptrs[pos])
        {
            return pos;
        }
    }

    return -1;
}


template<class T>
Foam::label Foam::Detail::PtrListDetail<T>::find_first_not(label pos)
const noexcept
{
    const label len = this->size();
    const auto ptrs = this->cdata();

    if (pos < 0) pos = 0;

    for (; pos < len; ++pos)
    {
        if (!ptrs[pos])
        {
            return pos;
        }
    }

    return -1;
}


template<class T>
Foam::label Foam::Detail::PtrListDetail<T>::find_next(label pos) const
noexcept
{
    return this->find_first(pos+1);
}


template<class T>
Foam::label Foam::Detail::PtrListDetail<T>::find_next_not(label pos) const
noexcept
{
    return this->find_first_not(pos+1);
}


template<class T>
void Foam::Detail::PtrListDetail<T>::free()
{
    // Presume they were allocated from front to back...
    List<T*>& ptrs = *this;

    for (auto i = this->size()-1; i >= 0; --i)
    {
        delete ptrs[i];
        ptrs[i] = nullptr;
    }
}


template<class T>
template<class... Args>
Foam::Detail::PtrListDetail<T>
Foam::Detail::PtrListDetail<T>::clone(Args&&... args) const
{
    const List<T*>& ptrs = *this;
    const label len = ptrs.size();

    PtrListDetail<T> cloned(len);

    for (label i = 0; i < len; ++i)
    {
        if (const T* ptr = ptrs[i]; ptr)
        {
            cloned[i] = ptr->clone(std::forward<Args>(args)...).ptr();
        }
    }

    return cloned;
}


// ************************************************************************* //
