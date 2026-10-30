/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::CompactListList<T>::enforceSizeSanity()
{
    // CSRAdressing should avoid size=1, but be certain
    if (addr_.offsets().size() == 1)
    {
        addr_.offsets().clear();
    }
    if (addr_.offsets().empty())
    {
        values_.clear();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const CompactListList<T>& list
)
:
    addr_(list.addr_),
    values_(list.values_)
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    CompactListList<T>&& list
)
:
    addr_(std::move(list.addr_)),
    values_(std::move(list.values_))
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    CompactListList<T>& list,
    bool reuse
)
:
    addr_(list.addr_, reuse),
    values_(list.values_, reuse)
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label mRows,
    const label nVals
)
:
    addr_(mRows),
    values_(nVals)
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label mRows,
    const label nVals,
    Foam::zero
)
:
    addr_(mRows),
    values_(nVals, Foam::zero{})
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label mRows,
    const label nVals,
    const T& val
)
:
    addr_(mRows),
    values_(nVals, val)
{}


template<class T>
inline Foam::autoPtr<Foam::CompactListList<T>>
Foam::CompactListList<T>::clone() const
{
    return autoPtr<CompactListList<T>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T* Foam::CompactListList<T>::cdata() const noexcept
{
    return values_.cdata();
}


template<class T>
inline T* Foam::CompactListList<T>::data() noexcept
{
    return values_.data();
}


template<class T>
inline const char* Foam::CompactListList<T>::cdata_bytes() const noexcept
{
    return values_.cdata_bytes();
}


template<class T>
inline char* Foam::CompactListList<T>::data_bytes() noexcept
{
    return values_.data_bytes();
}


template<class T>
inline std::streamsize Foam::CompactListList<T>::size_bytes() const noexcept
{
    return values_.size_bytes();
}


template<class T>
inline bool Foam::CompactListList<T>::empty() const noexcept
{
    // Not using addr_.empty() since that has extra logic
    return (addr_.length() == 0);
}


template<class T>
inline bool Foam::CompactListList<T>::single() const noexcept
{
    return addr_.single();
}


template<class T>
inline Foam::label Foam::CompactListList<T>::length() const noexcept
{
    // Same as size()
    return addr_.length();
}


template<class T>
inline Foam::label Foam::CompactListList<T>::size() const noexcept
{
    // Same as length()
    return addr_.size();
}


template<class T>
inline Foam::labelList Foam::CompactListList<T>::sizes() const
{
    return localSizes();
}


template<class T>
inline Foam::label Foam::CompactListList<T>::totalSize() const noexcept
{
    return addr_.totalSize();
}


template<class T>
inline Foam::label Foam::CompactListList<T>::maxSize() const noexcept
{
    return addr_.maxSize();
}


template<class T>
inline const Foam::labelUList
Foam::CompactListList<T>::localStarts() const
{
    return addr_.localStarts();
}


template<class T>
inline Foam::label Foam::CompactListList<T>::localStart(label i) const
{
    return addr_.localStart(i);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::localEnd(label i) const
{
    return addr_.localEnd(i);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::localSize(label i) const
{
    return addr_.localSize(i);
}


template<class T>
inline Foam::labelRange Foam::CompactListList<T>::range(label i) const noexcept
{
    return addr_.range(i);
}


template<class T>
inline const Foam::SubList<T>
Foam::CompactListList<T>::localList(label i) const
{
    return values_.slice(addr_.localStart(i), addr_.localSize(i));
}


template<class T>
inline Foam::SubList<T>
Foam::CompactListList<T>::localList(label i)
{
    return values_.slice(addr_.localStart(i), addr_.localSize(i));
}


template<class T>
inline Foam::label Foam::CompactListList<T>::toGlobal(label row, label i) const
{
    return addr_.toGlobal(row, i);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::toLocal(label row, label i) const
{
    return addr_.toLocal(row, i);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::findRow(label i) const
{
    return addr_.findRow(i);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::whichRow(label i) const
{
    label foundRow = addr_.findRow(i);

    if (foundRow  < 0)
    {
        FatalErrorInFunction
            << "Index " << i << " outside of range" << nl
            << "Offsets:" << addr_.offsets()
            << abort(FatalError);
    }

    return foundRow;
}


template<class T>
inline void Foam::CompactListList<T>::clear()
{
    addr_.clear();
    values_.clear();
}


template<class T>
inline void Foam::CompactListList<T>::resize(const label mRows)
{
    if (mRows < 1)
    {
        // Clear
        addr_.clear();
        values_.clear();
    }
    else if (mRows < size())
    {
        // Shrink
        addr_.resize(mRows);
        values_.resize(addr_.totalSize());
    }
    else if (mRows > size())
    {
        // Extend number of rows, each with local size of 0
        addr_.resize(mRows);
    }
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label mRows,
    const label nVals
)
{
    if (mRows < 1)
    {
        // Enforce sizing sanity
        addr_.clear();
        values_.clear();
    }
    else
    {
        addr_.offsets().resize(mRows+1, 0);
        values_.resize(nVals);
    }
}


template<class T>
inline void Foam::CompactListList<T>::resize_nocopy
(
    const label mRows,
    const label nVals
)
{
    if (mRows < 1)
    {
        // Enforce sizing sanity
        addr_.clear();
        values_.clear();
    }
    else
    {
        addr_.offsets().resize(mRows+1, 0);
        values_.resize_nocopy(nVals);
    }
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label mRows,
    const label nVals,
    const T& val
)
{
    if (mRows < 1)
    {
        // Enforce sizing sanity
        addr_.clear();
        values_.clear();
    }
    else
    {
        addr_.offsets().resize(mRows+1, 0);
        values_.resize(nVals, val);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::CompactListList<T>::operator=
(
    const CompactListList<T>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    addr_ = list.addr_;
    values_ = list.values_;
}


template<class T>
inline void Foam::CompactListList<T>::operator=
(
    CompactListList<T>&& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    addr_.offsets().transfer(list.addr_.offsets());
    values_.transfer(list.values_);
}


template<class T>
inline void Foam::CompactListList<T>::operator=(const T& val)
{
    values_ = val;
}


template<class T>
inline void Foam::CompactListList<T>::operator=(Foam::zero)
{
    values_ = Foam::zero{};
}


template<class T>
inline const Foam::SubList<T>
Foam::CompactListList<T>::operator[](label i) const
{
    return this->localList(i);
}


template<class T>
inline Foam::SubList<T>
Foam::CompactListList<T>::operator[](label i)
{
    return this->localList(i);
}


template<class T>
inline T& Foam::CompactListList<T>::operator()(label i, label j)
{
    return values_[toGlobal(i, j)];
}


template<class T>
inline const T& Foam::CompactListList<T>::operator()(label i, label j) const
{
    return values_[toGlobal(i, j)];
}


// ************************************************************************* //
