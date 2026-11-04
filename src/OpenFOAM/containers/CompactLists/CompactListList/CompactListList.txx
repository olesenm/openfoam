/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "CompactListList.H"
#include "Istream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
template<class ListListType>
Foam::CompactListList<T> Foam::CompactListList<T>::pack_impl
(
    const ListListType& lists,
    const bool checkOverflow
)
{
    CompactListList<T> compact;

    compact.addr_ = CSRAddressing<label>::pack(lists, checkOverflow);

    auto& newValues = compact.values_;

    const label total = compact.addr_.totalSize();

    if (total)
    {
        // Make a deepCopy of data
        newValues.resize(total);

        auto iter = newValues.begin();

        // NB: operator[] for sub-list read access (eg, an indirect list)
        // cannot replace with std::copy

        for (const auto& list : lists)
        {
            const label sublen = list.size();

            for (label i = 0; i < sublen; (void)++i, (void)++iter)
            {
                *iter = list[i];
            }
        }
    }

    return compact;
}


template<class T>
template<class SubListType>
Foam::CompactListList<T> Foam::CompactListList<T>::pack
(
    const UList<SubListType>& lists,
    const bool checkOverflow
)
{
    return CompactListList<T>::pack_impl<UList<SubListType>>
    (
        lists,
        checkOverflow
    );
}


template<class T>
template<class SubListType, class Addr>
Foam::CompactListList<T> Foam::CompactListList<T>::pack
(
    const IndirectListBase<SubListType, Addr>& lists,
    const bool checkOverflow
)
{
    return CompactListList<T>::pack_impl<IndirectListBase<SubListType, Addr>>
    (
        lists,
        checkOverflow
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::CompactListList<T>::CompactListList(const labelUList& listSizes)
:
    #ifdef FULLDEBUG
    addr_(listSizes, true),
    #else
    addr_(listSizes, false),
    #endif
    values_(addr_.totalSize())
{}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const labelUList& listSizes,
    const T& val
)
:
    #ifdef FULLDEBUG
    addr_(listSizes, true),
    #else
    addr_(listSizes, false),
    #endif
    values_(addr_.totalSize(), val)
{}


template<class T>
Foam::CompactListList<T>::CompactListList(Istream& is)
:
    addr_(is),
    values_(is)
{
    // Optionally: enforceSizeSanity();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
std::streamsize Foam::CompactListList<T>::byteSize() const
{
    if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }
    return this->size_bytes();
}


template<class T>
Foam::List<Foam::labelRange>
Foam::CompactListList<T>::ranges() const
{
    List<labelRange> values(addr_.length());

    forAll(values, i)
    {
        values[i].reset(addr_.localStart(i), addr_.localSize(i));
    }

    return values;
}


template<class T>
void Foam::CompactListList<T>::resize(const labelUList& listSizes)
{
    #ifdef FULLDEBUG
    addr_.reset(listSizes, true);
    #else
    addr_.reset(listSizes, false);
    #endif
    values_.resize(addr_.totalSize());
}


template<class T>
void Foam::CompactListList<T>::resize_nocopy(const labelUList& listSizes)
{
    #ifdef FULLDEBUG
    addr_.reset(listSizes, true);
    #else
    addr_.reset(listSizes, false);
    #endif
    values_.resize_nocopy(addr_.totalSize());
}


template<class T>
void Foam::CompactListList<T>::setLocalSize(label row, label len)
{
    addr_.setLocalSize(row, len);
}


template<class T>
Foam::labelList Foam::CompactListList<T>::localSizes() const
{
    return addr_.localSizes();
}


template<class T>
void Foam::CompactListList<T>::swap
(
    CompactListList<T>& other
)
{
    if (this == &other)
    {
        return;  // Self-swap is a no-op
    }

    addr_.offsets().swap(other.addr_.offsets());
    values_.swap(other.values_);
}


template<class T>
void Foam::CompactListList<T>::transfer
(
    CompactListList<T>& list
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
template<class SubListType, class OutputIter>
OutputIter Foam::CompactListList<T>::copy_unpack
(
    OutputIter d_iter,
    const label pos,
    label len
) const
{
    if (pos >= 0 && pos < this->size())
    {
        // Change sub-length to (one-past) end position
        // len == -1 (like std::string::npos) - search until end

        if (len > 0) len += pos;
        if (len < 0 || len > this->size())
        {
            len = this->size();
        }

        for (label i = pos; i < len; ++i)
        {
            *d_iter = SubListType(this->localList(i));
            ++d_iter;
        }
    }

    return d_iter;
}


template<class T>
template<class SubListType, class OutputIter>
OutputIter Foam::CompactListList<T>::copy_unpack
(
    OutputIter d_iter,
    const labelRange& range
) const
{
    return this->copy_unpack<SubListType>(d_iter, range.start(), range.size());
}


template<class T>
template<class SubListType, class OutputIter>
OutputIter Foam::CompactListList<T>::copy_unpack
(
    OutputIter d_iter,
    const labelUList& indices
) const
{
    for (label i : indices)
    {
        *d_iter = SubListType(this->localList(i));
        ++d_iter;
    }

    return d_iter;
}


// Could also support copy_unpack() with IndirectListBase, as required...
// or the caller can also just use copy_unpack with len = 1 and the
// desired position


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack() const
{
    List<SubListType> lists(size());

    this->copy_unpack<SubListType>(lists.begin());

    return lists;
}


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack(const labelRange& range) const
{
    List<SubListType> lists(range.size());

    this->copy_unpack<SubListType>(lists.begin(), range.start(), range.size());

    return lists;
}


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack(const labelUList& indices) const
{
    List<SubListType> lists(indices.size());

    this->copy_unpack<SubListType>(lists.begin(), indices);

    return lists;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::CompactListList<T>::readList(Istream& is)
{
    addr_.offsets().readList(is);
    values_.readList(is);
    // Optionally: enforceSizeSanity();
    return is;
}


template<class T>
Foam::Ostream& Foam::CompactListList<T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    addr_.offsets().writeList(os, shortLen);
    values_.writeList(os, shortLen);
    return os;
}


template<class T>
Foam::Ostream& Foam::CompactListList<T>::writeMatrix
(
    Ostream& os,
    const label shortLen
) const
{
    const CompactListList<T>& mat = *this;

    const auto oldFmt = os.format(IOstreamOption::ASCII);

    // Write like multi-line matrix of rows (ASCII)
    // TBD: with/without indenting?

    const label nRows = mat.length();

    os  << nRows << nl << token::BEGIN_LIST << nl;

    for (label i = 0; i < nRows; ++i)
    {
        mat.localList(i).writeList(os, shortLen) << nl;
    }

    os  << token::END_LIST << nl;
    os.format(oldFmt);

    return os;
}


// ************************************************************************* //
