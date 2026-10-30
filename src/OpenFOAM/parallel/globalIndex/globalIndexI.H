/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
    Copyright (C) 2026 Mark Olesen
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

template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    const UList<IntType>& listOffsets
)
:
    globalIndexBase(),
    addr_type(Foam::direct_construct_t{}, listOffsets)
{}


template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    List<IntType>&& listOffsets
)
:
    globalIndexBase(),
    addr_type(Foam::direct_construct_t{}, std::move(listOffsets))
{}


template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    globalIndexBase::gatherNone,
    const labelUList& localSizes
)
:
    globalIndexBase(),
    addr_type()
{
    reset(localSizes);
}


template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    label localSize,
    int communicator,
    bool parallel
)
:
    globalIndexBase(),
    addr_type()
{
    reset(localSize, communicator, parallel);
}


template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    globalIndexBase::gatherNone,
    label localSize,
    int communicator  /* unused: (not applicable) */
)
:
    globalIndexBase(),
    addr_type()
{
    addr_type::single(localSize);
}


template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    globalIndexBase::gatherOnly,
    label localSize,
    int communicator,
    bool parallel
)
:
    globalIndexBase(),
    addr_type()
{
    reset(globalIndexBase::gatherOnly{}, localSize, communicator, parallel);
}


template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    globalIndexBase::gatherOffset,
    label localSize,
    int communicator,
    bool parallel
)
:
    globalIndexBase(),
    addr_type()
{
    reset(globalIndexBase::gatherOffset{}, localSize, communicator, parallel);
}


template<class IntType>
inline Foam::GlobalIndex<IntType>::GlobalIndex
(
    globalIndexBase::gatherTotal,
    label localSize,
    int communicator,
    bool parallel
)
:
    globalIndexBase(),
    addr_type()
{
    reset(globalIndexBase::gatherTotal{}, localSize, communicator, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
inline Foam::OffsetRange<IntType>
Foam::GlobalIndex<IntType>::slice(label proci) const noexcept
{
    if
    (
        const auto len = (offsets().size() - 1);
        (len > 0 && proci >= 0 && proci < len)
    )
    {
        const auto* off = offsets().cdata();
        return { off[proci], (off[proci+1] - off[proci]), off[len] };
    }
    else
    {
        return {};
    }
}


template<class IntType>
inline void Foam::GlobalIndex<IntType>::reset
(
    globalIndexBase::gatherNone,
    label localSize,
    int communicator  /* unused: (not applicable) */
)
{
    addr_type::single(localSize);
}


template<class IntType>
inline Foam::label Foam::GlobalIndex<IntType>::whichProcID
(
    label proci,
    IntType i
) const
{
    label foundProc = addr_type::findRow(proci, i);

    if (foundProc < 0)
    {
        FatalErrorInFunction
            << "Index:" << i
            << " does not belong on any processor."
            << "\nOffsets:" << addr_type::offsets()
            << Foam::abort(FatalError);
    }

    return foundProc;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::labelRange Foam::globalIndex::allProcs() const noexcept
{
    // Proc 0 -> nProcs
    const label len = (offsets().size() - 1);
    return (len < 1) ? labelRange() : labelRange(0, len);
}


inline Foam::labelRange Foam::globalIndex::subProcs() const noexcept
{
    // Proc 1 -> nProcs
    const label len = (offsets().size() - 2);
    return (len < 1) ? labelRange() : labelRange(1, len);
}


inline Foam::labelList Foam::globalIndex::toGlobal
(
    label proci,
    const labelUList& labels
) const
{
    // Or using std::transform

    //std::transform(labels.begin(), labels.end(), result.begin(),
    //    [=](label id) { return id += off });

    labelList result(labels);
    inplaceToGlobal(proci, result);

    return result;
}


inline void Foam::globalIndex::inplaceToGlobal
(
    label proci,
    labelUList& labels
) const
{
    if (const auto beg = localStart(proci); beg)
    {
        for (auto& val : labels)
        {
            val += beg;
        }
    }
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::globalIndex::const_iterator::
const_iterator(const globalIndex* gi, label i) noexcept
:
    parent_(gi),
    index_(i)
{}


inline Foam::label Foam::globalIndex::const_iterator::
index() const noexcept
{
    return index_;
}


inline Foam::label Foam::globalIndex::const_iterator::
start() const
{
    return (*parent_).localStart(index_);
}


inline Foam::label Foam::globalIndex::const_iterator::
size() const
{
    return (*parent_).localSize(index_);
}


inline Foam::labelRange Foam::globalIndex::const_iterator::
range() const
{
    return (*parent_).range(index_);
}


inline Foam::labelRange Foam::globalIndex::const_iterator::
operator*() const
{
    return this->range();
}


inline Foam::globalIndex::const_iterator&
Foam::globalIndex::const_iterator::
operator++() noexcept
{
    ++index_;
    return *this;
}


inline Foam::globalIndex::const_iterator
Foam::globalIndex::const_iterator::
operator++(int) noexcept
{
    const_iterator old(*this);
    ++index_;
    return old;
}


inline Foam::globalIndex::const_iterator&
Foam::globalIndex::const_iterator::
operator--() noexcept
{
    --index_;
    return *this;
}


inline Foam::globalIndex::const_iterator
Foam::globalIndex::const_iterator::
operator--(int) noexcept
{
    const_iterator old(*this);
    --index_;
    return old;
}


inline bool
Foam::globalIndex::const_iterator::
operator==
(
    const const_iterator& iter
) const noexcept
{
    return (index_ == iter.index_);
}


inline bool
Foam::globalIndex::const_iterator::
operator!=
(
    const const_iterator& iter
) const noexcept
{
    return (index_ != iter.index_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::globalIndex::const_iterator
Foam::globalIndex::cbegin() const noexcept
{
    return const_iterator(this);
}


inline const Foam::globalIndex::const_iterator
Foam::globalIndex::cend() const noexcept
{
    return const_iterator(this, this->length());
}


inline Foam::globalIndex::const_iterator
Foam::globalIndex::begin() const noexcept
{
    return const_iterator(this);
}


inline const Foam::globalIndex::const_iterator
Foam::globalIndex::end() const noexcept
{
    return const_iterator(this, this->length());
}


inline Foam::globalIndex::const_iterator
Foam::globalIndex::cbegin(label i) const noexcept
{
    const label len = this->length();
    return const_iterator(this, (i < 0 ? 0 : len < i ? len : i));
}


inline Foam::globalIndex::const_iterator
Foam::globalIndex::begin(label i) const noexcept
{
    return this->cbegin(i);
}


// ************************************************************************* //
