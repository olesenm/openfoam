/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

#include "FixedList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
template<class ListType>
inline void Foam::DynamicList<T, SizeMin>::doAssignDynList
(
    const ListType& list
)
{
    const label len = list.size();

    if (capacity_ < len)
    {
        // Needs more space for the copy operation
        List<T>::setAddressableSize(capacity_);
        List<T>::resize_nocopy(len);
        capacity_ = List<T>::size();
    }

    // Perform copy into addressable portion
    List<T>::setAddressableSize(len);
    List<T>::operator=(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::doCapacity_copy
(
    label count,
    const label newCapacity
)
{
    if (newCapacity == capacity_)
    {
        return;
    }

    // Addressable length, possibly truncated by new capacity
    const label currLen = std::min(List<T>::size(), newCapacity);

    // The count truncated by the new addressable range
    count = std::min(count, currLen);

    // Consistent allocated sizing
    List<T>::setAddressableSize(capacity_);

    if (count > 0)
    {
        List<T>::resize_copy(count, newCapacity);
    }
    else
    {
        List<T>::resize_nocopy(newCapacity);
    }

    capacity_ = List<T>::size();
    List<T>::setAddressableSize(currLen);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::doReserve_copy
(
    label count,
    const label len
)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = List<T>::size();

        // The count truncated by the addressable range
        count = std::min(count, currLen);

        // Consistent allocated sizing
        List<T>::setAddressableSize(capacity_);

        // Increase capacity (eg, doubling)
        // - this may need better handling for when lists become very large
        capacity_ =
            Foam::ListPolicy::reserve_size<SizeMin, 2>(len, capacity_);

        if (count > 0)
        {
            List<T>::resize_copy(count, capacity_);
        }
        else
        {
            List<T>::resize_nocopy(capacity_);
        }

        capacity_ = List<T>::size();
        List<T>::setAddressableSize(currLen);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline constexpr Foam::DynamicList<T, SizeMin>::DynamicList() noexcept
:
    List<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList(label initialCapacity)
:
    List<T>(),
    capacity_(0)
{
    reserve_nocopy(initialCapacity);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    std::pair<label,label> sizing
)
:
    List<T>(std::max(sizing.first, sizing.second)),
    capacity_(List<T>::size())
{
    List<T>::setAddressableSize(sizing.first);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label len,
    const T& val
)
:
    List<T>(len, val),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label len,
    Foam::zero
)
:
    List<T>(len, Foam::zero{}),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, SizeMin>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, AnySizeMin>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const UList<T>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const UList<T>& list,
    const labelUList& indices
)
:
    List<T>(list, indices),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<unsigned N>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const FixedList<T, N>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    std::initializer_list<T> lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<class Addr>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const IndirectListBase<T, Addr>& lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicList<T, SizeMin>&& list
) noexcept
:
    List<T>(std::move(static_cast<List<T>&>(list))),
    capacity_(list.capacity())
{
    list.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicList<T, AnySizeMin>&& list
) noexcept
:
    List<T>(std::move(static_cast<List<T>&>(list))),
    capacity_(list.capacity())
{
    list.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicField<T, AnySizeMin>&& list
) noexcept
:
    List<T>(std::move(static_cast<List<T>&>(list))),
    capacity_(list.capacity())
{
    list.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    List<T>&& list
) noexcept
:
    List<T>(std::move(list)),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    List<T>& content,
    bool reuse
)
:
    List<T>(content, reuse),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicList<T, AnySizeMin>& content,
    bool reuse
)
:
    List<T>(static_cast<List<T>&>(content), reuse),
    capacity_(List<T>::size())
{
    if (reuse)
    {
        capacity_ = content.capacity();
        content.setCapacity_unsafe(0);  // All contents moved
    }
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicField<T, AnySizeMin>& content,
    bool reuse
)
:
    List<T>(static_cast<List<T>&>(content), reuse),
    capacity_(List<T>::size())
{
    if (reuse)
    {
        capacity_ = content.capacity();
        content.setCapacity_unsafe(0);  // All contents moved
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline std::streamsize
Foam::DynamicList<T, SizeMin>::capacity_bytes() const noexcept
{
    return std::streamsize(capacity_)*sizeof(T);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setCapacity
(
    const label len
)
{
    // Preserve all content
    doCapacity_copy(List<T>::size(), len);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setCapacity_nocopy
(
    const label len
)
{
    // Preserve 0 content
    doCapacity_copy(0, len);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::reserve
(
    const label len
)
{
    // Preserve all content
    doReserve_copy(List<T>::size(), len);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::reserve_nocopy
(
    const label len
)
{
    // Preserve 0 content
    doReserve_copy(0, len);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::reserve_exact
(
    const label len
)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = List<T>::size();

        // Consistent allocated sizing
        List<T>::setAddressableSize(capacity_);

        // Exact length
        capacity_ = len;

        // if (!nocopy)
        {
            List<T>::resize_copy(currLen, capacity_);
        }

        capacity_ = List<T>::size();
        List<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize
(
    const label len
)
{
    // Preserve all content
    resize_copy(List<T>::size(), len);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize_copy
(
    label count,
    const label len
)
{
    // Preserve count elements
    doReserve_copy(count, len);
    List<T>::setAddressableSize(len);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize_fill
(
    const label len,
    const T& val
)
{
    resize_nocopy(len);
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize_nocopy
(
    const label len
)
{
    // Preserve 0 content
    resize_copy(0, len);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize
(
    const label len,
    const T& val
)
{
    const label oldLen = List<T>::size();

    // Preserve all content
    resize(len);

    // Fill newly exposed with constant value
    if (oldLen < List<T>::size())
    {
        std::fill
        (
            this->begin(oldLen), this->end(), val
        );
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::clear() noexcept
{
    List<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::clearStorage()
{
    // Extra safety...?
    if (!List<T>::cdata() && List<T>::empty()) capacity_ = 0;

    // Consistent allocated sizing
    List<T>::setAddressableSize(capacity_);
    List<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::shrink_to_fit()
{
    const label currLen = List<T>::size();
    if (currLen < capacity_)
    {
        List<T>::setAddressableSize(capacity_);
        List<T>::resize(currLen);
        capacity_ = List<T>::size();
    }
}


template<class T, int SizeMin>
inline void
Foam::DynamicList<T, SizeMin>::swap(List<T>& other)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&other)
        )
    )
    {
        return;  // Self-swap is a no-op
    }

    // Remove unused storage
    this->shrink_to_fit();

    // Swap storage and addressable size
    UList<T>::swap(other);

    // Update capacity
    capacity_ = List<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::swap
(
    DynamicList<T, AnySizeMin>& other
) noexcept
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&other)
        )
    )
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UList<T>::swap(other);

    // Swap capacity
    const auto oldCapacity = capacity_;
    this->setCapacity_unsafe(other.capacity());
    other.setCapacity_unsafe(oldCapacity);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::swap
(
    DynamicField<T, AnySizeMin>& other
) noexcept
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&other)
        )
    )
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UList<T>::swap(other);

    // Swap capacity
    const auto oldCapacity = capacity_;
    this->setCapacity_unsafe(other.capacity());
    other.setCapacity_unsafe(oldCapacity);
}


template<class T, int SizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer(List<T>& list)
{
    // No check for self-assignment (different types)

    // Consistent allocated sizing
    List<T>::setAddressableSize(capacity_);
    List<T>::transfer(list);
    capacity_ = List<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer
(
    DynamicList<T, AnySizeMin>& other
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&other)
        )
    )
    {
        return;  // Self-assignment is a no-op
    }

    // Consistent allocated sizing
    List<T>::setAddressableSize(capacity_);
    List<T>::transfer(static_cast<List<T>&>(other));

    // Update capacity information
    capacity_ = other.capacity();
    other.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer
(
    DynamicField<T, AnySizeMin>& other
)
{
    // Different types, but check for self-assignment anyhow
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&other)
        )
    )
    {
        return;  // Self-assignment is a no-op
    }

    // Consistent allocated sizing
    List<T>::setAddressableSize(capacity_);
    List<T>::transfer(static_cast<List<T>&>(other));

    // Update capacity information
    capacity_ = other.capacity();
    other.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::DynamicList<T, SizeMin>::emplace_back(Args&&... args)
{
    // This could/should be better with inplace construction
    // (as per std::vector), but currently lacking the methods for that
    // so resize and move assign

    const label idx = List<T>::size();
    resize(idx + 1);

    // move assign element
    UList<T>::operator[](idx) = T(std::forward<Args>(args)...);
    return UList<T>::operator[](idx);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const T& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    UList<T>::operator[](idx) = val;  // copy element
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    T&& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    UList<T>::operator[](idx) = std::move(val);  // move assign element
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const UList<T>& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const UList<T>*>(this)
         == static_cast<const UList<T>*>(&list)
        )
    )
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T, int SizeMin>
template<unsigned N>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const FixedList<T, N>& list
)
{
    const label idx = List<T>::size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    std::initializer_list<T> list
)
{
    const label idx = List<T>::size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T, int SizeMin>
template<class Addr>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const IndirectListBase<T, Addr>& list
)
{
    // Note: push_back will still work even if the indirect list
    // actually references *this, since its source elements will not
    // overlap the new destinations.

    const label idx = this->size();
    const label n = list.size();
    resize(idx + n);

    auto iter = this->begin(idx);

    for (label i = 0; i < n; (void)++i, (void)++iter)
    {
        *iter = list[i];  // copy element
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    List<T>&& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&list)
        )
    )
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    // Move the elements
    std::move(list.begin(), list.end(), this->begin(idx));

    list.clear();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    DynamicList<T, AnySizeMin>&& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&list)
        )
    )
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    // Move the elements
    std::move(list.begin(), list.end(), this->begin(idx));

    list.clearStorage();  // Deletion, capacity=0 etc.
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    DynamicField<T, AnySizeMin>&& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&list)
        )
    )
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    // Move the elements
    std::move(list.begin(), list.end(), this->begin(idx));

    list.clearStorage();  // Deletion, capacity=0 etc.
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::push_uniq(const T& val)
{
    if (this->contains(val))
    {
        return 0;
    }
    else
    {
        this->push_back(val);
        return 1;  // Increased list length by one
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::pop_back(label n)
{
    if (n >= this->size())
    {
        clear();
    }
    else if (n > 0)
    {
        resize(this->size() - n);
    }
}


template<class T, int SizeMin>
inline T Foam::DynamicList<T, SizeMin>::remove()
{
    // Location of last element and simultaneously the new size
    const label idx = List<T>::size() - 1;

    if (idx < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    const T& val = List<T>::operator[](idx);

    List<T>::setAddressableSize(idx);

    return val;
}


template<class T, int SizeMin>
inline T Foam::DynamicList<T, SizeMin>::remove
(
    const label idx,
    const bool fast
)
{
    if (fast)
    {
        // Simply swap idx <=> last
        this->swapLast(idx);
    }
    else
    {
        // Move element to the end and move everything down
        this->moveLast(idx);
    }

    // Element to remove is now at the end
    return this->remove();
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::remove
(
    const labelRange& range
)
{
    return this->removeElements(this->validateRange(range));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline T& Foam::DynamicList<T, SizeMin>::operator()
(
    const label i
)
{
    if (i >= List<T>::size())
    {
        resize(i + 1);
    }

    return UList<T>::operator[](i);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const T& val
)
{
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    Foam::zero
)
{
    UList<T>::operator=(Foam::zero{});
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const UList<T>& lst
)
{
    doAssignDynList(lst);
}


template<class T, int SizeMin>
template<unsigned N>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const FixedList<T, N>& lst
)
{
    doAssignDynList(lst);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const DynamicList<T, SizeMin>& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&list)
        )
    )
    {
        return;  // Self-assignment is a no-op
    }

    doAssignDynList(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const DynamicList<T, AnySizeMin>& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const List<T>*>(this)
         == static_cast<const List<T>*>(&list)
        )
    )
    {
        return;  // Self-assignment is a no-op
    }

    doAssignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    std::initializer_list<T> lst
)
{
    doAssignDynList(lst);
}


template<class T, int SizeMin>
template<class Addr>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const IndirectListBase<T, Addr>& list
)
{
    // NOTE: Self-assignment needs special handling
    // if
    // (
    //     FOAM_UNLIKELY
    //     (
    //         static_cast<const UList<T>*>(this)
    //      == static_cast<const UList<T>*>(&list.values())
    //     )
    // )
    // { ... }

    doAssignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    List<T>&& list
)
{
    this->transfer(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    DynamicList<T, SizeMin>&& list
)
{
    this->transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    DynamicList<T, AnySizeMin>&& list
)
{
    this->transfer(list);
}


// ************************************************************************* //
