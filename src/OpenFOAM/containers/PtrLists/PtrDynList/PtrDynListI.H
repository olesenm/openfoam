/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "autoPtr.H"
#include "refPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline constexpr Foam::PtrDynList<T, SizeMin>::PtrDynList() noexcept
:
    PtrList<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList(label initialCapacity)
:
    PtrList<T>(),
    capacity_(0)
{
    reserve(initialCapacity);
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    std::pair<label,label> sizing
)
:
    PtrList<T>(std::max(sizing.first, sizing.second)),
    capacity_(PtrList<T>::size())
{
    PtrList<T>::setAddressableSize(sizing.first);
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    const PtrDynList<T, SizeMin>& list
)
:
    PtrList<T>(list),
    capacity_(PtrList<T>::size())
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrDynList<T, SizeMin>&& list
) noexcept
:
    PtrList<T>(std::move(static_cast<PtrList<T>&>(list))),
    capacity_(list.capacity())
{
    list.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrDynList<T, AnySizeMin>&& list
) noexcept
:
    PtrList<T>(std::move(static_cast<PtrList<T>&>(list))),
    capacity_(list.capacity())
{
    list.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrList<T>&& list
) noexcept
:
    PtrList<T>(std::move(list)),
    capacity_(PtrList<T>::size())
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrList<T>& list,
    bool reuse
)
:
    PtrList<T>(list, reuse),
    capacity_(PtrList<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrDynList<T, AnySizeMin>& list,
    bool reuse
)
:
    PtrList<T>(static_cast<PtrList<T>&>(list), reuse),
    capacity_(PtrList<T>::size())
{
    if (reuse)
    {
        capacity_ = list.capacity();
        list.setCapacity_unsafe(0);  // All contents moved
    }
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList(UList<T*>& list)
:
    PtrList<T>(list),
    capacity_(PtrList<T>::size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::reserve(const label len)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = PtrList<T>::size();

        // Consistent allocated sizing
        PtrList<T>::setAddressableSize(capacity_);

        // Increase capacity (eg, doubling)
        capacity_ =
            Foam::ListPolicy::reserve_size<SizeMin, 2>(len, capacity_);

        // No PtrList<T>::resize_nocopy(...) -> copying nullptr is cheap
        PtrList<T>::resize(capacity_);

        capacity_ = PtrList<T>::size();
        PtrList<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::reserve_exact(const label len)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = PtrList<T>::size();

        // Consistent allocated sizing
        PtrList<T>::setAddressableSize(capacity_);

        // No PtrList<T>::resize_nocopy(...) -> copying nullptr is cheap
        PtrList<T>::resize(len);

        capacity_ = PtrList<T>::size();
        PtrList<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::resize(const label newLen)
{
    auto& ptrs = this->ptrs_;
    const label oldLen = ptrs.size();

    if (capacity_ < newLen)
    {
        // Extend list
        this->reserve(newLen);
    }
    else if (newLen != oldLen)
    {
        // Truncation frees old pointers
        for (label i = newLen; i < oldLen; ++i)
        {
            delete ptrs[i];
            ptrs[i] = nullptr;
        }
    }

    // Adjust addressed size
    PtrList<T>::setAddressableSize(newLen);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::resize_null(const label len)
{
    if (capacity_ < len)
    {
        // Consistent allocated sizing
        PtrList<T>::setAddressableSize(capacity_);

        // Increase capacity (eg, doubling)
        capacity_ =
            Foam::ListPolicy::reserve_size<SizeMin, 2>(len, capacity_);

        PtrList<T>::resize_null(capacity_);
    }
    else
    {
        PtrList<T>::free();  // Free (and nullify) old pointers
    }

    // Adjust addressed size
    PtrList<T>::setAddressableSize(len);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clear()
{
    PtrList<T>::free();  // Free (and nullify) old pointers
    PtrList<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clearStorage()
{
    // Consistent allocated sizing
    PtrList<T>::setAddressableSize(capacity_);
    PtrList<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::shrink_to_fit()
{
    const label currLen = PtrList<T>::size();
    if (currLen < capacity_)
    {
        PtrList<T>::setAddressableSize(capacity_);
        PtrList<T>::resize(currLen);
        capacity_ = PtrList<T>::size();
    }
}


template<class T, int SizeMin>
inline Foam::label Foam::PtrDynList<T, SizeMin>::squeezeNull()
{
    const label newLen = UPtrList<T>::squeezeNull();
    PtrList<T>::setAddressableSize(newLen);
    return newLen;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::swap(PtrList<T>& list)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const PtrList<T>*>(this)
         == static_cast<const PtrList<T>*>(&list)
        )
    )
    {
        return;  // Self-swap is a no-op
    }

    // Remove unused storage
    this->shrink_to_fit();

    // Swap storage and addressable size
    UPtrList<T>::swap(list);

    // Update capacity
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::swap
(
    PtrDynList<T, AnySizeMin>& other
) noexcept
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const PtrList<T>*>(this)
         == static_cast<const PtrList<T>*>(&other)
        )
    )
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UPtrList<T>::swap(other);

    // Swap capacity
    std::swap(this->capacity_, other.capacity_);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::transfer(PtrList<T>& list)
{
    // No check for self-assignment (different types)

    // Consistent allocated sizing
    PtrList<T>::setAddressableSize(capacity_);
    PtrList<T>::transfer(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::transfer
(
    PtrDynList<T, AnySizeMin>& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const PtrList<T>*>(this)
         == static_cast<const PtrList<T>*>(&list)
        )
    )
    {
        return;  // Self assignment is a no-op
    }

    // Consistent allocated sizing
    PtrList<T>::setAddressableSize(capacity_);
    PtrList<T>::transfer(static_cast<PtrList<T>&>(list));

    capacity_ = list.capacity();
    list.setCapacity_unsafe(0);  // All contents moved
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::emplace_back(Args&&... args)
{
    T* ptr = new T(std::forward<Args>(args)...);
    this->push_back(ptr);
    return *ptr;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(T* ptr)
{
    const label idx = this->size();
    resize(idx + 1);
    this->ptrs_[idx] = ptr;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(std::unique_ptr<T>&& ptr)
{
    this->push_back(ptr.release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(autoPtr<T>&& ptr)
{
    this->push_back(ptr.release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(const refPtr<T>& ptr)
{
    this->push_back(ptr.ptr());  // release or clone
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(const tmp<T>& ptr)
{
    this->push_back(ptr.ptr());  // release or clone
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(PtrList<T>&& other)
{
    const label idx = this->size();
    const label len = other.size();

    resize(idx + len);

    for (label i = 0; i < len; ++i)
    {
        set(idx + i, other.release(i));  // Take pointer ownership
    }

    other.clear();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back
(
    PtrDynList<T, AnySizeMin>&& other
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const PtrList<T>*>(this)
         == static_cast<const PtrList<T>*>(&other)
        )
    )
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = this->size();
    const label len = other.size();

    resize(idx + len);

    for (label i = 0; i < len; ++i)
    {
        set(idx + i, other.release(i));  // Take pointer ownership
    }

    other.clearStorage();  // Ensure capacity=0
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::pop_back(label n)
{
    if (n >= this->size())
    {
        this->clear();
    }
    else if (n > 0)
    {
        this->resize(this->size() - n);
    }
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::emplace_set
(
    const label i,
    Args&&... args
)
{
    if (i >= this->size())
    {
        resize(i+1);
    }
    return PtrList<T>::emplace_set(i, std::forward<Args>(args)...);
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::emplace
(
    const label i,
    Args&&... args
)
{
    return this->emplace_set(i, std::forward<Args>(args)...);
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::try_emplace
(
    const label i,
    Args&&... args
)
{
    if (i >= this->size())
    {
        resize(i+1);
    }
    return PtrList<T>::try_emplace(i, std::forward<Args>(args)...);
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    T* ptr
)
{
    if (i >= this->size())
    {
        resize(i+1);
    }

    return autoPtr<T>(UPtrList<T>::set(i, ptr));
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    std::unique_ptr<T>&& ptr
)
{
    return this->set(i, ptr.release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    autoPtr<T>&& ptr
)
{
    return this->set(i, ptr.release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    const refPtr<T>& ptr
)
{
    return this->set(i, ptr.ptr());  // release or clone

}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    const tmp<T>& ptr
)
{
    return this->set(i, ptr.ptr());  // release or clone
}


template<class T, int SizeMin>
Foam::Ostream& Foam::PtrDynList<T, SizeMin>::printAddresses
(
    Ostream& os,
    const bool full
) const
{
    if (full)
    {
        return this->ptrs_.printAddresses(os, capacity_);
    }
    else
    {
        return UPtrList<T>::printAddresses(os);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const UPtrList<T>& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const UPtrList<T>*>(this)
         == static_cast<const UPtrList<T>*>(&list)
        )
    )
    {
        return;  // Self-assignment is a no-op
    }

    this->resize(list.size());
    PtrList<T>::operator=(list);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrDynList<T, SizeMin>& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const UPtrList<T>*>(this)
         == static_cast<const UPtrList<T>*>(&list)
        )
    )
    {
        return;  // Self-assignment is a no-op
    }

    this->resize(list.size());
    PtrList<T>::operator=(static_cast<UPtrList<T>&>(list));
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrDynList<T, AnySizeMin>& list
)
{
    if
    (
        FOAM_UNLIKELY
        (
            static_cast<const PtrList<T>*>(this)
         == static_cast<const PtrList<T>*>(&list)
        )
    )
    {
        return;  // Self-assignment is a no-op
    }

    this->resize(list.size());
    PtrList<T>::operator=(static_cast<UPtrList<T>&>(list));
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrList<T>&& list
)
{
    this->transfer(list);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrDynList<T, SizeMin>&& list
)
{
    this->transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrDynList<T, AnySizeMin>&& list
)
{
    this->transfer(list);
}


// ************************************************************************* //
