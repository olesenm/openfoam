/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

#include "error.H"

// <algorithm> already included by stdFoam.H
#include <execution>

template <typename Func>
void Foam::apply_policy(size_t loop_length, Func&& f)
{
    #ifdef FOAM_OFFLOAD
    if (ListPolicy::use_offload(loop_length))
    {
        f(std::execution::par_unseq);
    }
    else
    #endif
    {
        f(std::execution::seq);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::UList<T>::UList() noexcept
:
    size_(0),
    v_(nullptr)
{}


template<class T>
inline constexpr Foam::UList<T>::UList(T* __restrict__ ptr, label len) noexcept
:
    size_(len),
    v_(ptr)
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::UList<T>::fill_uniform(const T& val)
{
    #ifdef FOAM_OFFLOAD
    if constexpr(is_contiguous_v<T> && !std::is_same_v<T,bool>)
    {
        // For arithmetic types, we can use a more efficient fill
        // operation that avoids the overhead of the generic algorithm.
        apply_policy(this->size_, [&](auto policy)
        {
            std::fill_n(policy, this->v_, this->size_, val);
        });
    }
    else
    #endif
    {
        std::fill_n(this->v_, this->size_, val);
    }
}


template<class T>
inline void Foam::UList<T>::fill_uniform(Foam::zero)
{
    // Avoid ambiguous casting from Foam::zero
    #ifdef FOAM_OFFLOAD
    if constexpr (std::is_arithmetic_v<T>)
    {
        apply_policy(this->size_, [&](auto policy)
        {
            std::fill_n(policy, this->v_, this->size_, T(0));
        });
    }
    else if constexpr (is_contiguous_v<T>)
    {
        constexpr auto val = char(0);

        // Use byte-wise filling for contiguous data
        apply_policy(this->size_bytes(), [&](auto policy)
        {
            std::fill_n(policy, this->data_bytes(), this->size_bytes(), val);
        });
    }
    #else
    if constexpr (std::is_scalar_v<T>)
    {
        // Avoid ambiguous casting from Foam::zero
        std::fill_n(this->v_, this->size_, T(0));
    }
    else if constexpr (is_contiguous_v<T>)
    {
        constexpr auto val = char(0);

        // Use byte-wise filling for contiguous data
        std::fill_n(this->data_bytes(), this->size_bytes(), val);
    }
    #endif
    else
    {
        // May also have special triggers when assigning non-contiguous
        // from zero...

        const auto last = (this->v_ + this->size_);

        for (auto first = this->v_; (first != last); (void)++first)
        {
            *first = Foam::zero{};
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::UList<T>::fcIndex(const label i) const noexcept
{
    return (i == size()-1 ? 0 : i+1);
}


template<class T>
inline Foam::label Foam::UList<T>::rcIndex(const label i) const noexcept
{
    return (i ? i-1 : size()-1);
}


template<class T>
inline const T& Foam::UList<T>::fcValue(const label i) const
{
    return this->operator[](this->fcIndex(i));
}


template<class T>
inline T& Foam::UList<T>::fcValue(const label i)
{
    return this->operator[](this->fcIndex(i));
}


template<class T>
inline const T& Foam::UList<T>::rcValue(const label i) const
{
    return this->operator[](this->rcIndex(i));
}


template<class T>
inline T& Foam::UList<T>::rcValue(const label i)
{
    return this->operator[](this->rcIndex(i));
}


template<class T>
inline void Foam::UList<T>::checkStart(const label start) const
{
    if (start < 0 || (start && start >= size_))
    {
        // Note: accept start=0 for zero-sized lists
        FatalErrorInFunction
            << "start " << start << " out of range [0,"
            << size_ << "]\n"
            << abort(FatalError);
    }
}


template<class T>
inline void Foam::UList<T>::checkSize(const label size) const
{
    if (size < 0 || size > size_)
    {
        FatalErrorInFunction
            << "size " << size << " out of range [0,"
            << size_ << "]\n"
            << abort(FatalError);
    }
}


template<class T>
inline void Foam::UList<T>::checkRange
(
    const label start,
    const label len
) const
{
    // Artificially allow the start of a zero-sized subList to be
    // one past the end of the original list.
    if (len)
    {
        if (len < 0)
        {
            FatalErrorInFunction
                << "size " << len << " is negative, out of range [0,"
                << size_ << "]\n"
                << abort(FatalError);
        }
        this->checkStart(start);
        this->checkSize(start + len);
    }
    else
    {
        // Start index needs to fall between 0 and size.  One position
        // behind the last element is allowed
        this->checkSize(start);
    }
}


template<class T>
inline void Foam::UList<T>::checkIndex(const label i) const
{
    if (!size_)
    {
        FatalErrorInFunction
            << "attempt to access element " << i << " from zero sized list"
            << abort(FatalError);
    }
    else if (i < 0 || i >= size_)
    {
        FatalErrorInFunction
            << "index " << i << " out of range [0,"
            << size_ << "]\n"
            << abort(FatalError);
    }
}


template<class T>
inline bool Foam::UList<T>::uniform() const
{
    // This convenience method is mostly just used for compact output
    // notation in the writeList() method.
    // - requires 'operator!=' to be defined for the type.

    if constexpr (Foam::ListPolicy::has_uniformity<T>::value)
    {
        return
        (
            Foam::ListPolicy::uniformity::UNIFORM
         == Foam::ListPolicy::check_uniformity(this->v_, this->v_ + size_)
        );
    }
    else
    {
        return false;
    }
}


template<class T>
inline T& Foam::UList<T>::front()
{
    return this->operator[](0);
}


template<class T>
inline const T& Foam::UList<T>::front() const
{
    return this->operator[](0);
}


template<class T>
inline T& Foam::UList<T>::back()
{
    return this->operator[](this->size()-1);
}


template<class T>
inline const T& Foam::UList<T>::back() const
{
    return this->operator[](this->size()-1);
}


template<class T>
inline const T* Foam::UList<T>::cdata() const noexcept
{
    return v_;
}


template<class T>
inline T* Foam::UList<T>::data() noexcept
{
    return v_;
}


template<class T>
inline const char* Foam::UList<T>::cdata_bytes() const noexcept
{
    return reinterpret_cast<const char*>(v_);
}


template<class T>
inline char* Foam::UList<T>::data_bytes() noexcept
{
    return reinterpret_cast<char*>(v_);
}


template<class T>
inline std::streamsize Foam::UList<T>::size_bytes() const noexcept
{
    return std::streamsize(size_)*sizeof(T);
}


template<class T>
inline bool Foam::UList<T>::contains(const T& val) const
{
    const auto iter = std::find(this->begin(), this->end(), val);
    return (iter != this->end());
}


template<class T>
inline bool Foam::UList<T>::contains(const T& val, label pos, label len) const
{
    return (this->find(val, pos, len) >= 0);
}


template<class T>
inline void Foam::UList<T>::shallowCopy
(
    T* __restrict__ ptr,
    const label len
) noexcept
{
    size_ = len;
    v_ = ptr;
}


template<class T>
inline void Foam::UList<T>::shallowCopy(std::nullptr_t) noexcept
{
    size_ = 0;
    v_ = nullptr;
}


template<class T>
inline void Foam::UList<T>::shallowCopy(const UList<T>& list) noexcept
{
    size_ = list.size_;
    v_ = list.v_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::UList<T>::operator=(const T& val)
{
    this->fill_uniform(val);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::UList<T>::operator=(Foam::zero)
{
    this->fill_uniform(Foam::zero{});
}


template<class T>
inline T& Foam::UList<T>::operator[](const label i)
{
    if constexpr (std::is_same_v<bool, std::remove_cv_t<T>>)
    {
        // Lazy evaluation - return false for out-of-range
        // Note: strictly speaking should not be modifiable but we cannot
        // alway control which signature (const or non-const) is called
        if (i < 0 || i >= size_)
        {
            return const_cast<bool&>(Foam::pTraits<bool>::null());
        }
    }
    else
    {
        #ifdef FULLDEBUG
        checkIndex(i);
        #endif
    }
    return v_[i];
}


template<class T>
inline const T& Foam::UList<T>::operator[](const label i) const
{
    if constexpr (std::is_same_v<bool, std::remove_cv_t<T>>)
    {
        // Lazy evaluation - return false for out-of-range
        if (i < 0 || i >= size_)
        {
            return Foam::pTraits<bool>::null();
        }
    }
    else
    {
        #ifdef FULLDEBUG
        checkIndex(i);
        #endif
    }
    return v_[i];
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T>
inline typename Foam::UList<T>::iterator
Foam::UList<T>::begin() noexcept
{
    return v_;
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::begin() const noexcept
{
    return v_;
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::cbegin() const noexcept
{
    return v_;
}


template<class T>
inline typename Foam::UList<T>::iterator
Foam::UList<T>::begin(const label i) noexcept
{
    return (v_ + (i < 0 ? 0 : size_ < i ? size_ : i));
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::begin(const label i) const noexcept
{
    return (v_ + (i < 0 ? 0 : size_ < i ? size_ : i));
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::cbegin(const label i) const noexcept
{
    return (v_ + (i < 0 ? 0 : size_ < i ? size_ : i));
}


template<class T>
inline typename Foam::UList<T>::iterator
Foam::UList<T>::end() noexcept
{
    return (v_ + size_);
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::end() const noexcept
{
    return (v_ + size_);
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::cend() const noexcept
{
    return (v_ + size_);
}

template<class T>
inline typename Foam::UList<T>::reverse_iterator
Foam::UList<T>::rbegin()
{
    return reverse_iterator(end());
}

template<class T>
inline typename Foam::UList<T>::const_reverse_iterator
Foam::UList<T>::rbegin() const
{
    return const_reverse_iterator(end());
}

template<class T>
inline typename Foam::UList<T>::const_reverse_iterator
Foam::UList<T>::crbegin() const
{
    return const_reverse_iterator(end());
}

template<class T>
inline typename Foam::UList<T>::reverse_iterator
Foam::UList<T>::rend()
{
    return reverse_iterator(begin());
}

template<class T>
inline typename Foam::UList<T>::const_reverse_iterator
Foam::UList<T>::rend() const
{
    return const_reverse_iterator(begin());
}

template<class T>
inline typename Foam::UList<T>::const_reverse_iterator
Foam::UList<T>::crend() const
{
    return const_reverse_iterator(begin());
}


template<class T>
inline void Foam::UList<T>::setAddressableSize(const label n) noexcept
{
    size_ = n;
}


template<class T>
inline void Foam::UList<T>::swap(UList<T>& list) noexcept
{
    if (&list == this)
    {
        return;  // Self-swap is a no-op
    }

    std::swap(size_, list.size_);
    std::swap(v_, list.v_);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::reverse(UList<T>& list, const label n)
{
    const label nBy2 = n/2;

    for (label i = 0; i < nBy2; ++i)
    {
        Foam::Swap(list[i], list[n-1-i]);
    }
}


template<class T>
inline void Foam::reverse(UList<T>& list)
{
    Foam::reverse(list, list.size());
}


// ************************************************************************* //
