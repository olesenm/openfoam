/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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
inline void Foam::List<T>::doAlloc(const label len)
{
    if (len > 0)
    {
        // With sign-check to avoid spurious -Walloc-size-larger-than
        this->v_ = ListPolicy::allocate<T>(len);
        this->size_ = len;
    }
}


template<class T>
template<class ListType>
inline void Foam::List<T>::copyList(const ListType& list)
{
    // NB: operator[] for list read access (eg, an indirect list)
    // cannot necessarily replace with std::copy
 
    const label len = this->size_;
 
    // For noncontiguous element types or std::string, use sequential loop
    #ifdef FOAM_OFFLOAD
    if constexpr (!is_contiguous_v<T>)
    #endif
    {
        auto iter = this->v_;

        for (label i = 0; i < len; (void)++i, (void)++iter)
        {
            *iter = list[i];
        }
    }
    #ifdef FOAM_OFFLOAD
    else
    {
        T* dst = this->v_;

        apply_policy(len, [&](auto policy)
        {
            std::for_each
            (
                policy,
                dst,
                dst + len,
                [&, base = dst](T& x)
                {
                    const label idx = &x - base;
                    x = list[idx];
                }
            );
        });
    }
    #endif
}


template<class T>
template<class ListType, class ListIndices>
inline void Foam::List<T>::copyList
(
    const ListType& list,
    const ListIndices& indices
)
{
    const label n = this->size_;
    if (n == 0) return;
 
 
    T* dst = this->v_;
 
    // For noncontiguous element types or std::string, use sequential loop
    #ifdef FOAM_OFFLOAD
    if constexpr (!is_contiguous_v<T>)
    #endif
    {
        for (label i = 0; i < n; ++i)
        {
            dst[i] = list[indices[i]];
        }
        return;
    }
    #ifdef FOAM_OFFLOAD
    else
    {
        labelList addr(n);
        std::iota(addr.begin(), addr.end(), label(0));

        apply_policy(n, [&](auto policy)
        {
            std::transform
            (
                policy,
                addr.begin(),
                addr.end(),
                dst,
                [&](const label i)
                {
                    return list[indices[i]];
                }
            );
        });
    }
    #endif
}
 

template<class T>
template<class InputIterator>
inline Foam::List<T>::List
(
    InputIterator input,
    InputIterator inputEnd,  // (unused)
    const label len
)
:
    UList<T>()
{
    if (len > 0)
    {
        doAlloc(len);

        // Like std::copy() or std::copy_n()
        // but without any requirements on the iterator category

        auto iter = this->v_;
        const auto endIter = (iter + len);

        for (; iter != endIter; (void)++iter, (void)++input)
        {
            *iter = *input;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::List<T>::List() noexcept
{}


template<class T>
inline Foam::List<T>::List(List<T>&& list) noexcept
:
    UList<T>(list.data(), list.size())
{
    // Stole content
    list.size_ = 0;
    list.v_ = nullptr;
}


template<class T>
inline Foam::autoPtr<Foam::List<T>> Foam::List<T>::clone() const
{
    return autoPtr<List<T>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::clear()
{
    ListPolicy::deallocate(this->v_, this->size_);
    this->v_ = nullptr;
    this->size_ = 0;
}


namespace Foam
{
    // Template specialization for bool. Fills new entries with false
    template<>
    inline void List<bool>::resize(const label len)
    {
        this->resize(len, false);
    }
}


template<class T>
inline void Foam::List<T>::resize(const label len)
{
    if (this->size_ != len)
    {
        resize_copy(this->size_, len);
    }
}


template<class T>
inline void Foam::List<T>::resize_fill(const label len, const T& val)
{
    resize_nocopy(len);
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::List<T>::resize_nocopy(const label len)
{
    // Same as resize_copy(0, len);
    if (this->size_ != len)
    {
        clear();
        doAlloc(len);
    }
}


template<class T>
inline T& Foam::List<T>::newElmt(const label i)
{
    label n = this->size();

    if (i >= n)
    {
        if (!n) n = 1;  // Avoid dead-lock when starting from zero-sized

        do
        {
            n *= 2;
        }
        while (i >= n);

        resize(n);
    }

    return UList<T>::operator[](i);
}


template<class T>
template<class... Args>
inline T& Foam::List<T>::emplace_back(Args&&... args)
{
    // This could/should be better with inplace construction
    // (as per std::vector), but currently lacking the methods for that
    // so resize and move assign

    const label idx = this->size();
    resize(idx + 1);

    UList<T>::operator[](idx) = T(std::forward<Args>(args)...);
    return UList<T>::operator[](idx);
}


template<class T>
inline void Foam::List<T>::push_back(const T& val)
{
    const label idx = this->size();
    resize(idx + 1);

    UList<T>::operator[](idx) = val;  // copy element
}


template<class T>
inline void Foam::List<T>::push_back(T&& val)
{
    const label idx = this->size();
    resize(idx + 1);

    UList<T>::operator[](idx) = std::move(val);  // move assign element
}


template<class T>
inline void Foam::List<T>::push_back(const UList<T>& list)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted push_back to self" << abort(FatalError);
    }

    const label idx = this->size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T>
template<class Addr>
inline void Foam::List<T>::push_back(const IndirectListBase<T, Addr>& list)
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


template<class T>
inline bool Foam::List<T>::push_uniq(const T& val)
{
    if (!this->contains(val))
    {
        this->push_back(val);
        return true;
    }

    return false;
}


template<class T>
inline void Foam::List<T>::pop_back(label n)
{
    if (n >= this->size())
    {
        this->clear();
    }
    else if (n > 0)
    {
        resize(this->size() - n);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::List<T>::operator=(Foam::zero)
{
    UList<T>::operator=(Foam::zero{});
}


// ************************************************************************* //
