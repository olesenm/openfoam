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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned N>
inline constexpr Foam::FixedList<T, N>::FixedList(const T& val)
:
    v_{}
{
    // Raw loop: std::fill() not constexpr until C++20.
    // #if (__cplusplus >= 202002L)
    // std::fill_n(v_, N, val);
    // #endif
    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = val;
    }
}


template<class T, unsigned N>
inline constexpr Foam::FixedList<T, N>::FixedList(Foam::zero)
:
    v_{}
{
    // For scalar types, already value-initialized to zero (by v_{} above)
    // so no further assignment is needed.

    if constexpr (!std::is_scalar_v<T>)
    {
        // Raw loop: std::fill() not constexpr until C++20.
        for (unsigned i = 0; i < N; ++i)
        {
            v_[i] = Foam::zero{};
        }
    }
}


template<class T, unsigned N>
template<class... Args>
inline constexpr Foam::FixedList<T, N>::FixedList
(
    Foam::direct_construct_t,
    Args&&... vals
)
:
    v_{ static_cast<T>(std::forward<Args>(vals))... }
{
    static_assert
    (
        (sizeof...(Args) == N),
        "Number of arguments != number of FixedList elements"
    );
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(std::initializer_list<T> list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList(const UList<T>& list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


// Not constexpr (unknown bounds on the indices)
template<class T, unsigned N>
template<unsigned AnyNum>
inline Foam::FixedList<T, N>::FixedList
(
    const FixedList<T, AnyNum>& list,
    const FixedList<label, N>& indices
)
{
    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = list[indices[i]];
    }
}


template<class T, unsigned N>
inline Foam::FixedList<T, N>::FixedList
(
    const UList<T>& list,
    const FixedList<label, N>& indices
)
{
    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = list[indices[i]];
    }
}


template<class T, unsigned N>
inline Foam::autoPtr<Foam::FixedList<T, N>>
Foam::FixedList<T, N>::clone() const
{
    return autoPtr<FixedList<T, N>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned N>
inline constexpr const T*
Foam::FixedList<T, N>::cdata() const noexcept
{
    return v_;
}


template<class T, unsigned N>
inline constexpr T*
Foam::FixedList<T, N>::data() noexcept
{
    return v_;
}


template<class T, unsigned N>
inline constexpr const char*
Foam::FixedList<T, N>::cdata_bytes() const noexcept
{
    return reinterpret_cast<const char*>(v_);
}


template<class T, unsigned N>
inline constexpr char*
Foam::FixedList<T, N>::data_bytes() noexcept
{
    return reinterpret_cast<char*>(v_);
}


template<class T, unsigned N>
inline std::streamsize Foam::FixedList<T, N>::size_bytes() noexcept
{
    return N*sizeof(T);
}


template<class T, unsigned N>
template<unsigned I>
inline constexpr T& Foam::FixedList<T, N>::get() noexcept
{
    static_assert(I < N, "Address outside FixedList range");
    return v_[I];
}


template<class T, unsigned N>
template<unsigned I>
inline constexpr const T& Foam::FixedList<T, N>::get() const noexcept
{
    static_assert(I < N, "Address outside FixedList range");
    return v_[I];
}


template<class T, unsigned N>
inline constexpr T& Foam::FixedList<T, N>::front() noexcept
{
    return v_[0];
}


template<class T, unsigned N>
inline constexpr const T& Foam::FixedList<T, N>::front() const noexcept
{
    return v_[0];
}


template<class T, unsigned N>
inline constexpr T& Foam::FixedList<T, N>::back() noexcept
{
    return v_[N-1];
}


template<class T, unsigned N>
inline constexpr const T& Foam::FixedList<T, N>::back() const noexcept
{
    return v_[N-1];
}


template<class T, unsigned N>
inline constexpr Foam::label
Foam::FixedList<T, N>::fcIndex(label i) const noexcept
{
    return (i == N-1 ? 0 : i+1);
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::fcValue(label i) const
{
    return this->operator[](this->fcIndex(i));
}


template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::fcValue(label i)
{
    return this->operator[](this->fcIndex(i));
}


template<class T, unsigned N>
inline constexpr Foam::label
Foam::FixedList<T, N>::rcIndex(label i) const noexcept
{
    return (i ? i-1 : N-1);
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::rcValue(label i) const
{
    return this->operator[](this->rcIndex(i));
}


template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::rcValue(label i)
{
    return this->operator[](this->rcIndex(i));
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::checkStart(const label start) const
{
    if (start < 0 || (start && unsigned(start) >= N))
    {
        // Note: always accept start=0, even for zero-sized lists
        FatalErrorInFunction
            << "start " << start << " out of range [0," << N << ")"
            << abort(FatalError);
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::checkSize(const label size) const
{
    if (unsigned(size) != N)
    {
        FatalErrorInFunction
            << "size " << size << " != " << N
            << abort(FatalError);
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::checkIndex(const label i) const
{
    if (i < 0 || unsigned(i) >= N)
    {
        FatalErrorInFunction
            << "index " << i << " out of range [0," << N << ")"
            << abort(FatalError);
    }
}


template<class T, unsigned N>
inline bool Foam::FixedList<T, N>::uniform() const
{
    // This convenience method is not even used for compact output
    // notation in the writeList() method.
    // - requires 'operator!=' to be defined for the type.

    return
    (
        Foam::ListPolicy::uniformity::UNIFORM
     == Foam::ListPolicy::check_uniformity(this->v_, this->v_ + N)
    );
}


template<class T, unsigned N>
inline bool Foam::FixedList<T, N>::contains(const T& val) const
{
    const auto iter = std::find(this->cbegin(), this->cend(), val);
    return (iter != this->cend());
}


template<class T, unsigned N>
inline bool Foam::FixedList<T, N>::contains
(
    const T& val,
    label pos,
    label len
) const
{
    return (this->find(val, pos, len) >= 0);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::resize(const label n)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::resize_fill(const label n, const T& val)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
    this->fill(val);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::resize_nocopy(const label n)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::fill(const T& val)
{
    std::fill_n(v_, N, val);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::fill(Foam::zero)
{
    if constexpr (std::is_scalar_v<T>)
    {
        // Raw loop: std::fill() not constexpr until C++20.
        for (unsigned i = 0; i < N; ++i)
        {
            // Avoid ambiguous casting from Foam::zero
            v_[i] = T(0);
        }
    }
    else if constexpr (is_contiguous_v<T>)
    {
        // Use byte-wise filling for contiguous data
        std::fill_n
        (
            this->data_bytes(), this->size_bytes(), char(0)
        );
    }
    else
    {
        // Raw loop: std::fill() not constexpr until C++20.
        for (unsigned i = 0; i < N; ++i)
        {
            v_[i] = Foam::zero{};
        }
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::swap(FixedList<T, N>& other)
{
    if (this == &other)
    {
        return;  // Self-swap is a no-op
    }

    // Essentially std::swap_ranges
    for (unsigned i=0; i<N; ++i)
    {
        Foam::Swap(v_[i], other.v_[i]);
    }
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::transfer(FixedList<T, N>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // std::execution::unsequenced_policy
    std::move(list.begin(), list.end(), v_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned N>
inline T& Foam::FixedList<T, N>::operator[](label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return v_[i];
}


template<class T, unsigned N>
inline const T& Foam::FixedList<T, N>::operator[](label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return v_[i];
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(const UList<T>& list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(std::initializer_list<T> list)
{
    checkSize(list.size());
    std::copy_n(list.begin(), N, v_);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(const T& val)
{
    this->fill(val);
}


template<class T, unsigned N>
inline void Foam::FixedList<T, N>::operator=(Foam::zero)
{
    this->fill(Foam::zero{});
}


template<class T, unsigned N>
inline constexpr void Foam::FixedList<T, N>::operator=(const FixedList& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // Raw loop: std::copy() not constexpr until C++20.
    // #if (__cplusplus >= 202002L)
    // std::copy(list.begin(), list.end(), v_);
    // #endif
    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = list.v_[i];
    }
}


template<class T, unsigned N>
inline constexpr void Foam::FixedList<T, N>::operator=(FixedList&& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // Raw loop: std::move() not constexpr until C++20.
    // #if (__cplusplus >= 202002L)
    // std::move(list.begin(), list.end(), v_);
    // #endif
    for (unsigned i = 0; i < N; ++i)
    {
        v_[i] = std::move(list.v_[i]);
    }
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::iterator
Foam::FixedList<T, N>::begin() noexcept
{
    return v_;
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::begin() const noexcept
{
    return v_;
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::cbegin() const noexcept
{
    return v_;
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::iterator
Foam::FixedList<T, N>::begin(const int i) noexcept
{
    return (v_ + (i < 0 ? 0 : int(N) < i ? int(N) : i));
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::begin(const int i) const noexcept
{
    return (v_ + (i < 0 ? 0 : int(N) < i ? int(N) : i));
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::cbegin(const int i) const noexcept
{
    return (v_ + (i < 0 ? 0 : int(N) < i ? int(N) : i));
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::iterator
Foam::FixedList<T, N>::end() noexcept
{
    return (v_ + N);
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::end() const noexcept
{
    return (v_ + N);
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_iterator
Foam::FixedList<T, N>::cend() const noexcept
{
    return (v_ + N);
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::reverse_iterator
Foam::FixedList<T, N>::rbegin()
{
    return reverse_iterator(end());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::rbegin() const
{
    return const_reverse_iterator(end());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::crbegin() const
{
    return const_reverse_iterator(end());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::reverse_iterator
Foam::FixedList<T, N>::rend()
{
    return reverse_iterator(begin());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::rend() const
{
    return const_reverse_iterator(begin());
}


template<class T, unsigned N>
inline typename Foam::FixedList<T, N>::const_reverse_iterator
Foam::FixedList<T, N>::crend() const
{
    return const_reverse_iterator(begin());
}


// ************************************************************************* //
