/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "IndirectListBase.H"
#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Addr>
Foam::label Foam::IndirectListBase<T, Addr>::find
(
    const T& val,
    label pos,
    label len
) const
{
    if (pos >= 0 && pos < addr_.size())
    {
        // Change sub-length to (one-past) end position
        // len == -1 (like std::string::npos) - search until end

        if (len > 0) len += pos;
        if (len < 0 || len > addr_.size())
        {
            len = addr_.size();
        }

        const T* const vals = values_.begin();

        while (pos < len)
        {
            if (vals[addr_[pos]] == val)
            {
                return pos;
            }

            ++pos;
        }
    }

    return -1;
}


template<class T, class Addr>
Foam::label Foam::IndirectListBase<T, Addr>::rfind
(
    const T& val,
    label pos
) const
{
    // pos == -1 (like std::string::npos) - search from end

    if (pos < 0 || pos >= addr_.size())
    {
        pos = addr_.size()-1;
    }

    const T* const vals = values_.begin();

    while (pos >= 0)
    {
        if (vals[addr_[pos]] == val)
        {
            return pos;
        }

        --pos;
    }

    return -1;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Addr>
Foam::Ostream& Foam::IndirectListBase<T, Addr>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const IndirectListBase<T, Addr>& list = *this;

    const label len = list.size();

    // Unlikely to have IndirectList<char> ...

    constexpr bool with_contiguous
    (
        is_contiguous_v<T> && ListPolicy::allow_contiguous<T>::value
    );
    constexpr bool use_uniformity
    (
        is_contiguous_v<T> && ListPolicy::has_uniformity<T>::value
    );
    constexpr bool no_linebreak
    (
        is_contiguous_v<T> || ListPolicy::no_linebreak<T>::value
    );

    if constexpr (with_contiguous)
    {
        if (os.format() == IOstreamOption::BINARY)
        {
            // Binary and contiguous
            os << nl << len << nl;

            if (len)
            {
                // The TOTAL number of bytes to be written.
                // - possibly add start delimiter
                os.beginRawWrite(len*sizeof(T));

                // Contents
                for (label i=0; i < len; ++i)
                {
                    os.writeRaw
                    (
                        reinterpret_cast<const char*>(&(list[i])),
                        sizeof(T)
                    );
                }

                // End delimiter and/or cleanup.
                os.endRawWrite();
            }

            os.check(FUNCTION_NAME);
            return os;
        }
    }

    if constexpr (use_uniformity)
    {
        if (len > 1 && list.uniform())
        {
            // Two or more entries, all with identical values.
            os << len << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;

            os.check(FUNCTION_NAME);
            return os;
        }
    }

    // Generic output
    // (non-binary or non-contiguous, non-uniform)

    if
    (
        (len <= 1 || !shortLen)
     || (len <= shortLen && no_linebreak)
    )
    {
        // Single-line output

        // Size and start delimiter
        os << len << token::BEGIN_LIST;

        // Contents
        for (label i=0; i < len; ++i)
        {
            if (i) os << token::SPACE;
            os << list[i];
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Multi-line output

        // Size and start delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (label i=0; i < len; ++i)
        {
            os << list[i] << nl;
        }

        // End delimiter
        os << token::END_LIST << nl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
