/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "CompactIOList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
bool Foam::CompactIOList<T>::readIOcontents()
{
    typedef IOList<T> plain_type;

    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        Istream& is = readStream(word::null);

        if (isHeaderClass(typeName))
        {
            // Compact form
            is >> *this;  // or: this->readCompact(is);
        }
        else if (isHeaderClass<plain_type>())
        {
            // Non-compact form
            is >> static_cast<content_type&>(*this);
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Unexpected class name " << headerClassName()
                << " expected " << typeName
                << " or " << plain_type::typeName << endl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }

        close();
        return true;
    }

    return false;
}


template<class T>
Foam::label Foam::CompactIOList<T>::readIOsize()
{
    typedef IOList<T> plain_type;

    label count(-1);

    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        Istream& is = readStream(word::null);

        token tok(is);

        if (tok.isLabel())
        {
            // The majority of files will have lists with sizing prefix
            count = tok.labelToken();

            if (isHeaderClass(typeName))
            {
                // Compact form: read offsets, not content
                if (--count < 0)
                {
                    count = 0;
                }
            }
        }
        else
        {
            is.putBack(tok);

            if (isHeaderClass(typeName))
            {
                // Compact form: can just read the offsets
                labelList offsets(is);
                count = Foam::max(0, (offsets.size()-1));
            }
            else if (isHeaderClass<plain_type>())
            {
                // Non-compact form: need to read everything
                List<T> list(is);
                count = list.size();
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "Unexpected class name " << headerClassName()
                    << " expected " << typeName
                    << " or " << plain_type::typeName << endl
                    << "    while reading object " << name()
                    << exit(FatalIOError);
            }
        }
        close();
    }

    return count;
}


template<class T>
bool Foam::CompactIOList<T>::overflows() const
{
    // Can safely assume that int64 will not overflow
    if constexpr (sizeof(label) < sizeof(int64_t))
    {
        const UList<T>& lists = *this;

        label total = 0;
        for (const auto& list : lists)
        {
            const label prev = total;
            total += list.size();
            if (total < prev)
            {
                return true;
            }
        }
    }
    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::CompactIOList<T>::CompactIOList(const IOobject& io)
:
    regIOobject(io)
{
    (void) readIOcontents();
}


template<class T>
Foam::CompactIOList<T>::CompactIOList
(
    const IOobject& io,
    Foam::zero
)
:
    regIOobject(io)
{
    (void) readIOcontents();
}


template<class T>
Foam::CompactIOList<T>::CompactIOList
(
    const IOobject& io,
    const label len
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        List<T>::resize(len);
    }
}


template<class T>
Foam::CompactIOList<T>::CompactIOList
(
    const IOobject& io,
    const label len,
    const T& val
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        List<T>::resize(len, val);
    }
}


template<class T>
Foam::CompactIOList<T>::CompactIOList
(
    const IOobject& io,
    const label len,
    Foam::zero
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        List<T>::resize(len, Foam::zero{});
    }
}


template<class T>
Foam::CompactIOList<T>::CompactIOList
(
    const IOobject& io,
    const UList<T>& content
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        List<T>::operator=(content);
    }
}


template<class T>
Foam::CompactIOList<T>::CompactIOList
(
    const IOobject& io,
    List<T>&& content
)
:
    regIOobject(io)
{
    // Transfer first to reclaim memory before continuing
    List<T>::transfer(content);

    (void) readIOcontents();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
Foam::label Foam::CompactIOList<T>::readContentsSize(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }
    rio.resetHeader();

    // Construct NO_READ, changing after construction
    const auto rOpt = rio.readOpt(IOobjectOption::NO_READ);

    CompactIOList<T> reader(rio);
    reader.readOpt(rOpt);

    return reader.readIOsize();
}


template<class T>
Foam::List<T> Foam::CompactIOList<T>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }
    rio.resetHeader();

    CompactIOList<T> reader(rio);

    return List<T>(std::move(static_cast<List<T>&>(reader)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::CompactIOList<T>::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if
    (
        streamOpt.format() == IOstreamOption::BINARY
     && overflows()
    )
    {
        streamOpt.format(IOstreamOption::ASCII);

        WarningInFunction
            << "Overall number of elements of CompactIOList (size:"
            << this->size() << ") overflows a label (int"
            << (8*sizeof(label)) << ')' << nl
            << "    Switching to ascii writing" << endl;
    }

    if (streamOpt.format() != IOstreamOption::BINARY)
    {
        // Change type to be non-compact format type
        const word oldTypeName(typeName);

        const_cast<word&>(typeName) = IOList<T>::typeName;

        bool good = regIOobject::writeObject(streamOpt, writeOnProc);

        // Change type back
        const_cast<word&>(typeName) = oldTypeName;

        return good;
    }

    return regIOobject::writeObject(streamOpt, writeOnProc);
}


template<class T>
bool Foam::CompactIOList<T>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::CompactIOList<T>::readCompact(Istream& is)
{
    List<T>& lists = *this;

    // The base type for packed values
    typedef typename T::value_type base_type;

    // Read compact: offsets + packed values
    const labelList offsets(is);
    List<base_type> values(is);

    // Transcribe
    const label len = Foam::max(0, (offsets.size()-1));
    lists.resize_nocopy(len);

    auto iter = values.begin();

    for (label i = 0; i < len; ++i)
    {
        auto& list = lists[i];
        const label count = (offsets[i+1] - offsets[i]);

        list.resize_nocopy(count);

        std::move(iter, iter + count, list.begin());

        iter += count;
    }

    return is;
}


template<class T>
Foam::Ostream& Foam::CompactIOList<T>::writeCompact(Ostream& os) const
{
    const List<T>& lists = *this;

    // The base type for packed values
    typedef typename T::value_type base_type;

    // Convert to compact format
    label total = 0;
    const label len = lists.size();

    // offsets
    {
        labelList offsets(len+1);

        for (label i = 0; i < len; ++i)
        {
            offsets[i] = total;
            total += lists[i].size();

            if (total < offsets[i])
            {
                FatalIOErrorInFunction(os)
                    << "Overall number of elements of CompactIOList (size:"
                    << len
                    << ") overflows the representation of a label" << nl
                    << "Please recompile with a larger representation"
                    << " for label" << exit(FatalIOError);
            }
        }
        offsets[len] = total;
        os << offsets;
    }

    // packed values: make deepCopy for writing
    {
        List<base_type> values(total);

        auto iter = values.begin();

        for (const auto& list : lists)
        {
            iter = std::copy_n(list.begin(), list.size(), iter);

            // With IndirectList? [unlikely]
            // const label count = list.size();
            // for (label i = 0; i < count; (void)++i, (void)++iter)
            // {
            //     *iter = list[i];
            // }
        }
        os << values;
    }

    return os;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::CompactIOList<T>& lists
)
{
    return lists.readCompact(is);
}


template<class T>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::CompactIOList<T>& lists
)
{
    // Keep ASCII writing same
    if (os.format() != IOstreamOption::BINARY)
    {
        os << static_cast<const List<T>&>(lists);
    }
    else
    {
        lists.writeCompact(os);
    }

    return os;
}


// ************************************************************************* //
