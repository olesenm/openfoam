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

#include "IOList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
bool Foam::IOList<T>::readIOcontents()
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        readStream(typeName) >> *this;
        close();
        return true;
    }

    return false;
}


template<class T>
Foam::label Foam::IOList<T>::readIOsize()
{
    label count(-1);

    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        Istream& is = readStream(typeName);

        token tok(is);

        const bool quick = tok.isLabel();

        if (quick)
        {
            // The majority of files will have lists with sizing info
            count = tok.labelToken();
        }
        is.putBack(tok);

        if (!quick)
        {
            List<T> list(is);
            close();
            count = list.size();
        }
    }

    return count;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::IOList<T>::IOList(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOList<T>>();

    (void) readIOcontents();
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, Foam::zero)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOList<T>>();

    (void) readIOcontents();
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, const label len)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOList<T>>();

    if (!readIOcontents())
    {
        List<T>::resize(len);
    }
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, const label len, const T& val)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOList<T>>();

    if (!readIOcontents())
    {
        List<T>::resize(len, val);
    }
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, const label len, Foam::zero)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOList<T>>();

    if (!readIOcontents())
    {
        List<T>::resize(len, Foam::zero{});
    }
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, const UList<T>& content)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOList<T>>();

    if (!readIOcontents())
    {
        List<T>::operator=(content);
    }
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, List<T>&& content)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOList<T>>();

    // Transfer first to reclaim memory before continuing
    List<T>::transfer(content);

    (void) readIOcontents();
}


template<class T>
Foam::IOListRef<T>::IOListRef
(
    const IOobject& io,
    const UList<T>& content
)
:
    regIOobject(io),
    contentRef_(content)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
Foam::label Foam::IOList<T>::readContentsSize(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }
    rio.resetHeader();

    // Construct NO_READ, changing after construction
    const auto rOpt = rio.readOpt(IOobjectOption::NO_READ);

    IOList<T> reader(rio);
    reader.readOpt(rOpt);

    return reader.readIOsize();
}


template<class T>
Foam::List<T> Foam::IOList<T>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }
    rio.resetHeader();

    IOList<T> reader(rio);

    return List<T>(std::move(static_cast<List<T>&>(reader)));
}


template<class T>
void Foam::IOList<T>::writeContents
(
    const IOobject& io,
    const UList<T>& content
)
{
    IOListRef<T> writer
    (
        IOobject(io, IOobjectOption::NO_REGISTER),
        content
    );

    writer.write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOList<T>::writeData(Ostream& os) const
{
    os << static_cast<const List<T>&>(*this);
    return os.good();
}


template<class T>
bool Foam::IOListRef<T>::writeData(Ostream& os) const
{
    os << contentRef_;
    return os.good();
}


// ************************************************************************* //
