/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "IOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::IOField<Type>::readIOcontents(bool readOnProc)
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // Do reading
        Istream& is = readStream(typeName, readOnProc);

        if (readOnProc)
        {
            is >> *this;
        }
        close();
        return true;
    }

    return false;
}


template<class T>
Foam::label Foam::IOField<T>::readIOsize()
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    (void) readIOcontents();
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const bool readOnProc)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    (void) readIOcontents(readOnProc);
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, Foam::zero)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    (void) readIOcontents();
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const label len)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if (!readIOcontents())
    {
        Field<Type>::resize(len);
    }
}


template<class Type>
Foam::IOField<Type>::IOField
(
    const IOobject& io,
    const label len,
    const Type& val
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if (!readIOcontents())
    {
        Field<Type>::resize(len, val);
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const label len, Foam::zero)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if (!readIOcontents())
    {
        Field<Type>::resize(len, Foam::zero{});
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const UList<Type>& content)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if (!readIOcontents())
    {
        Field<Type>::operator=(content);
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, Field<Type>&& content)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    // Transfer first to reclaim memory before continuing
    Field<Type>::transfer(content);

    readIOcontents();
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const tmp<Field<Type>>& tfld)
:
    regIOobject(io)
{
    const bool reuse = tfld.movable();

    if (reuse)
    {
        Field<Type>::transfer(tfld.ref());
    }

    if (!readIOcontents() && !reuse)
    {
        Field<Type>::operator=(tfld());
    }

    tfld.clear();
}


template<class Type>
Foam::IOFieldRef<Type>::IOFieldRef
(
    const IOobject& io,
    const Field<Type>& content
)
:
    regIOobject(io),
    contentRef_(content)  // cref
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::IOField<Type>::readContentsSize(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }
    rio.resetHeader();

    // Construct NO_READ, changing after construction
    const auto rOpt = rio.readOpt(IOobjectOption::NO_READ);

    IOField<Type> reader(rio);
    reader.readOpt(rOpt);

    return reader.readIOsize();
}


template<class Type>
Foam::Field<Type> Foam::IOField<Type>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }
    rio.resetHeader();

    IOField<Type> reader(rio);

    return Field<Type>(std::move(static_cast<Field<Type>&>(reader)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOField<Type>::writeData(Ostream& os) const
{
    os << static_cast<const Field<Type>&>(*this);
    return os.good();
}


template<class Type>
bool Foam::IOFieldRef<Type>::writeData(Ostream& os) const
{
    os << contentRef_.cref();
    return os.good();
}


// ************************************************************************* //
