/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "surfZoneIdentifier.H"
#include "dictionary.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZoneIdentifier::surfZoneIdentifier()
:
    index_(0)
{}


Foam::surfZoneIdentifier::surfZoneIdentifier(const label index)
:
    index_(index)
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    const word& name,
    const label index
)
:
    name_(name),
    index_(index)
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    const word& name,
    const label index,
    const word& geometricType
)
:
    name_(name),
    index_(index),
    geometricType_(geometricType)
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    label index,
    const Pair<word>& name_type
)
:
    name_(name_type.first()),
    index_(index),
    geometricType_(name_type.second())
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    surfZoneIdentifier(name, index)
{
    dict.readIfPresent("geometricType", geometricType_);
}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    const surfZoneIdentifier& ident,
    const label newIndex
)
:
    surfZoneIdentifier(ident)
{
    if (newIndex >= 0)
    {
        index_ = newIndex;
    }
}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    surfZoneIdentifier&& ident,
    const label newIndex
)
:
    surfZoneIdentifier(std::move(ident))
{
    if (newIndex >= 0)
    {
        index_ = newIndex;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfZoneIdentifier::write(Ostream& os) const
{
    if (!geometricType_.empty())
    {
        os.writeEntry("geometricType", geometricType_);
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator==
(
    const surfZoneIdentifier& a,
    const surfZoneIdentifier& b
)
{
    return
    (
        (a.index() == b.index())
     && (a.name()  == b.name())
     && (a.geometricType() == b.geometricType())
    );
}


bool Foam::operator!=
(
    const surfZoneIdentifier& a,
    const surfZoneIdentifier& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, surfZoneIdentifier& obj)
{
    // Also read "" for empty words
    obj.name() = word::validate(token(is).stringToken());
    obj.geometricType() = word::validate(token(is).stringToken());

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const surfZoneIdentifier& obj)
{
    // Force unconditional line-breaks on list output.
    // We otherwise risk extremely unreadable entries
    os << nl;

    // Empty words are double-quoted so they are treated as 'string'
    // and not simply lost

    os.writeQuoted(obj.name(), obj.name().empty()) << token::SPACE;
    os.writeQuoted(obj.geometricType(), obj.geometricType().empty());

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
