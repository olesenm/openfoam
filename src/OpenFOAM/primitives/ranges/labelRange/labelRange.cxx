/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "labelRange.H"
#include "List.H"
#include "MinMax.H"
#include "Pair.H"
#include "IntRangeIO.H"  // For IO operations

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    int labelRange::debug(debug::debugSwitch("labelRange", 0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelRange::labelRange(const MinMax<label>& range) noexcept
:
    labelRange()
{
    if (range.min() <= range.max())
    {
        // max is inclusive, so size with +1. Hope for no overflow
        start() = range.min();
        size()  = 1 + (range.max() - range.min());
    }
}


Foam::labelRange::labelRange(const Pair<label>& start_end) noexcept
:
    labelRange()
{
    if (start_end.first() <= start_end.second())
    {
        // second is exclusive, so size directly. Hope for no overflow
        start() = start_end.first();
        size()  = (start_end.second() - start_end.first());
    }
}


Foam::labelRange::labelRange(Istream& is)
:
    labelRange()
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::labelRange::adjust() noexcept
{
    if (this->start() < 0)
    {
        if (this->size() > 0)
        {
            // Decrease size accordingly
            this->size() += this->start();
        }
        this->start() = 0;
    }
    clampSize();
}


bool Foam::labelRange::overlaps(const labelRange& range, bool touches) const
{
    const label extra = touches ? 1 : 0;

    return
    (
        this->size() && range.size()
     &&
        (
            (
                range.min() >= this->min()
             && range.min() <= this->max() + extra
            )
         ||
            (
                this->min() >= range.min()
             && this->min() <= range.max() + extra
            )
        )
    );
}


Foam::labelRange Foam::labelRange::join(const labelRange& range) const
{
    // Trivial cases first
    if (!this->size())
    {
        return *this;
    }
    else if (!range.size())
    {
        return range;
    }

    const label lower = std::min(this->begin_value(), range.begin_value());
    const label upper = std::max(this->end_value(), range.end_value());

    if (lower < upper)
    {
        return labelRange(lower, (upper-lower));
    }
    else
    {
        return labelRange(lower, 0);
    }
}


Foam::labelRange Foam::labelRange::subset(const labelRange& range) const
{
    const label lower = std::max(this->begin_value(), range.begin_value());
    const label upper = std::min(this->end_value(), range.end_value());

    if (lower < upper)
    {
        return labelRange(lower, (upper-lower));
    }
    else
    {
        return labelRange();
    }
}


Foam::labelRange Foam::labelRange::subset
(
    label start,
    label size
) const
{
    const label lower = Foam::max(this->begin_value(), start);
    const label upper = Foam::min(this->end_value(), start+size);

    if (lower < upper)
    {
        return labelRange(lower, (upper-lower));
    }
    else
    {
        return labelRange();
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::labelRange& range
)
{
    return is >> static_cast<IntRange<label>&>(range);
}


Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::labelRange& range
)
{
    return os << static_cast<const IntRange<label>&>(range);
}


// ************************************************************************* //
