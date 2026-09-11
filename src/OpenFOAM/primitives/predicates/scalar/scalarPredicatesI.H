/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::predicates::scalars::unary
Foam::predicates::scalars::operation
(
    const word& opName,
    scalar operand,
    scalar tolerance
)
{
    return operation(opNames[opName], operand, tolerance);
}


inline Foam::predicates::scalars::unary
Foam::predicates::scalars::operation
(
    const Tuple2<word, scalar>& op,
    scalar tolerance
)
{
    return operation(opNames[op.first()], op.second(), tolerance);
}


inline Foam::predicates::scalars::unary
Foam::predicates::scalars::operation
(
    const std::pair<word, scalar>& op,
    scalar tolerance
)
{
    return operation(opNames[op.first], op.second, tolerance);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::predicates::scalars::contains
(
    scalar value,
    label pos
) const
{
    return (this->find(value, pos) >= 0);
}


inline bool Foam::predicates::scalars::match(scalar value) const
{
    return matchAny(value);
}


inline bool Foam::predicates::scalars::matchAny(scalar value) const
{
    for (const unary& pred : *this)
    {
        if (pred(value))
        {
            return true;
        }
    }

    return false;
}


inline bool Foam::predicates::scalars::matchAll(scalar value) const
{
    for (const unary& pred : *this)
    {
        if (!pred(value))
        {
            return false;
        }
    }

    return (!this->empty());
}


inline Foam::labelList Foam::predicates::scalars::matching
(
    scalar value
) const
{
    labelList indices(this->size());

    label i = 0, count = 0;
    for (const unary& pred : *this)
    {
        if (pred(value))
        {
            indices[count] = i;
            ++count;
        }
        ++i;
    }
    indices.resize(count);

    return indices;
}


inline Foam::labelList Foam::predicates::scalars::matching
(
    const UList<scalar>& input,
    const bool invert
) const
{
    const label len = input.size();

    labelList indices(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (matchAny(input[i]) ? !invert : invert)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


// ************************************************************************* //
