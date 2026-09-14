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

#include "scalarPredicates.H"
#include "HashSet.H"
#include "Tuple2.H"
#include <iterator>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::predicates::scalars::opType
>
Foam::predicates::scalars::opNames
({
    { opType::EQUAL, "eq" },
    { opType::EQUAL, "equal" },
    { opType::NOT_EQUAL,  "neq" },
    { opType::NOT_EQUAL,  "notEqual" },
    { opType::LESS, "lt" },
    { opType::LESS, "less" },
    { opType::LESS_EQUAL, "le" },
    { opType::LESS_EQUAL, "lessEqual" },
    { opType::LESS_EQUAL, "lessEq" },
    { opType::GREATER, "gt" },
    { opType::GREATER, "greater" },
    { opType::GREATER_EQUAL, "ge" },
    { opType::GREATER_EQUAL, "greaterEqual" },
    { opType::GREATER_EQUAL, "greaterEq" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::predicates::scalars::unary
Foam::predicates::scalars::operation(opType op, scalar operand, scalar tol)
{
    switch (op)
    {
        case opType::EQUAL          : return equalOp(operand, tol); break;
        case opType::NOT_EQUAL      : return notEqualOp(operand, tol); break;
        case opType::LESS           : return lessOp(operand); break;
        case opType::LESS_EQUAL     : return lessEqualOp(operand); break;
        case opType::GREATER        : return greaterOp(operand); break;
        case opType::GREATER_EQUAL  : return greaterEqualOp(operand); break;
        case opType::ALWAYS         : return trueOp(); break;
        default : break;
    }

    return falseOp();
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Append list of (opName operand) tuples to target after verifying their
// validity.
template<class ListType, class IteratorType>
void append_impl(ListType& target, IteratorType first, IteratorType last)
{
    const auto len = std::distance(first, last);

    using namespace Foam;
    using std::get;   // for ADL lookups

    labelHashSet badIndices;

    // Check for bad/unknown operations
    {
        label idx = 0;
        for (IteratorType iter = first; iter != last; ++iter)
        {
            const auto& entry = *iter;

            if (!Foam::predicates::scalars::opNames.contains(get<0>(entry)))
            {
                badIndices.insert(idx);
            }
            ++idx;
        }
    }

    // Report bad/unknown operations and exit (Fatal)
    if (!badIndices.empty())
    {
        FatalErrorInFunction
            << "Entries with unknown operations:\n"
            << len << "\n(\n";

        label idx = 0;
        for (IteratorType iter = first; iter != last; ++iter)
        {
            const auto& entry = *iter;

            const bool bad = badIndices.contains(idx);
            ++idx;

            if (bad)
            {
                FatalError << ">>> ";
            }
            else
            {
                FatalError << "    ";
            }
            FatalError << '(' << get<0>(entry) << ' ' << get<1>(entry) << ')';

            if (bad)
            {
                FatalError << " <<<";
            }
            FatalError << '\n';
        }

        FatalError << ")\n"<< exit(FatalError);
    }

    // Appears to be good, fill the list

    const auto start = target.size();

    target.resize(start + len);

    auto output = (target.begin() + start);

    for (IteratorType iter = first; iter != last; ++iter)
    {
        *output = Foam::predicates::scalars::operation(*iter);
        ++output;
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::predicates::scalars::scalars
(
    std::initializer_list<std::pair<const char*, scalar>> entries
)
{
    append_impl(*this, entries.begin(), entries.end());
}


Foam::predicates::scalars::scalars
(
    const UList<Tuple2<word, scalar>>& entries
)
{
    append_impl(*this, entries.begin(), entries.end());
}


Foam::predicates::scalars::scalars(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::predicates::scalars::assign
(
    std::initializer_list<std::pair<const char*, scalar>> entries
)
{
    clear();
    append_impl(*this, entries.begin(), entries.end());
}


void Foam::predicates::scalars::assign
(
    const UList<Tuple2<word, scalar>>& entries
)
{
    clear();
    append_impl(*this, entries.begin(), entries.end());
}


void Foam::predicates::scalars::push_back
(
    opType op,
    scalar operand,
    scalar tolerance
)
{
    List<unary>::push_back
    (
        predicates::scalars::operation(op, operand, tolerance)
    );
}


void Foam::predicates::scalars::push_back
(
    const word& opName,
    scalar operand,
    scalar tolerance
)
{
    // Check for bad/unknown operations
    if (opType op = opNames.lookup(opName, opType::WARN); op != opType::WARN)
    {
        this->push_back(op, operand, tolerance);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown operation: <" << opName << "> with value:"
            << operand
            << exit(FatalError);
    }
}


Foam::label Foam::predicates::scalars::find(scalar value, label pos) const
{
    const label len = this->size();

    if (pos >= 0 && len)
    {
        // auto iter = this->cbegin();

        while (pos < len)
        {
            if ((*this)[pos](value))
            {
                return pos;
            }

            ++pos;
        }
    }

    return -1;
}


Foam::label Foam::predicates::scalars::rfind(scalar value, label pos) const
{
    // pos == -1 has same meaning as std::string::npos - search from end
    if (pos < 0 || pos >= this->size())
    {
        pos = this->size()-1;
    }

    while (pos >= 0)
    {
        if ((*this)[pos](value))
        {
            return pos;
        }

        --pos;
    }

    return -1;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, Foam::predicates::scalars& list)
{
    // Read tuples
    List<Tuple2<word, scalar>> entries(is);

    list.assign(entries);

    return is;
}


// ************************************************************************* //
