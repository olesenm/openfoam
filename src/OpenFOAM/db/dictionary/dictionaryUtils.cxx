/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Mark Olesen
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

#include "dictionary.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Simple type codes to sort streams first (delaying recursion).
// At the moment we only have two known categories (stream/dictionary)
// - no idea what the 'unknown' fallback would be, but provision for it
//   anyhow
//
enum class typeCode : uint8_t { stream = 0, dict, unknown };

using tuple_type = std::pair<const Foam::entry*, typeCode>;

// Extract entries as (pointer, typeCode) pairs and sort by
// (typeCode, keyword) so that comparisons are order-independent and
// primitiveEntry types are checked before dictionaryEntry types
// (delaying recursion).
Foam::List<tuple_type> extract_sortedEntries(const Foam::dictionary& dict)
{
    using namespace Foam;

    const label nEntries = dict.size();

    List<tuple_type> ptrs(nEntries);

    auto iter = dict.begin();
    for (label i = 0; i < nEntries; ++i, ++iter)
    {
        const entry& e = *iter;

        typeCode code =
        (
            e.isStream() ? typeCode::stream
          : e.isDict()   ? typeCode::dict
          : typeCode::unknown
        );

        ptrs[i] = tuple_type(&e, code);
    }

    std::sort
    (
        ptrs.begin(), ptrs.end(),
        [](const tuple_type& a, const tuple_type& b)
        {
            if (a.second != b.second)
            {
                // primitive before dictionary entry
                return a.second < b.second;
            }

            // sort by name (ignore regex vs non-regex)
            // - this is the same as what sortedToc() does.
            return
            (
                static_cast<const word&>(a.first->keyword())
              < static_cast<const word&>(b.first->keyword())
            );
        }
    );

    return ptrs;
}


// Compare two pre-sorted (pointer, typeCode) lists (equal size assumed)
// for identical types and keywords, without looking at the content.
bool test_sameKeys
(
    const Foam::UList<tuple_type>& a,
    const Foam::UList<tuple_type>& b
)
{
    using namespace Foam;

    if (a.size() != b.size())
    {
        // Sizes must match
        return false;
    }

    const label nEntries = a.size();

    // First sweep - require identical types (cheap)
    for (label i = 0; i < nEntries; ++i)
    {
        if (a[i].second != b[i].second)
        {
            return false;
        }
    }

    // Second sweep - require identical names
    for (label i = 0; i < nEntries; ++i)
    {
        if
        (
            static_cast<const word&>(a[i].first->keyword())
         != static_cast<const word&>(b[i].first->keyword())
        )
        {
            return false;
        }
    }

    return true;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::dictionary::depth(const dictionary& dict)
{
    label maxDepth = 0;

    for (const entry& e : dict)
    {
        if (const auto* dictptr = e.dictPtr(); dictptr)
        {
            label subDepth = 1 + depth(*dictptr);

            if (maxDepth < subDepth)
            {
                maxDepth = subDepth;
            }
        }
    }

    return maxDepth;
}


bool Foam::dictionary::equal_keys(const dictionary& a, const dictionary& b)
{
    if (a.size() != b.size())
    {
        // Sizes must match
        return false;
    }
    else if (a.empty())
    {
        // Both empty
        return true;
    }

    const List<tuple_type> ptrs1(extract_sortedEntries(a));
    const List<tuple_type> ptrs2(extract_sortedEntries(b));

    return test_sameKeys(ptrs1, ptrs2);
}


bool Foam::dictionary::equal_contents(const dictionary& a, const dictionary& b)
{
    if (a.size() != b.size())
    {
        // Sizes must match
        return false;
    }
    else if (a.empty())
    {
        // Both empty
        return true;
    }

    const List<tuple_type> ptrs1(extract_sortedEntries(a));
    const List<tuple_type> ptrs2(extract_sortedEntries(b));

    if (!test_sameKeys(ptrs1, ptrs2))
    {
        return false;
    }

    // - The sorted extraction (above) ensures that primitive entries will
    //   be tested first and dictionaries (with recursion) checked later.
    // - After test_sameKeys() it is known that types are in lockstep

    const auto nEntries = ptrs1.size();

    // Compare content per-entry
    for (label i = 0; i < nEntries; ++i)
    {
        const auto& entry1 = *(ptrs1[i].first);
        const auto& entry2 = *(ptrs2[i].first);

        switch (ptrs1[i].second)  // Both have same typeCode (already checked)
        {
            case typeCode::stream :
            {
                if
                (
                    static_cast<const tokenList&>(entry1.stream())
                 != static_cast<const tokenList&>(entry2.stream())
                )
                {
                    return false;
                }
                break;
            }

            case typeCode::dict :
            {
                if (!equal_contents(entry1.dict(), entry2.dict()))
                {
                    return false;
                }
                break;
            }

            default :
            {
                break;
            }
        }
    }

    // Everything compared
    return true;
}


// ************************************************************************* //
