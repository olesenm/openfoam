/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

Namespace
    Foam::predicates

Description
    Constant predicate types.

SourceFiles
    predicates.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_predicates_H
#define Foam_predicates_H

#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace predicates
{

/*---------------------------------------------------------------------------*\
                           Class always Declaration
\*---------------------------------------------------------------------------*/

//- Unary, binary, n-ary predicates that \b always returning true,
//- useful for templating.
struct always : std::true_type
{
    //- Unary, binary (or any other n-ary arity >= 1) predicate.
    //  \return true
    template<class T, class... Args>
    constexpr bool operator()(const T&, const Args&...) const noexcept
    {
        return value;
    }

    //- String matching (for example)
    //  \return true
    template<class T, class... Args>
    constexpr bool match(const T&, const Args&...) const noexcept
    {
        return value;
    }
};


/*---------------------------------------------------------------------------*\
                            Class never Declaration
\*---------------------------------------------------------------------------*/

//- Unary, binary, n-ary predicates that \b never return true (only false),
//- useful for templating.
struct never : std::false_type
{
    //- Unary, binary (or any other n-ary arity >= 1) predicate
    //  \return false
    template<class T, class... Args>
    constexpr bool operator()(const T&, const Args&...) const noexcept
    {
        return value;
    }

    //- String matching (for example)
    //  \return false
    template<class T, class... Args>
    constexpr bool match(const T&, const Args&...) const noexcept
    {
        return value;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace predicates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
