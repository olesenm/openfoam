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

Class
    Foam::predicates::scalars

Description
    A list of unary predicates (tests) on scalars.
    Includes a number of standard comparison predicates
    (eg, "less", "greater", ...)

Note
    This class is still subject to larger changes (2018-11) as it matures.

SeeAlso
    Foam::scalarRange, Foam::scalarRanges

SourceFiles
    scalarPredicates.cxx
    scalarPredicatesI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_scalarPredicates_H
#define Foam_scalarPredicates_H

#include "scalar.H"
#include "predicates.H"
#include "Enum.H"
#include "List.H"
#include "Tuple2.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace predicates
{

/*---------------------------------------------------------------------------*\
                           Class scalars Declaration
\*---------------------------------------------------------------------------*/

class scalars
:
    public List<std::function<bool(Foam::scalar)>>
{
public:

    // Public Types

        //- The unary function type for testing a scalar
        using unary = std::function<bool(Foam::scalar)>;

        //- The default tolerance to use is \c VSMALL
        static constexpr scalar defaultTol() noexcept { return VSMALL; }

        //- Enumerations for some standard comparison predicates
        enum class opType
        {
            EQUAL,              //!< "eq", "equal"
            NOT_EQUAL,          //!< "neq", "notEqual"
            LESS,               //!< "lt", "less"
            LESS_EQUAL,         //!< "le", "lessEqual", "lessEq"
            GREATER,            //!< "gt", "greater"
            GREATER_EQUAL,      //!< "ge", "greaterEqual", "greaterEq"
            ALWAYS,             //!< Always matches
            NEVER,              //!< Never matches
            WARN,               //!< A warning operation (false)
        };

        //- Names for the opType enumeration.
        static const Enum<opType> opNames;


    // Standard comparison methods

        //- Compare for equality to \p target,
        //- within specified \p tolerance (non-negative)
        static unary equalOp(scalar target, scalar tolerance = defaultTol())
        {
            return [target,tolerance](scalar x)
            {
                return (std::fabs(target - x) <= tolerance);
            };
        }

        //- Compare for inequality to \p target,
        //- within specified \p tolerance (non-negative)
        static unary notEqualOp(scalar target, scalar tolerance = defaultTol())
        {
            return [target,tolerance](scalar x)
            {
                return (std::fabs(target - x) > tolerance);
            };
        }

        //- Test if value is 'less' than prescribed \p threshold
        static unary lessOp(scalar threshold)
        {
            return [threshold](scalar x) { return (x < threshold); };
        }

        //- Test if value is 'less_equal' to prescribed \p threshold
        static unary lessEqualOp(scalar threshold)
        {
            return [threshold](scalar x) { return (x <= threshold); };
        }

        //- Test if value is 'greater' than prescribed \p threshold
        static unary greaterOp(scalar threshold)
        {
            return [threshold](scalar x) { return (x > threshold); };
        }

        //- Test if value is 'greater_equal' to prescribed \p threshold
        static unary greaterEqualOp(scalar threshold)
        {
            return [threshold](scalar x) { return (x >= threshold); };
        }


    // Special purpose comparison methods

        //- Predicate that always returns true
        static unary trueOp()
        {
            return [](scalar) { return true; };
        }

        //- Predicate that always returns false
        static unary falseOp()
        {
            return [](scalar) { return false; };
        }

        //- Predicate that triggers a warning when used.
        //  \returns false
        static unary warnOp()
        {
            return [](scalar)
            {
                WarningIn("predicates::scalars::warnOp")
                    << "Using an invalid scalar unary predicate\n";
                return false;
            };
        }

        //- Predicate that triggers a warning when used.
        //- Reports the context in the message.
        //  \returns false
        static unary warnOp(std::string context)
        {
            return [context = std::move(context)](scalar)
            {
                WarningIn(context)
                    << "Using an invalid scalar unary predicate\n";
                return false;
            };
        }

        //- Predicate that triggers a warning when used.
        //- Reports the \p context, file/line in the message.
        //  \returns false
        static unary warnOp(std::string context, const char* file, int line)
        {
            return [context = std::move(context), file, line](scalar)
            {
                Foam::Warning(context, file, line)
                    << "Using an invalid scalar unary predicate\n";
                return false;
            };
        }


    // Combining operations

        //- Combine unary tests as an AND operation
        static unary andOp(const unary& test1, const unary& test2)
        {
            return [test1, test2](scalar x)
            {
                return (test1(x) && test2(x));
            };
        }

        //- Combine unary tests as an OR operation
        static unary orOp(const unary& test1, const unary& test2)
        {
            return [test1, test2](scalar x)
            {
                return (test1(x) || test2(x));
            };
        }


    // Factory for standard comparison methods

        //- Standard comparison method by type
        static unary operation
        (
            opType op,
            scalar operand,
            scalar tolerance = defaultTol()
        );

        //- Standard comparison method by name
        inline static unary operation
        (
            const word& opName,
            scalar operand,
            scalar tolerance = defaultTol()
        );

        //- Standard comparison method by (name, value)
        inline static unary operation
        (
            const Tuple2<word, scalar>& op,
            scalar tolerance = defaultTol()
        );

        //- Standard comparison method by (name, value)
        inline static unary operation
        (
            const std::pair<word, scalar>& op,
            scalar tolerance = defaultTol()
        );


    // Constructors

        //- Inherit constructors from List of unary functions
        using List<unary>::List;

        //- Construct from an initializer list of (opName operand) tuples
        explicit scalars
        (
            std::initializer_list<std::pair<const char*, scalar>> entries
        );

        //- Copy construct from a list of (opName operand) tuples
        explicit scalars(const UList<Tuple2<word, scalar>>& entries);

        //- Construct from Istream, from list of (opName operand) tuples
        explicit scalars(Istream& is);


    //- Destructor
    ~scalars() = default;


    // Member Functions

        //- Assign new predicates from an initializer list of
        //- (opName operand) tuples
        void assign
        (
            std::initializer_list<std::pair<const char*, scalar>> entries
        );

        //- Assign new predicates from a list of (opName operand) tuples
        void assign(const UList<Tuple2<word, scalar>>& entries);

        //- Inherit regular push_back methods (even if not particularly useful)
        using List<unary>::push_back;

        //- Append operation by operation \b type and operand value
        void push_back
        (
            opType op,
            scalar operand,
            scalar tolerance = defaultTol()
        );

        //- Append operation by operation \b name and operand value
        void push_back
        (
            const word& opName,
            scalar operand,
            scalar tolerance = defaultTol()
        );


    // Search

        //- Index of the first matching test for the value.
        //  Any occurrences before the start pos are ignored.
        //  Linear search.
        //  \return position in list or -1 if not found.
        label find(scalar value, label pos = 0) const;

        //- Index of the last match for the value.
        //  Any occurrences after the end pos are ignored.
        //  Linear search.
        //  \return position in list or -1 if not found.
        label rfind(scalar value, label pos = -1) const;

        //- Matches \em any condition in the list after the offset position.
        //  Linear search from start pos until the end of the list.
        //  Any occurrences before the start pos are ignored.
        //
        //  \return True if the value matches any condition in the (sub)list.
        inline bool contains(scalar value, label pos = 0) const;

        //- Matches \em any condition in the list.
        //
        //  \return True if the value matches any condition in the list.
        inline bool match(scalar value) const;

        //- Matches \em any condition in the list.
        //
        //  \return True if the value matches any condition in the list.
        inline bool matchAny(scalar value) const;

        //- Matches \em all conditions in the list.
        //
        //  \return True if the value matches all conditions in the list.
        inline bool matchAll(scalar value) const;

        //- Extract list indices for all matches.
        //
        //  \return The locations (indices) of the predicates that match
        inline labelList matching(scalar value) const;

        //- Extract list indices for all matches.
        //
        //  \return The locations (indices) in the input list
        //      where match() is true
        inline labelList matching
        (
            //! A list of scalar inputs to test
            const UList<scalar>& input,
            //! invert the matching logic
            const bool invert = false
        ) const;


    // Member Operators

        //- Identical to contains(), match(), for use as a predicate.
        bool operator()(scalar value) const { return matchAny(value); }


    // Housekeeping

        //- Same as contains()
        bool found(scalar value, label pos = 0) const
        {
            return this->contains(value, pos);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace predicates


// IOstream Operators

//- Suppressed write
Ostream& operator<<(Ostream& os, const predicates::scalars& list)
{
    return os;
}

//- Read list of (opName opValue) tuple
Istream& operator>>(Istream& is, predicates::scalars& list);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarPredicatesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
