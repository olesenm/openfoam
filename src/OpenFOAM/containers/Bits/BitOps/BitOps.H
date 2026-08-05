/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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
    Foam::BitOps

Description
    Various bit-wise operations and adaptor methods for containers
    that are somewhat similar to bitSet (eg, boolList, labelHashSet).

    The population count uses the Hamming weight
    (http://en.wikipedia.org/wiki/Hamming_weight).

Namespace
    Foam::BitSetOps

Description
    Factory and other methods for bitSet.
    Adaptor methods for other containers that are somewhat similar to
    bitSet (eg, boolList, labelHashSet).

\*---------------------------------------------------------------------------*/

#ifndef Foam_BitOps_H
#define Foam_BitOps_H

#include "UList.H"
#include "HashSet.H"
#include "Ostream.H"
#include <limits>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;
template<class T> class List;

/*---------------------------------------------------------------------------*\
                         Namespace BitOps Declaration
\*---------------------------------------------------------------------------*/

namespace BitOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Count number of 'true' entries.
//  \param val can be set to false to count the number of false values instead
//  For compatibility with bitSet::count()
inline unsigned int count(const UList<bool>& bools, const bool val=true)
{
    return std::count(bools.begin(), bools.end(), val);
}

//- True if all entries are 'true' or if the list is empty.
//  For compatibility with bitSet::all()
inline bool all(const UList<bool>& bools) { return bools.all(); }

//- True if any entries are 'true'.
//  For compatibility with bitSet::any()
inline bool any(const UList<bool>& bools) { return bools.any(); }

//- True if no entries are 'true'.
//  For compatibility with bitSet::none()
inline bool none(const UList<bool>& bools) { return bools.none(); }


//- Set the listed locations (assign 'true').
//  Does auto-vivify for non-existent entries.
//
//  For compatibility with bitSet::set(labelUList)
void set(List<bool>& bools, const labelUList& locations);

//- Set the specified range 'on' in a boolList.
//  For compatibility with bitSet::set(labelRange)
void set(List<bool>& bools, const labelRange& range);

//- Set the specified range in a labelHashSet.
//  For compatibility with bitSet::set(labelRange)
void set(labelHashSet& hashset, const labelRange& range);

//- Forward to bitSet::set(labelRange)
void set(bitSet& bitset, const labelRange& range);


//- Unset the listed locations (assign 'false').
//  No auto-vivify of non-existent entries.
//
//  For compatibility with bitSet::set(labelUList)
void unset(UList<bool>& bools, const labelUList& locations);

//- Unset the specified range 'on' in a boolList.
//  No auto-vivify of non-existent entries.
//
//  For compatibility with bitSet::unset(labelRange)
void unset(UList<bool>& bools, const labelRange& range);

//- Unset the specified range in a labelHashSet.
//  No auto-vivify of non-existent entries.
//
//  For compatibility with bitSet::unset(labelRange)
void unset(labelHashSet& hashset, const labelRange& range);

//- Forward to bitSet::unset(labelRange)
void unset(bitSet& bitset, const labelRange& range);


//- Construct a selection list of bools (all false) with the given pre-size,
//- subsequently add specified locations as true,
//- auto-vivify entries if needed.
//  Similar to bitSet construction from locations
//
//  \return a List of bools
List<bool> select(const label n, const labelUList& locations);

//- Construct an auto-sized selection list of bools (all false),
//- and populate the specified locations as true.
//  Similar to bitSet construction from locations
//
//  \return a List of bools
List<bool> select(const labelUList& locations);

//- Forward to bitSet::count()
unsigned int count(const bitSet& bitset, const bool on=true);


//- Return the (sorted) values corresponding to 'true' entries.
//  Similar to bitSet::toc()
//
//  \return a List of labels
List<label> toc(const UList<bool>& bools);

//- Return the (sorted) values corresponding to 'true' entries.
//  Similar to bitSet::sortedToc() and labelHashSet::sortedToc()
//
//  \return a List of labels
List<label> sortedToc(const UList<bool>& bools);


//- Forward to bitSet::toc(), the sorted values of the 'on' entries.
//
//  \return a List of labels
List<label> toc(const bitSet& bitset);

//- Forward to bitSet::sortedToc(), the sorted values of the 'on' entries.
//
//  \return a List of labels
List<label> sortedToc(const bitSet& bitset);


//- Forward to labelHashSet::sortedToc(), the sorted values of the 'on' entries.
//
//  \return a List of labels
List<label> toc(const labelHashSet& hashset);

//- Forward to labelHashSet::sortedToc(), the sorted values of the 'on' entries.
//
//  \return a List of labels
List<label> sortedToc(const labelHashSet& hashset);


//- Count number of bits of an unsigned integral type
template<class UIntType>
inline unsigned int bit_count(UIntType x) noexcept
{
    static_assert(std::is_unsigned_v<UIntType>, "Expecting unsigned type");

    if constexpr (sizeof(UIntType) == sizeof(uint32_t))
    {
        // Count bits in a 32-bit value (Hamming weight method)

        x -= (x >> 1) & 0x55555555;
        x = (x & 0x33333333) + ((x >> 2) & 0x33333333);

        return
        (
            (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24
        );
    }
    else if constexpr (sizeof(UIntType) == sizeof(uint64_t))
    {
        // Count bits in a 64-bit value (Hamming weight method)

        x -= (x >> 1) & 0x5555555555555555;
        x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);

        return unsigned
        (
            (((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101) >> 56
        );
    }
    else
    {
        // Manual bit counting
        unsigned int n = 0u;

        for (; x; ++n) { x &= (x-1); }
        return n;
    }
}


//- Repeat a value of the given BitWidth into the destination output type.
//
//  \note when BitWidth is 1, it is better to do directly.
//  \code
//     (val ? ~0u : 0u)
//  \endcode
template<class UIntType, unsigned BitWidth>
inline constexpr UIntType repeat_value(unsigned val) noexcept
{
    static_assert
    (
        BitWidth && std::numeric_limits<UIntType>::digits >= BitWidth,
        "BitWidth too large for target output"
    );
    static_assert(std::is_unsigned_v<UIntType>, "Expecting unsigned type");

    if constexpr (BitWidth <= 1)
    {
        // Width == single bit
        return (val ? ~UIntType(0) : UIntType(0));
    }
    else
    {
        if (!val)
        {
            // Fast path for zero value
            return UIntType(0);
        }

        // Max value for a single element
        constexpr UIntType mask = ((UIntType(1) << BitWidth) - UIntType(1));

        // The first occurrence
        UIntType fillval = ((val >= mask) ? mask : val);

        // How many fit into the target
        constexpr unsigned nrepeat
            = (std::numeric_limits<UIntType>::digits / BitWidth);

        // Repeated
        for (unsigned i = 1; i < nrepeat; ++i)
        {
            fillval |= (fillval << BitWidth);
        }

        return fillval;
    }
}


//- Print 0/1 bits in the (unsigned) integral type
template<class UIntType>
inline Ostream& print(Ostream& os, UIntType value, char off='0', char on='1')
{
    if (os.format() == IOstreamOption::BINARY)
    {
        // Perhaps not the most sensible, but the only thing we currently have.
        os << label(value);
    }
    else
    {
        // Starting from most significant bit - makes for easy reading.
        for
        (
            unsigned test = (1u << (std::numeric_limits<UIntType>::digits-1));
            test;
            test >>= 1u
        )
        {
            os << ((value & test) ? on : off);
        }
    }

    return os;
}


//- An (unsigned) integral type adapter, for output of bit values
template<class UIntType>
struct bitInfo
{
    typedef UIntType value_type;
    value_type value;

    //- Null constructible as zero
    constexpr bitInfo() noexcept : value(0) {}

    //- Value construct
    explicit bitInfo(UIntType val) noexcept : value(val) {}

    //- Conversion to base type
    operator UIntType () const noexcept { return value; }

    //- Conversion to base type
    operator UIntType& () noexcept { return value; }
};

} // End namespace BitOps


/*---------------------------------------------------------------------------*\
                        Namespace BitSetOps Declaration
\*---------------------------------------------------------------------------*/

namespace BitSetOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Create a bitSet with length n with the specified \a on locations.
//  The resulting bitSet is guaranteed to have \b exactly the specified length,
//  any values or positions larger than n-1 are silently ignored.
//
//  \param n the size of the output bitSet
//  \param locations the list of positions corresponding to an \a on bit.
//  \param on the value for on. Set as false to invert the logic.
//
//  \return a bitset
bitSet create
(
    const label n,
    const labelHashSet& locations,
    const bool on = true
);


//- Create a bitSet with length n with the specified \a on locations.
//  The resulting bitSet is guaranteed to have \b exactly the specified length,
//  any values or positions larger than n-1 are silently ignored.
//
//  \param n the size of the output bitSet
//  \param locations the list of positions corresponding to an \a on bit.
//  \param on the value for on. Set as false to invert the logic.
//
//  \return a bitset
bitSet create
(
    const label n,
    const labelUList& locations,
    const bool on = true
);


//- Create a bitSet with length n with the specified \a on locations
//- when the list values are equal to the select value.
//
//  The resulting bitSet is guaranteed to have \b exactly the specified length,
//  any values or positions larger than n-1 are silently ignored.
//
//  \param n the size of the output bitSet
//  \param select the value to select as 'on'
//  \param values the values to scan for 'select'
//  \param on the value for on. Set as false to invert the logic.
//
//  \return a bitset
bitSet create
(
    const label n,
    const label select,
    const labelUList& values,
    const bool on = true
);

} // End namespace BitSetOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- Print 0/1 bits of an (unsigned) integral type via an adapter
template<class UIntType>
inline Ostream& operator<<(Ostream& os, const BitOps::bitInfo<UIntType>& info)
{
    BitOps::print(os, info.value);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
