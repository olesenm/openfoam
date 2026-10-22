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

Class
    Foam::bitSet

Description
    A bitSet stores bits (elements with only two states) in packed internal
    format and supports a variety of bit-set operations.
    Its behaviour is largely list-like, with some HashSet features.

SourceFiles
    bitSetI.H
    bitSet.C
    bitSetIO.C
    bitSet.txx

See also
    Foam::BitOps
    Foam::PackedList

\*---------------------------------------------------------------------------*/

#ifndef Foam_bitSet_H
#define Foam_bitSet_H

#include "PackedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;

/*---------------------------------------------------------------------------*\
                           Class bitSet Declaration
\*---------------------------------------------------------------------------*/

class bitSet
:
    public PackedList<1>
{
protected:

    // Logic/Set Operations

        //- The set difference
        //  \code
        //      A = (A - B)
        //      A = (A & !B)
        //      A = (A & ~B)
        //  \endcode
        //  A and B can have different sizes.
        //  Never changes the original set size.
        //  Non-overlapping parts are considered \em off.
        void minusEq(const bitSet& other);

        //- The set logical AND
        //  \code
        //      A = (A & B)
        //
        //  \endcode
        //  A and B can have different sizes.
        //  Never changes the original set size.
        //  Non-overlapping parts are considered \em off.
        void andEq(const bitSet& other);

        //- The set logical OR
        //  \code
        //      A = (A | B)
        //  \endcode
        //  A and B can have different sizes
        //
        //  The size grows to accommodate new \em on bits.
        void orEq(const bitSet& other);

        //- The set logical XOR
        //  \code
        //      A = (A ^ B)
        //  \endcode
        //  A and B can have different sizes.
        //
        //  The size grows to accommodate new \em on bits.
        void xorEq(const bitSet& other);

        //- Shift bits left by \p n positions,
        //- discarding bits pushed past the end and filling vacated
        //- bits with 0. The size is unchanged.
        void left_shiftEq(label n);

        //- Shift bits right by \p n positions,
        //- discarding bits pushed past the end and filling vacated
        //- bits with 0. The size is unchanged.
        void right_shiftEq(label n);

        //- Extract/copy a sub-section from \p source into \p target bitSet,
        //- resizing and overwriting any existing content.
        static void extract_copy
        (
            const bitSet& source,
            bitSet& target,
            const labelRange& range
        );

        //- Extract/copy a sub-section from \p source into \p target bitSet,
        //- resizing and overwriting any existing content.
        static void extract_copy
        (
            const bitSet& source,
            bitSet& target,
            label pos,
            label len = -1
        );


    // Constructors

        //- Move construct from the base PackedList<1>
        inline explicit bitSet(PackedList<1>&& list) noexcept;


public:

    //- A bitSet acts like a packed boolList
    typedef bool value_type;


    // Forward Declarations

        class reference;
        class const_iterator;
        typedef unsigned int const_reference;


    // Static Member Functions

        //- Return a null bitSet (reference to a nullObject).
        static const bitSet& null() noexcept
        {
            return NullObjectRef<bitSet>();
        }


    //- Declare type-name (with debug switch)
    ClassName("bitSet");


    // Constructors

        //- Default construct an empty, zero-sized bitSet
        inline constexpr bitSet() noexcept;

        //- Construct from Istream
        explicit bitSet(Istream& is);

        //- Construct with given size, with all bits set to 0
        inline explicit bitSet(const label n);

        //- Construct with given size and value for all elements
        inline bitSet(const label n, const bool val);

        //- Copy construct
        inline bitSet(const bitSet& bitset);

        //- Move construct
        inline bitSet(bitSet&& bitset);

        //- Construct a new bitSet by extracting the specified (unique)
        //- locations of an existing bitSet.
        bitSet(const bitSet& bitset, const labelUList& addr);

        //- Construct a new set by extracting the specified (unique)
        //- locations of an existing bitSet.
        template<class Addr>
        bitSet
        (
            const bitSet& bitset,
            const IndirectListBase<label, Addr>& addr
        );

        //- Construct a new bitSet by extracting the specified range
        //- locations of an existing bitSet.
        bitSet(const bitSet& bitset, const labelRange& range);

        //- Construct from a list of bools
        inline explicit bitSet(const UList<bool>& bools);

        //- Copy construct with given min size
        explicit bitSet(const label minSize, const bitSet& bitset);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        explicit bitSet(const label n, const labelRange& range);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        inline bitSet(const label n, const labelUList& locations);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        template<class Addr>
        inline bitSet
        (
            const label n,
            const IndirectListBase<label, Addr>& locations
        );

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        template<unsigned N>
        bitSet(const label n, const FixedList<label, N>& locations);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        inline bitSet(const label n, std::initializer_list<label> locations);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        explicit bitSet(const labelRange& range);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        inline explicit bitSet(const labelUList& locations);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        template<class Addr>
        inline explicit bitSet(const IndirectListBase<label, Addr>& locations);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        template<unsigned N>
        explicit bitSet(const FixedList<label, N>& locations);

        //- Clone
        inline autoPtr<bitSet> clone() const;


    // Member Functions

    // Query

        //- True if all bits in this bitset are set or if the set is \b empty.
        //  Returning true for an empty set may not seem intuitive, but
        //  conforms with boost definitions and std::all_of behaviour.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline bool all() const noexcept;

        //- True if any bits in this bitset are set.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline bool any() const noexcept;

        //- True if no bits in this bitset are set.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline bool none() const noexcept;

        //- Count number of bits set.
        //  \param on can be set to false to count the number of unset bits
        //     instead.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline unsigned int count(const bool on=true) const noexcept;

        //- True if any bits in the other bitset intersect (are the same).
        //
        //  \note Method name compatibility with boost::dynamic_bitset
        bool intersects(const bitSet& other) const;

        //- Test for \c true value at specified position.
        //- A no-op and returns \c false for out-of-range positions
        //- (safely/silently ignores negative positions).
        //  Never auto-vivify entries.
        //
        //  \note Method name compatibility with std::bitset
        bool test(label pos) const noexcept { return this->get(pos); }

        //- Test for \c true value at specified position.
        //- A no-op and returns \c false for out-of-range positions
        //- (safely/silently ignores negative positions).
        //  Never auto-vivify entries.
        //
        //  \note Method name compatibility with HashSet
        bool contains(label pos) const noexcept { return this->get(pos); }

        //- True if the bitset is unused (empty) or if the given position
        //- is \em on.
        //
        //  The combined test can be convenient when a bitSet is being used
        //  as an optional filter mechanism with an empty (zero-size)
        //  implying no filtering.
        inline bool empty_or_test(label i) const noexcept;

        //- Locate the first bit that is set, starting at (and including)
        //- the specified position (default: 0).
        //  \return the location or -1 if there are no further bits set.
        //
        //  \note Method name compatibility with boost::dynamic_bitset
        inline label find_first(label pos = 0) const noexcept;

        //- Locate the first bit that is unset.
        //  \return the location or -1 if the set is empty or all bits are on.
        inline label find_first_not(label pos = 0) const noexcept;

        //- Same as find_first_not()
        //  \note Method name compatibility with boost::dynamic_bitset
        label find_first_off(label pos = 0) const noexcept
        {
            return find_first_not(pos);
        }

        //  \return the location or -1 if there are no bits set.
        //
        //  \note Provided for symmetry with find_first()
        inline label find_last() const noexcept;

        //- Locate the last bit that is unset.
        //  \return the location or -1 if the set is empty or all bits are on.
        //
        //  \note Provided for symmetry with find_last(), find_first_not()
        inline label find_last_not() const noexcept;

        //- Locate the next bit set, starting one beyond the specified position
        //  \return the location or -1 if there are no further bits set.
        //
        //  \note Method name compatibility with boost::dynamic_bitset
        inline label find_next(label pos) const noexcept;

        //- Locate the next unset bit,
        //- starting one beyond the specified position
        //  \return the location or -1 if there are no further bits set.
        inline label find_next_not(label pos) const noexcept;

        //- Same as find_next_not()
        //  \note Method name compatibility with boost::dynamic_bitset
        label find_next_off(label pos) const noexcept
        {
            return find_next_not(pos);
        }

        //- The indices of the \a on bits as a sorted labelList.
        //  \param offset apply an optional offset to the values
        //
        //  \note Method name compatibility with HashSet
        labelList toc(const label offset = 0) const;

        //- Identical to toc(), which is always sorted.
        //  \param offset apply an optional offset to the values
        //
        //  \note Method name compatibility with HashSet
        inline labelList sortedToc(const label offset = 0) const;

        //- Return the bitset values as a boolList.
        //  When the output is a bool, this is more efficient than unpack()
        List<bool> values() const;

        //- Extract a sub-section \b copy as a new bitSet
        inline bitSet extract(label pos, label len = -1) const;

        //- Extract a sub-section \b copy as a new bitSet
        inline bitSet extract(const labelRange& range) const;

        //- Extract/copy a sub-section into the \p target bitSet,
        //- overwriting any existing content in \p target.
        inline void extract(bitSet& target, label pos, label len = -1) const;

        //- Extract/copy a sub-section into the \p target bitSet,
        //- overwriting any existing content in \p target.
        inline void extract(bitSet& target, const labelRange& range) const;


    // Assignment

        //- Assign all entries to the given value.
        inline void fill(const bool val);

        //- Copy assign all entries from a list of bools.
        void assign(const UList<bool>& bools);

        //- Copy/assign all values of \p source to the current set
        //- at the given \p offset, overwriting existing content in that range.
        //  If the offset + source ends exceeds the current size,
        //  the contents will grow, even if filled with zero bits.
        void assign(const label offset, const bitSet& source);


    // Setting single or multiple values

        //- Assign a single index/value, or a range with value
        using PackedList<1>::set;

        //- Set specified bits from another bitset.
        //  The current set size may grow to accommodate any new bits
        //  (auto-vivifies).
        inline void set(const bitSet& bitset);

        //- Test for \em True value at specified position
        //- and change the value at that position.
        //  Does auto-vivify for non-existent, non-zero entries.
        //
        //  \note Method name compatibility with std::bitset
        inline bool test_set(const label i, const bool val = true);


    // Unsetting single or multiple values

        //- Unset a single index, or a range
        using PackedList<1>::unset;

        //- Unset (subtract) the bits specified in the other bitset, which is
        //- a set difference corresponds to the logical operation
        //  \code
        //      A = (A & !B)
        //  \endcode
        //  The result is comparable to 'operator-='
        //  \endcode
        //
        //  A and B can have different sizes.
        //  Does not change the original set size.
        inline void unset(const bitSet& other);

        //- Unset (subtract) the bits specified in the other bitset,
        //- at the specified offset of the current.
        void unset(const label offset, const bitSet& other);


    // Edit

        //- Invert all bits in the addressable region
        inline void flip();

        //- Invert bit at the specified position.
        //  A no-op if the position is out-of-range
        inline void flip(const label pos);

        //- Resize to include the last \em on bit only.
        //  Functionally identical to resize(find_last()+1)
        inline void resize_last();

        //- Swap contents
        inline void swap(bitSet& bitset);

        //- Transfer the contents of the argument list into this list
        //- and annul the argument list.
        inline void transfer(bitSet& bitset);


    // Convenience methods

        //- Set the locations listed by the iterator range,
        //- auto-vivify entries if needed.
        //
        //  \return number of locations changed
        template<class InputIter>
        label setMany(InputIter first, InputIter last);

        //- Set the listed locations to 1.
        //  Does auto-vivify for non-existent entries.
        //
        //  \return number of locations changed
        inline label set(const labelUList& locations);

        //- Set the listed locations to 1.
        //  Does auto-vivify for non-existent entries.
        //
        //  \return number of locations changed
        template<class Addr>
        inline label set(const IndirectListBase<label, Addr>& locations);

        //- Set the listed locations to 1.
        //  Does auto-vivify for non-existent entries.
        //
        //  \return number of locations changed
        template<unsigned N>
        label set(const FixedList<label, N>& locations);

        //- Unset the locations listed by the iterator range,
        //- never auto-vivify entries.
        //
        //  \return number of locations changed
        template<class InputIter>
        label unset(InputIter first, InputIter last);

        //- Unset the listed locations, never auto-vivifies.
        //
        //  \return number of locations changed
        inline label unset(const labelUList& locations);

        //- Unset the listed locations, never auto-vivifies.
        //
        //  \return number of locations changed
        template<class Addr>
        inline label unset(const IndirectListBase<label, Addr>& locations);

        //- Unset the listed locations, never auto-vivifies.
        //
        //  \return number of locations changed
        template<unsigned N>
        label unset(const FixedList<label, N>& locations);


    // Access helpers

        //- A reference supporting read/write access to an entry
        class reference
        :
            public PackedList<1>::reference
        {
        protected:

            friend class bitSet;        // Access for parent
            void operator&() = delete;  // Refuse to provide its address

            //- Construct by taking reference of block from within
            //- the list and the specified position
            inline reference(bitSet* parent, const label pos) noexcept;

        public:

            //- Copy construct
            reference(const reference&) noexcept = default;

            //- Move construct
            reference(reference&&) noexcept = default;

            //- Flip the bit at current position
            inline void flip() noexcept;

            //- Value assignment
            inline void operator=(const reference& other) noexcept;

            //- Value assignment
            inline void operator=(const unsigned int val) noexcept;

            //- Conversion operator.
            operator unsigned int () const noexcept { return this->get(); }
        };


    // Iteration

        //- A const_iterator for iterating across \a on values
        class const_iterator
        {
            friend class bitSet;

            //- The parent being iterated
            const bitSet* set_;

            //- Global position of the current \a on bit
            label pos_;

            //- Default construct - an end iterator
            inline constexpr const_iterator() noexcept;

            //- Construct begin iterator. Uses find_first
            inline const_iterator(const bitSet* parent) noexcept;

            //- Construct iterator, starting at or beyond the given position.
            //- Uses find_next()
            inline const_iterator(const bitSet* parent, label pos) noexcept;

        public:

            //- Return the current \a on position
            inline label operator*() const noexcept;

            //- Move to the next \a on position
            inline const_iterator& operator++() noexcept;

            inline bool operator==(const const_iterator& iter) const noexcept;
            inline bool operator!=(const const_iterator& iter) const noexcept;
        };


        //- Iterator set to the position of the first \a on bit
        inline const_iterator begin() const noexcept;

        //- Iterator set to the position of the first \a on bit
        inline const_iterator cbegin() const noexcept;

        //- Iterator set to the position of the first \a on bit that occurs
        //- at or beyond the given position
        inline const_iterator begin(label pos) const noexcept;

        //- Iterator set to the position of the first \a on bit that occurs
        //- at or beyond the given position
        inline const_iterator cbegin(label pos) const noexcept;

        //- Iterator beyond the end of the bitSet
        inline const_iterator end() const noexcept;

        //- Iterator beyond the end of the bitSet
        inline const_iterator cend() const noexcept;


    // Member Operators

        //- Test value at specified position, same as test()
        //  Enables use as a predicate
        inline bool operator()(const label pos) const noexcept;

        //- Identical to get() - get value at index.
        //  Never auto-vivify entries.
        inline unsigned int operator[](const label i) const noexcept;

        //- Non-const access to value at index.
        //  Fatal for out-of-range indices
        inline reference operator[](const label i);

        //- Copy assignment
        inline bitSet& operator=(const bitSet& bitset);

        //- Move assignment
        inline bitSet& operator=(bitSet&& bitset);

        //- Assign all entries to the given value. fill()
        inline bitSet& operator=(const bool val);

        //- Bitset complement, returns a copy with all its bits flipped
        inline bitSet operator~() const;

        //- Bitwise-AND all the bits in other with the bits in this bitset.
        //  The operands may have dissimilar sizes,
        //  never changes the original set size.
        //  Non-overlapping elements are considered \em off.
        inline bitSet& operator&=(const bitSet& other);

        //- Bitwise-OR operator - similar to the set() method.
        //  The operands may have dissimilar sizes,
        //  the set grows to accommodate new \em on bits.
        inline bitSet& operator|=(const bitSet& other);

        //- Bitwise-XOR operator - retains unique entries.
        //  The operands may have dissimilar sizes,
        //  the set grows to accommodate new \em on bits.
        inline bitSet& operator^=(const bitSet& other);

        //- Remove entries from this list - identical to the unset() method.
        //  The operands may have dissimilar sizes,
        //  never changes the original set size.
        inline bitSet& operator-=(const bitSet& other);

        //- Shift bits to the left by \p n positions.
        //  For each bit in the bitset, the bit at position pos takes on
        //  the previous value of the bit at position (pos - n), or zero if
        //  no such bit exists.
        //  Does not change the size.
        inline bitSet& operator<<=(label n);

        //- Shift bits to the right by \p n positions.
        //  For each bit in the bitset, the bit at position pos takes on
        //  the previous value of the bit at position (pos + n), or zero if
        //  no such bit exists.
        //  Does not change the size.
        inline bitSet& operator>>=(label n);

        //- Return a shifted copy of \c *this
        inline bitSet operator<<(label n) const;

        //- Return a shifted copy of \c *this
        inline bitSet operator>>(label n) const;


    // Reading/writing


        //- Read contents in "compact" format, which is the listing
        //- of the 'on' bits (toc)
        //  This function mostly useful for relatively sparse sets.
        //  \note currently ASCII only.
        Istream& readListToc(Istream& is);

        //- Write contents in "compact" format, which is the listing
        //- of the 'on' bits (toc), with line-breaks in ASCII guided
        //- by shortLen parameter ('0' suppresses line-breaks entirely).
        //  This function mostly useful for relatively sparse sets.
        //  \note currently ASCII only.
        Ostream& writeListToc(Ostream& os, label shortLen=0) const;

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<bitSet> info() const noexcept
        {
            return *this;
        }


    // Parallel Operations

        //- Broadcast the contents
        void broadcast
        (
            //! The UPstream communicator (default is worldComm)
            int communicator = -1,
            //! The broadcast root
            int root = 0
        );

        //- Inplace \c bit_and parallel reduction
        void reduceAnd
        (
            //! The UPstream communicator (default is worldComm)
            int communicator = -1,
            //! False: sizes already consistent. True: adjust (shrink) sizes.
            bool syncSizes = true
        );

        //- Inplace \c bit_or parallel reduction
        void reduceOr
        (
            //! The UPstream communicator (default is worldComm)
            int communicator = -1,
            //! False: sizes already consistent. True: adjust (extend) sizes.
            bool syncSizes = true
        );

        //- Gather individual values into bitSet locations.
        //  On master, the resulting bitSet has size == nProcs,
        //  otherwise zero length.
        //  \n
        //  For \b non-parallel :
        //  the length of the returned set is 1 with localValue.
        static bitSet gatherValues
        (
            //! The processor-local value
            bool localValue,
            //! The UPstream communicator (default is worldComm)
            int communicator = -1
        );

        //- Allgather individual values into bitSet locations.
        //  The resulting bitSet has size nProcs, identical on all ranks.
        //  Behaves like Pstream::allGatherValues() but returning a bitSet.
        static bitSet allGather
        (
            //! The processor-local value
            bool localValue,
            //! The UPstream communicator (default is worldComm)
            int communicator = -1
        );


    // Housekeeping

        //- Same as contains()
        bool found(const label pos) const { return this->contains(pos); }

        //- \deprecated(2020-11) use fill()
        FOAM_DEPRECATED_FOR(2020-11, "fill()")
        void assign(const unsigned int val) { this->fill(val); }

        //- \deprecated(2026-01) use two-parameter version
        FOAM_DEPRECATED_FOR(2026-01, "broadcast(comm, root)")
        void broadcast(std::pair<int,int> comm_root)
        {
            this->broadcast(comm_root.first, comm_root.second);
        }

        //- Ensure the addressable range does not exceed maxSize
        FOAM_DEPRECATED_FOR(2026-01, "<removed>")
        void bound(const label maxSize) = delete;

        FOAM_DEPRECATED_FOR(2026-01, "<removed>")
        void bound(const bitSet& other) = delete;

        //- Ensure that minSize is covered by the bitSet
        FOAM_DEPRECATED_FOR(2026-01, "<removed>")
        void extend(label minSize) = delete;

        FOAM_DEPRECATED_FOR(2026-01, "<removed>")
        void extend(const bitSet& other) = delete;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Hashing for bitSet data
template<> struct Hash<bitSet> : bitSet::hasher {};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Write bitset to Ostream with 40 items per line.
Ostream& operator<<(Ostream& os, const bitSet& bitset);

//- Output bitset information
Ostream& operator<<(Ostream& os, const InfoProxy<bitSet>& iproxy);

//- Bitwise-AND of two bitsets.
//  See bitSet::operator&= for more details.
inline bitSet operator&(const bitSet& a, const bitSet& b);

//- Bitwise-OR of two bitsets
//  See bitSet::operator|= for more details.
inline bitSet operator|(const bitSet& a, const bitSet& b);

//- Bitwise-XOR of two bitsets to form a unique bit-set
//  See bitSet::operator^= for more details.
inline bitSet operator^(const bitSet& a, const bitSet& b);

//- Bitwise difference (subset) of two bitsets to form a unique bit-set
//  See bitSet::operator-= for more details.
inline bitSet operator-(const bitSet& a, const bitSet& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bitSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "bitSet.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
