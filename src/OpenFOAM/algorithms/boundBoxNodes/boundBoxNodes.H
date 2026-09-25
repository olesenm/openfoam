/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Different binary tree addressing types.

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundBoxNodes_H
#define Foam_boundBoxNodes_H

#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace meshTools
{

// Forward Declarations
template<class BBoxType> class BoundBoxBinaryNode;

/*---------------------------------------------------------------------------*\
                    Class binaryNodeAddressing Declaration
\*---------------------------------------------------------------------------*/

//- Binary tree leaf/interior addressing.
//- Holds left/right addressing, which is also reused for managing
//- leaf subset indexing.
class binaryNodeAddressing
{
    // Private Data

        //- Left sub-node index (interior) or content start (leaf)
        label a_ = 0;

        //- Right sub-node index (interior) or content size (leaf),
        //- in which case it is the negative (size+1) value.
        label b_ = -1;


protected:

    // Protected Member Functions

        //- Ordering based on (a_, b_) values with interior nodes
        //- ordered before leaf nodes.
        //
        //  Follows the standard three-way comparison convention:
        //  negative when *this < other, positive when *this > other,
        //  or 0 when the derived class must break the tie
        //  (eg, by comparing bounding boxes).
        label compare_keys(const binaryNodeAddressing& other) const noexcept
        {
            if (this->is_interior() != other.is_interior())
            {
                // Mixed types. Sort interior before leaf.
                return (this->is_interior() ? -1 : 1);
            }

            label diff = (this->a_ - other.a_);

            if (diff == 0)
            {
                diff =
                (
                    this->is_leaf()
                  ? (this->size() - other.size())
                  : (this->b_ - other.b_)
                );
            }

            return diff;
        }

        //- Equality based on (a_, b_) alone
        bool equal_keys(const binaryNodeAddressing& other) const noexcept
        {
            return (this->a_ == other.a_ && this->b_ == other.b_);
        }


public:

    // Public Data Types

        //- The bound box type (ie, none)
        typedef void box_type;

        //- The indexing type
        typedef label index_type;


    // Methods

        //- Is this an interior node?
        bool is_interior() const noexcept { return (b_ >= 0); }

        //- Is this a leaf node?
        bool is_leaf() const noexcept { return (b_ < 0); }

        //- Left sub-node index (interior nodes only)
        label left() const noexcept { return a_; }

        //- Right sub-node index (interior nodes only)
        label right() const noexcept { return b_; }

        //- Start of the leaf content range into indices()
        //- (leaf nodes only)
        label start() const noexcept { return a_; }

        //- Size of the leaf content range into indices()
        //- (leaf nodes only)
        label size() const noexcept
        {
            return (b_ < 0 ? (-b_ - 1) : 0);
        }

        //- Set this node as an interior node with the given
        //- sub-node indices (both must be > 0)
        void set_interior(label leftIndex, label rightIndex) noexcept
        {
            a_ = leftIndex;
            b_ = rightIndex;
        }

        //- Set this node as a leaf with the given content range
        void set_leaf(label start, label len) noexcept
        {
            a_ = start;
            b_ = -(len + 1);
        }

        //- Return read-only slice into values (for leaf)
        template<class Type>
        const SubList<Type> slice(const UList<Type>& values) const
        {
            if (is_leaf())
            {
                return values.slice(start(), size());
            }
            else
            {
                return SubList<Type>();
            }
        }


    // IO Operations

        //- Read node addressing from Istream
        Istream& readNode(Istream& is)
        {
            is  >> a_ >> b_;
            return is;
        }

        //- Write node addressing to Ostream
        Ostream& writeNode(Ostream& os) const
        {
            os  << a_ << ' ' << b_;
            return os;
        }


    // Operators

        //- Less-than (for sorted output)
        bool operator<(const binaryNodeAddressing& other) const noexcept
        {
            return (this->compare_keys(other) < 0);
        }

        //- Compare equality
        bool operator==(const binaryNodeAddressing& other) const noexcept
        {
            return (this->equal_keys(other));
        }

        //- Compare inequality
        bool operator!=(const binaryNodeAddressing& other) const noexcept
        {
            return !(*this == other);
        }
};


/*---------------------------------------------------------------------------*\
                     Class BoundBoxBinaryNode Declaration
\*---------------------------------------------------------------------------*/

//- A tree node, which is either leaf or interior but not both.
//
//  Interior nodes have a left and right sub-node index.
//  Leaf nodes reference a contiguous range [start, start+size)
//  into the indices() addressing.
template<class BBoxType>
class BoundBoxBinaryNode
:
    public binaryNodeAddressing
{
    // Private Data

        //- Bounding box enclosing all content under this node
        BBoxType bbox_;


public:

    // Public Data Types

        //- The type of bounding box
        typedef BBoxType box_type;


    // Methods

        //- Bounding box enclosing all content under this node
        const box_type& bb() const noexcept { return bbox_; }

        //- Set bounding box enclosing all content under this node
        void set_bb(const box_type& bb) noexcept { bbox_ = bb; }


    // IO Operations

        //- Read a bound box node from Istream
        Istream& readNode(Istream& is)
        {
            is  >> bbox_;
            return this->binaryNodeAddressing::readNode(is);
        }

        //- Write a bound box node to Ostream
        Ostream& writeNode(Ostream& os) const
        {
            os  << bbox_ << ' ';
            return this->binaryNodeAddressing::writeNode(os);
        }


    // Operators

        //- Less-than (for sorted output)
        bool operator<(const BoundBoxBinaryNode& other) const
        {
            const auto diff = this->compare_keys(other);

            if (diff == 0)
            {
                // Tie-break: compare boxes
                return (this->bbox_ < other.bbox_);
            }

            return (diff < 0);
        }

        //- Compare equality
        bool operator==(const BoundBoxBinaryNode& other) const
        {
            return
            (
                this->equal_keys(other)
             && this->bbox_ == other.bbox_
            );
        }

        //- Compare inequality
        bool operator!=(const BoundBoxBinaryNode& other) const
        {
            return !(*this == other);
        }
};


//- A tree node, which is either leaf or interior but not both.
//
//  Interior nodes have a left and right sub-node index.
//  Leaf nodes reference a contiguous range [start, start+size)
//  into the indices() addressing.
template<>
class BoundBoxBinaryNode<void> : public binaryNodeAddressing
{
public:

    // Public Data Types

        //- No normal box type
        typedef void box_type;
};

} // End namespace meshTools


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Binary node addressing is contiguous (pair of labels)
template<>
struct is_contiguous<meshTools::binaryNodeAddressing> : std::true_type {};

//- Binary node (void) is contiguous
template<>
struct is_contiguous<meshTools::BoundBoxBinaryNode<void>> : std::true_type {};

//- Binary node is contiguous when boundBox is
template<class BBType>
struct is_contiguous<meshTools::BoundBoxBinaryNode<BBType>>
:
    is_contiguous<BBType>
{};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

inline Istream& operator>>
(
    Istream& is,
    meshTools::binaryNodeAddressing& node
)
{
    return node.readNode(is);
}

inline Ostream& operator<<
(
    Ostream& os,
    const meshTools::binaryNodeAddressing& node
)
{
    return node.writeNode(os);
}


template<class BBoxType>
Istream& operator>>
(
    Istream& is,
    meshTools::BoundBoxBinaryNode<BBoxType>& node
)
{
    return node.readNode(is);
}

template<class BBoxType>
Ostream& operator<<
(
    Ostream& os,
    const meshTools::BoundBoxBinaryNode<BBoxType>& node
)
{
    return node.writeNode(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
