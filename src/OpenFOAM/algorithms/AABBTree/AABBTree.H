/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::AABBTree

Description

    The class provides simple construction of axis-aligned bounding boxes
    (AABB) for objects (eg, Foam:faces, Foam:cells) or points only.

    Although the name suggests a tree structure, the current consumers
    of the output only require the leaf contents. The class is thus designed
    not to retain any tree navigation and only deliver the leaf elements.

    When constructed with objects such as Foam::face, Foam::cell etc, the
    object centres (or an estimated value of them) are used with a kd-tree
    to recursively split in the largest axis of current bounding box. After
    all leafs have been created, they are subsequently expanded to include
    all of the object points (eg, face vertices).

SourceFiles
    AABBTree.cxx
    AABBTree.txx
    AABBTreeI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_AABBTree_H
#define Foam_AABBTree_H

#include "treeBoundBox.H"
#include "boundBoxNodes.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class AABBTreeBase;

Istream& operator>>(Istream&, AABBTreeBase&);
Ostream& operator<<(Ostream&, const AABBTreeBase&);

/*---------------------------------------------------------------------------*\
                        Class AABBTreeBase Declaration
\*---------------------------------------------------------------------------*/

//- Template invariant part for AABBTree
class AABBTreeBase
{
public:

    // Public Data Types

        //- Different splitting approaches
        enum class algorithm
        {
            //! split longest bounding box axis into equal number of elements
            median_split,
            //! split longest bounding box axis into equal geometric sizes
            geometric_split
        };

        //- The type for bounding box binary tree nodes
        using node = meshTools::BoundBoxBinaryNode<treeBoundBox>;


    // Properties

        //- Default splitting approach
        static constexpr algorithm defaultAlgorithm() noexcept
        {
            return algorithm::median_split;
        }

        //- The default max tree depth/level
        static constexpr int defaultMaxLevel() noexcept { return 3; }

        //- The default min leaf size
        static constexpr int defaultMinLeafSize() noexcept { return 100; }

        //- Node caching used for queries?
        static constexpr bool cacheNodes() noexcept { return true; }


protected:

    // Protected Data Types

        //- Control parameters for building bounding boxes
        struct control
        {
            //! Maximum tree level
            int maxLevel = defaultMaxLevel();

            //! Minimum elements per leaf
            int minLeafSize = defaultMinLeafSize();

            //! Splitting method
            algorithm algo = defaultAlgorithm();

            //- Reset to defaults
            void reset() noexcept
            {
                maxLevel = defaultMaxLevel();
                minLeafSize = defaultMinLeafSize();
                algo = defaultAlgorithm();
            }
        };


    // Protected Member Functions

        //- True if type is not indexable
        template<class T>
        static constexpr bool nonIndexable() noexcept
        {
            // OR (!Foam::is_ulist_v<T> && !Foam::is_fixedlist_v<T>);
            return
            (
                std::is_void_v<T>
             || std::is_arithmetic_v<T>
             || is_vectorspace_v<T>
            );
        }

private:

    // Control Parameters

        //- Target occupancy (for auto-size mode)
        static constexpr int targetOccupancy = 64;

        //- Max depth level (for auto-size mode)
        static constexpr int maxDepthLimit = 12;


    // Private Data

        //- The bounding boxes for all leafs
        List<treeBoundBox> leafBoxes_;

        //- Permutation of original object indices, grouped by leaf.
        labelList indices_;

        //- Offset (addressing) table - for the leaves, slices into indices_
        labelList offsets_;

        //- Binary tree addressing with bounding boxes
        mutable std::unique_ptr<List<node>> nodes_;


    // Private Member Functions

        //- Reset (re-index) the nodes
        void calcNodeAddressing() const;

        //- Make reasonable guesses for level and occupancy limits
        static void autosize
        (
            control& ctrl,          //!< [in,out] control parameters
            label total,            //!< total number of objects
            int maxLevel = -1,      //!< requested max level
            int minLeafSize = -1    //!< requested min leaf size
        );

        //- Recursively split offsets_ by 2 to create addressing nodes.
        //- A leaf is reached when last=first+1.
        template<class NodeType>
        label reindex_impl
        (
            const label first,
            const label last,
            DynamicList<NodeType>& nodesToBuild
        ) const;


    // Queries

        //- Depth-first search for a matching leaf bounding box.
        //  Note: requires navigation nodes.
        template<class Predicate>
        bool found_first_of(label nodei, const Predicate& pred) const;


protected:

    // Protected Member Functions

        //- Inflate the leaf bounding boxes by given factor
        void inflate(scalar factor = 0.01);

        //- Inflate the leaf bounding boxes to also enclose the full
        //- extents of each object (as opposed to just its centre)
        template<class Type>
        void object_inflate
        (
            const UList<Type>& objects,
            const UList<point>& objectPoints,
            scalar factor = 0.01
        );

        //- Recursively kd-split [first,last) of indices_, accumulating
        //- leaf begin/end offsets. On the top-level (currLevel == 0) call,
        //- also transfers the result into offsets_ and builds leafBoxes_.
        void build_impl
        (
            const label first,
            const label last,
            const UList<point>& points,
            DynamicList<label>& offsetsBuilt,
            const int currLevel,
            const control& ctrl
        );

        //- Determine a representative centre point for each object.
        //- For non-indexable (primitive) types the objectPoints are
        //- simply reused directly (eg, when the type has no point loop).
        template<class Type>
        static tmp<pointField> objectCentres
        (
            const UList<Type>& objects,
            const pointField& objectPoints
        );

        //- Clear the tree and rebuild all addressing from a flat list
        //- of representative points (eg, object centres).
        //
        //  No inflation of the leaf boxes at all.
        void reset
        (
            const label nObjects,
            const UList<point>& points,
            bool medianSplit = true,                //!< median/geometric split
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

public:

    // Constructors

        //- Default construct
        constexpr AABBTreeBase() noexcept = default;

        //- Copy construct
        AABBTreeBase(const AABBTreeBase&);

        //- Move construct
        AABBTreeBase(AABBTreeBase&&);

        //- Construct by splitting a plain list of points
        explicit AABBTreeBase
        (
            const pointField& points,
            bool medianSplit = true,                //!< median/geometric split
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct with object list, centroids and vertices - using
        //- the object centres to perform the splitting.
        //
        //  After initial construction, the leaf boxes are inflated
        //  to the full object point extents.
        template<class Type>
        AABBTreeBase
        (
            //! The list of objects to scan
            const UList<Type>& objects,
            //! The reference centroid for each object
            const pointField& objectCentres,
            //! The referenced points (eg, vertices) for the objects
            const pointField& objectPoints,
            bool medianSplit = true,                //!< median/geometric split
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct with object list and vertices - using
        //- estimates of the object centres to perform the splitting.
        //
        //  After initial construction, the leaf boxes are inflated
        //  to the full object point extents.
        template<class Type>
        AABBTreeBase
        (
            //! The list of objects to scan
            const UList<Type>& objects,
            //! The referenced points (eg, vertices) for the objects
            const pointField& objectPoints,
            bool medianSplit = true,                //!< median/geometric split
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );


    // Member Functions

        //- Clear leafs and addressing
        void clear();

        //- Clear addressing (indices, offsets, etc)
        void clearAddressing();

        //- Clear node addressing
        void clearNodes();

        //- The number of leaves in the tree
        inline label nLeafs() const noexcept;

        //- Has binary tree addressing?
        inline bool hasNodes() const noexcept;

        //- List of all nodes (demand-driven)
        const List<node>& nodes() const;

        //- Overall bounding box
        treeBoundBox bb() const;

        //- The bounding boxes for each leaf (non-const access)
        List<treeBoundBox>& leafBoxes() noexcept
        {
            return leafBoxes_;
        }

        //- The bounding boxes for each leaf
        const List<treeBoundBox>& leafBoxes() const noexcept
        {
            return leafBoxes_;
        }

        //- Discard any tree addressing and
        //- return the leaf bounding boxes
        List<treeBoundBox> release_boxes()
        {
            clearAddressing();
            return std::move(leafBoxes_);
        }

        //- The permutation of original object indices, grouped by leaf
        labelList& indices() noexcept { return indices_; }

        //- The permutation of original object indices, grouped by leaf
        const labelList& indices() const noexcept { return indices_; }

        //- The begin/end offsets into indices() for each leaf
        labelList& offsets() noexcept { return offsets_; }

        //- The begin/end offsets into indices() for each leaf
        const labelList& offsets() const noexcept { return offsets_; }

        //- Return the object indices for the specified leaf
        inline const labelList::subList indices(label leaf) const;

        //- Return the range of object indices for the specified leaf
        inline labelRange range(label leaf) const;

        //- Build the binary tree node indexing
        void reindex();

        //- Transfer contents
        void transfer(AABBTreeBase&);


    // Evaluation

        //- Determine whether a point is inside the bounding boxes
        bool pointInside(const point& pt) const;

        //- Determine whether a bounding box overlaps the tree bounding boxes
        bool overlaps(const boundBox& queryBb) const;


    // Write

        //- Write all tree boxes (leaves) in OBJ format
        void writeOBJ(Ostream& os) const;


    // Operators

        //- Copy assignment
        void operator=(const AABBTreeBase&);

        //- Move assignment
        void operator=(AABBTreeBase&&);


    // Housekeeping

        //- The bounding boxes within the tree
        const List<treeBoundBox>& boundBoxes() const { return leafBoxes(); }
};


/*---------------------------------------------------------------------------*\
                          Class AABBTree Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class AABBTree
:
    public AABBTreeBase
{
public:

    // Constructors

        //- Default construct
        constexpr AABBTree() noexcept = default;

        //- Construct with points only (Type is ignored).
        explicit AABBTree
        (
            const pointField& points,
            bool medianSplit = true,                //!< median/geometric split
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct with object list, centroids and vertices
        AABBTree
        (
            //! The list of objects to scan
            const UList<Type>& objects,
            //! The reference centroid for each object
            const pointField& objectCentres,
            //! The referenced points (eg, vertices) for the objects
            const pointField& objectPoints,
            bool medianSplit = true,                //!< median/geometric split
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct with object list and vertices.
        //- The corresponding object centres are estimated on the fly
        AABBTree
        (
            //! The list of objects to scan
            const UList<Type>& objects,
            //! The referenced points (eg, vertices) for the objects
            const pointField& objectPoints,
            bool medianSplit = true,                //!< median/geometric split
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AABBTreeI.H"

#ifdef NoRepository
    #include "AABBTree.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
