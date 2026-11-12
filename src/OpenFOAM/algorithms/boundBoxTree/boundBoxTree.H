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

Class
    Foam::boundBoxTree

Description
    A simple binary bounding volume hierarchy (BVH) built directly from an
    existing list of bounding boxes.

    The tree is built top-down. At each node, the entries are split along
    the longest axis of the node bounding box using a median split on the
    entry centroids. Splitting stops when the number of entries falls to
    or below the specified maximum leaf size, or when a split fails to
    reduce the number of entries in a child (e.g. many coincident boxes).

    Once built, the tree supports rapid broad-phase queries.

    Example use\n
    \verbatim
        treeBoundBoxList boxes(...);
        boundBoxTree tree(std::move(boxes));

        i = tree.first_overlaps(queryBox);
    \endverbatim

SourceFiles
    boundBoxTree.cxx
    boundBoxTree.txx
    boundBoxTreeI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundBoxTree_H
#define Foam_boundBoxTree_H

#include "treeBoundBox.H"
#include "boundBoxNodes.H"
#include "DynamicList.H"
#include "InfoProxy.H"
#include "Tuple2.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class boundBoxTree;
Ostream& operator<<(Ostream&, const boundBoxTree&);
Ostream& operator<<(Ostream&, const InfoProxy<boundBoxTree>&);

/*---------------------------------------------------------------------------*\
                        Class boundBoxTree Declaration
\*---------------------------------------------------------------------------*/

class boundBoxTree
{
public:

    // Public Data Types

        //- Different splitting approaches
        enum class algorithm
        {
            //! split longest bounding box axis of the box centres
            median_split,
            //! split longest bounding box axis into equal geometric sizes
            geometric_split,
            //! split along the longest axis of the range bounding box
            box_split,
            //! split on a pre-balanced set of spatially grouped boxes
            pre_balanced
        };

        //- The type for bounding box binary tree nodes
        using node = meshTools::BoundBoxBinaryNode<treeBoundBox>;

        //- Reference to an empty leaf (placeholder)
        static const node& empty_leaf();


    // Properties

        //- Default splitting approach
        static constexpr algorithm defaultAlgorithm() noexcept
        {
            return algorithm::median_split;
        }

        //- The default max tree depth/level
        static constexpr int defaultMaxLevel() noexcept { return 8; }

        //- The default min leaf size
        static constexpr int defaultMinLeafSize() noexcept { return 8; }

        //- Index of the root node within nodes_ (always 0)
        static constexpr int root() noexcept { return 0; }


private:

    // Private Data Types

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

    // Control Parameters

        //- Target occupancy (for auto-size mode)
        static constexpr int targetOccupancy = 64;

        //- Max depth level (for auto-size mode)
        static constexpr int maxDepthLimit = 12;


    // Private Member Functions

        //- Make reasonable guesses for level and occupancy limits
        static void autosize
        (
            control& ctrl,          //!< [in,out] control parameters
            label total,            //!< total number of objects
            int maxLevel = -1,      //!< requested max level
            int minLeafSize = -1    //!< requested min leaf size
        );


    // Private Data Types

        //- Predicate: contains the specified box
        struct containsBoxOp
        {
            const boundBox& box;

            bool operator()(const treeBoundBox& node_bb) const
            {
                return node_bb.contains(box);
            }
        };

        //- Predicate: overlaps the specified box
        struct overlapsBoxOp
        {
            const boundBox& box;

            bool operator()(const treeBoundBox& node_bb) const
            {
                return node_bb.overlaps(box);
            }
        };

        //- Predicate: contains the specified point
        struct containsPointOp
        {
            const point& p;

            bool operator()(const treeBoundBox& node_bb) const
            {
                return node_bb.contains(p);
            }
        };

        //- Predicate: contains the specified bounding sphere
        struct containsSphereOp
        {
            const point& centre;
            const scalar radiusSqr;

            bool operator()(const treeBoundBox& node_bb) const
            {
                return node_bb.contains(centre, radiusSqr);
            }
        };

        //- Predicate: overlaps the specified bounding sphere
        struct overlapsSphereOp
        {
            const point& centre;
            const scalar radiusSqr;

            bool operator()(const treeBoundBox& node_bb) const
            {
                return node_bb.overlaps(centre, radiusSqr);
            }
        };

        //- Predicate: contains both ends of a line segment
        struct containsLineOp
        {
            const point& from;
            const point& to;

            bool operator()(const treeBoundBox& node_bb) const
            {
                return node_bb.contains_all(from, to);
            }
        };

        //- Predicate: intersected by a line segment
        struct intersectsLineOp
        {
            const point& from;
            const point& to;

            bool operator()(const treeBoundBox& node_bb) const
            {
                return node_bb.intersects(from, to);
            }
        };


    // Private Data

        //- Original (unordered) bounding boxes, as supplied on construction
        List<treeBoundBox> leafBoxes_;

        //- Permutation of original box indices, grouped by leaf.
        //  A leaf node covers the indices_ defined by start(),size() range
        labelList indices_;

        //- All tree nodes. nodes_[0] is the root (empty if no boxes)
        List<node> nodes_;

        //- Control parameters
        control ctrl_;


    // Private Member Functions

        //- Assign new min leaf size, with value clamping
        void setAlgorithm(algorithm algo) noexcept
        {
            ctrl_.algo = algo;
        }

        //- Assign new max level, with value clamping
        void setMaxLevel(int n) noexcept
        {
            if (n > 0) ctrl_.maxLevel = n;
        }

        //- Assign new min leaf size, with value clamping
        void setMinLeafSize(int n) noexcept
        {
            if (n > 0) ctrl_.minLeafSize = n;
        }

        //- Number of elements at or below given node
        label countElements(label nodei) const;

        //- Number of leafs at or below given node
        label countLeafs(label nodei) const;

        //- Write OBJ for node + sub-nodes
        void writeOBJ
        (
            const label nodeId,
            Ostream& os,
            label& vertIndex,
            bool leavesOnly,
            bool writeLinesOnly
        ) const;

        //- Recursively build the tree over indices_[first, last).
        //
        //  algorithm::box_split
        //
        //  * splits along the longest axis of the range bounding box
        //    using a median split on the entry centroids, reordering
        //    indices_ as needed (general-purpose).
        //
        //  algorithm::centroid_split
        //
        //  * like box_split, but only examining the bounding box centres
        //
        //  algorithm::pre_balanced
        //
        //  * trusts that indices_[first, last) are already
        //    spatially grouped into contiguous left/right halves
        //    (e.g. by a prior k-d tree split).
        //    The range is simply bisected
        //
        //  \returns the index of the newly created node.
        label build_impl
        (
            const label first,
            const label last,
            DynamicList<node>& nodesToBuild,
            const int currLevel
        );


    // Queries

        //- Depth-first traversal of the sub-tree rooted at nodeId,
        //- dispatching to the given visitor. See the public visit()
        //- for the required Visitor interface.
        template<class Visitor>
        void visit(label nodei, Visitor& visitor) const;

        //- Depth-first traversal of the sub-tree rooted at nodeId,
        //- dispatching to the given visitor. See the public visit()
        //- for the required Visitor interface.
        template<class Visitor>
        void visit_leaf(label nodei, Visitor& visitor) const;

        //- Depth-first walking, returning first leaf object index
        //- matching the predicate (-1 for no matches).
        template<class Predicate>
        label find_first_of(label nodei, const Predicate& pred) const;

        //- Depth-first walking, counting leaf items matches
        //- The predicate expects a treeBoundBox as argument
        template<class Predicate>
        label count_if(label nodei, const Predicate& pred) const;

        //- Depth-first search collecting leaf indices for predicate matches.
        //- The predicate expects a treeBoundBox as argument
        template<class Predicate>
        void find_all
        (
            label nodei,
            const Predicate& pred,
            DynamicList<label>& hits
        ) const;


public:

    // Generated Methods

        //- Default construct
        boundBoxTree() noexcept = default;


    // Constructors

        //- Construct empty with specified controls
        inline explicit boundBoxTree
        (
            boundBoxTree::algorithm algo,
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        ) noexcept;

        //- Construct from a list of bounding boxes,
        //- building the tree with the specified algorithm
        inline boundBoxTree
        (
            boundBoxTree::algorithm algo,
            List<treeBoundBox>&& boxes,
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct from a list of bounding boxes,
        //- building the tree with the specified algorithm
        inline explicit boundBoxTree
        (
            boundBoxTree::algorithm algo,
            const UList<treeBoundBox>& boxes,
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct from a list of bounding boxes,
        //- building the tree with the specified algorithm
        inline explicit boundBoxTree
        (
            boundBoxTree::algorithm algo,
            const UList<boundBox>& boxes,
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct from a list of bounding boxes,
        //- building the tree with the default algorithm
        inline explicit boundBoxTree
        (
            List<treeBoundBox>&& boxes,
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct from a list of bounding boxes,
        //- building the tree with the default algorithm
        inline explicit boundBoxTree
        (
            const UList<treeBoundBox>& boxes,
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );

        //- Construct from a list of bounding boxes,
        //- building the tree with the default algorithm
        inline explicit boundBoxTree
        (
            const UList<boundBox>& boxes,
            int maxLevel = defaultMaxLevel(),       //!< max tree levels
            int minLeafSize = defaultMinLeafSize()  //!< min elements per leaf
        );


    // Member Functions

        //- Discard the tree and list of bounding boxes.
        void clear();

        //- Discard the tree addressing but not the bounding boxes
        void clearAddressing();

        //- Reset the list of bounding boxes, clearing the tree addressing
        void setBoxes(List<treeBoundBox>&& boxes);

        //- Reset the list of bounding boxes, clearing the tree addressing
        void setBoxes(const UList<treeBoundBox>& boxes);

        //- Reset the list of bounding boxes, clearing the tree addressing
        void setBoxes(const UList<boundBox>& boxes);

        //- Rebuild the tree from the existing set of bounding boxes
        void rebuild();

        //- True if the tree contains built nodes
        bool good() const noexcept { return !nodes_.empty(); }

        //- True if the tree contains no boxes
        bool empty() const noexcept { return leafBoxes_.empty(); }

        //- Number of boxes used to build the tree
        label size() const noexcept { return leafBoxes_.size(); }

        //- True if the tree contains built nodes
        bool hasNodes() const noexcept { return !nodes_.empty(); }

        //- Number of nodes (interior and data leafs)
        label nNodes() const noexcept { return nodes_.size(); }

        //- Number of interior nodes
        label nInterior() const;

        //- Number of data leafs
        label nLeafs() const;

        //- The original (unordered) bounding boxes
        const List<treeBoundBox>& boxes() const noexcept { return leafBoxes_; }

        //- Discard any tree addressing and
        //- return the original (unordered) bounding boxes
        List<treeBoundBox> release_boxes()
        {
            clearAddressing();
            return std::move(leafBoxes_);
        }

        //- The permutation of box indices, grouped by leaf
        const labelList& indices() const noexcept { return indices_; }

        //- All tree nodes (nodes_[root()] is the root)
        const List<node>& nodes() const noexcept { return nodes_; }

        //- The ranges into indices() for each leaf,
        //- constructed from the node information.
        List<labelRange> ranges() const;

        //- The overall bounding box (root node), or an inverted box
        //- if the tree is empty
        inline const treeBoundBox& bb() const;


    // Member Operators

        //- Return the original (unordered) bounding box for the given leaf
        inline const treeBoundBox& operator[](label leafi) const;


    // Information (fragile API)

       int maxLevel() const noexcept { return ctrl_.maxLevel; }
       int minLeafSize() const noexcept { return ctrl_.minLeafSize; }


    // Write

        //- Write all tree boxes (leaves) in OBJ format
        void writeOBJ(Ostream& os) const;

        //- Return info proxy for printing information to a stream
        InfoProxy<boundBoxTree> info() const noexcept { return *this; }


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const boundBoxTree&);
        friend Ostream& operator<<(Ostream&, const InfoProxy<boundBoxTree>&);


    // General Queries

        //- True if the predicate matches on any leaf object.
        //- Same as (find_first_of(pred) >= 0). Uses depth-first-search.
        template<class Predicate>
        bool any_of(const Predicate& pred) const;

        //- Return the first leaf object index for which the predicate matches
        //- (-1 for no matches). Uses depth-first-search.
        template<class Predicate>
        label find_first_of(const Predicate& pred) const;

        //- Return the number of leaf objects matching the predicate.
        //- Uses depth-first-search.
        template<class Predicate>
        label count_if(const Predicate& pred) const;

        //- Collect leaf object indices matching the predicate.
        //- Uses depth-first-search.
        template<class Predicate>
        bool find_all
        (
            const Predicate& pred,
            //! [out] The matching object indices
            DynamicList<label>& hits
        ) const;


    // Queries

        //- Index of first indexed box that overlaps \p queryBb
        inline label first_overlaps(const boundBox& queryBb) const;

        //- True if any indexed box overlaps \p queryBb
        bool overlaps(const boundBox& queryBb) const
        {
            return (first_overlaps(queryBb) >= 0);
        }

        //- Indices of the boxes that overlap \p queryBb
        inline bool find_overlaps
        (
            const boundBox& queryBb,
            //! [out] indices of boxes that overlap
            DynamicList<label>& hits
        ) const;

        //- Indices of the boxes that overlap \p queryBb
        //- returned as a range (the supplied storage) for range-for use.
        inline const labelUList& all_overlaps
        (
            const boundBox& queryBb,
            //! [out] indices of boxes that overlap
            DynamicList<label>& hits
        ) const;

        //- Index of first indexed box that overlaps the bounding sphere
        inline label first_overlaps
        (
            const point& centre,
            scalar radiusSqr
        ) const;

        //- True if any indexed box overlaps the bounding sphere
        bool overlaps(const point& centre, scalar radiusSqr) const
        {
            return (first_overlaps(centre, radiusSqr) >= 0);
        }

        //- Indices of the boxes that overlap the bounding sphere
        inline bool find_overlaps
        (
            const point& centre,
            scalar radiusSqr,
            //! [out] indices of boxes that overlap
            DynamicList<label>& hits
        ) const;

        //- Indices of the boxes that overlap the bounding sphere
        //- returned as a range (the supplied storage) for range-for use.
        inline const labelUList& all_overlaps
        (
            const point& centre,
            scalar radiusSqr,
            //! [out] indices of boxes that overlap
            DynamicList<label>& hits
        ) const;

        //- Index of first box that contains point (inside or on edge)
        inline label first_contains(const point& p) const;

        //- Contains point? (inside or on edge)
        bool contains(const point& p) const
        {
            return (first_contains(p) >= 0);
        }

        //- Indices of the boxes that contain point (inside or on edge)
        inline bool find_contains
        (
            const point& p,
            //! [out] indices of boxes that contain point
            DynamicList<label>& hits
        ) const;

        //- Index of first box containing the bounding sphere
        //- (centre + sqr(radius))
        inline label first_contains
        (
            const point& centre,
            scalar radiusSqr
        ) const;

        //- Contains the bounding sphere (centre + sqr(radius))
        inline bool contains(const point& centre, scalar radiusSqr) const;

        //- Index of first box that contains both line ends
        label first_contains(const point& p0, const point& p1) const;

        //- Contains both line ends?
        bool contains(const point& p0, const point& p1) const
        {
            return (first_contains(p0, p1) >= 0);
        }

        //- Indices of the boxes that contain both line ends
        inline bool find_contains
        (
            const point& p0,
            const point& p1,
            //! [out] indices of boxes that contain both line ends
            DynamicList<label>& hits
        ) const;

        //- Indices of the boxes that contain both line ends,
        //- returned as a range (the supplied storage) for range-for use.
        inline const labelUList& all_contains
        (
            const point& p0,
            const point& p1,
            //! [out] indices of boxes that contain both line ends
            DynamicList<label>& hits
        ) const;

        //- Index of first box that contains other boundingBox
        label first_contains(const boundBox& other) const;

        //- Fully contains other boundingBox?
        bool contains(const boundBox& other) const
        {
            return (first_contains(other) >= 0);
        }

        //- Indices of the boxes that fully contain other boundingBox
        inline bool find_contains
        (
            const boundBox& other,
            //! [out] indices of boxes that fully contain \p other
            DynamicList<label>& hits
        ) const;


        //- Index of first box that the from/to ray intersects
        inline label first_intersects(const point& from, const point& to) const;

        //- Index of first box that the from/to ray intersects
        inline label first_intersects(const linePointRef& ln) const;

        //- True if the from/to ray intersects any box
        bool intersects(const point& from, const point& to) const
        {
            return (first_intersects(from, to) >= 0);
        }

        //- True if the from/to ray intersects any box
        bool intersects(const linePointRef& ln) const
        {
            return (first_intersects(ln) >= 0);
        }

        //- Indices of the boxes that intersect from/to ray
        inline bool find_intersects
        (
            const point& from,
            const point& to,
            //! [out] indices of boxes that intersect from/to ray
            DynamicList<label>& hits
        ) const;

        //- Indices of the boxes that intersect from/to ray
        inline bool find_intersects
        (
            const linePointRef& ln,
            //! [out] indices of boxes that intersect from/to ray
            DynamicList<label>& hits
        ) const;

        //- Indices of the boxes that intersect from/to ray,
        //- returned as a range (the supplied storage) for range-for use.
        inline const labelUList& all_intersects
        (
            const point& from,
            const point& to,
            //! [out] indices of boxes that intersect from/to ray
            DynamicList<label>& hits
        ) const;

        //- Indices of the boxes that intersect from/to ray,
        //- returned as a range (the supplied storage) for range-for use.
        inline const labelUList& all_intersects
        (
            const linePointRef& ln,
            //! [out] indices of boxes that intersect from/to ray
            DynamicList<label>& hits
        ) const;


    // More search functions

        //- Find leaf boxes that are inside or overlap the bounding sphere
        label findSphere
        (
            const point& centre,
            const scalar radiusSqr,
            DynamicList<label>& hits
        ) const;

        //- Find the leaf box nearest to centre, plus all leaf boxes
        //- overlapping the tightened search bound (hits).
        //
        //  Returns (boundSqr, volumeType):
        //  - boundSqr : tightened search bound (radius squared) guaranteed
        //    to reach a hit, never larger than radiusSqr.
        //    - for INSIDE: the smallest sphere enclosing a containing box;
        //    - for MIXED: the far distance of the nearest box
        //
        //  - INSIDE  : centre is inside the boxes in hits
        //  - MIXED   : no containing box; nearest overlap only
        //  - OUTSIDE : search completed, nothing within bound (hits empty)
        //  - UNKNOWN : empty tree (hits empty)
        Tuple2<scalar, volumeType>
        findNearest
        (
            //! The bounding sphere centre
            const point& centre,
            //! The initial search bound (radius squared)
            const scalar radiusSqr,
            //! [out] indices of the containing (INSIDE) or overlapping (MIXED)
            //! boxes. Sorted \em nearest \em first. Empty for OUTSIDE, UNKNOWN.
            DynamicList<label>& hits
        ) const;


    // Iteration

        //- Depth-first tree traversal, dispatching to given visitor.
        //
        //  The Visitor type must provide:
        //  \code
        //      // Called for each leaf node.
        //      // indices addresses the (reordered) object indices
        //      // within the leaf.
        //      void leaf(const treeBoundBox&, const labelUList& indices);
        //
        //      // Called for each interior node, before descending into
        //      // its children. Return true to continue descending into
        //      // the sub-tree, false to prune it.
        //      // (e.g. after a bounding box rejection test, or once a
        //      // visitor-specific search has already been satisfied).
        //      bool interior(const treeBoundBox&);
        //  \endcode
        template<class Visitor>
        void visit(Visitor& visitor) const;

        //- Depth-first tree traversal of leafs, dispatching to given visitor.
        //
        //  The Visitor type must provide a caller:
        //  \code
        //      void operator()(const node& nd);
        //  \endcode
        template<class Visitor>
        void visit_leaf(Visitor&& visitor) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundBoxTreeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "boundBoxTree.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
