/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::treeDataEdge

Description
    Holds data for octree to work on an edges subset.

SourceFiles
    treeDataEdge.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_treeDataEdge_H
#define Foam_treeDataEdge_H

#include "labelRange.H"
#include "line.H"
#include "treeBoundBox.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                           Class treeDataEdge Declaration
\*---------------------------------------------------------------------------*/

class treeDataEdge
{
    // Private Data

        //- Reference to the underlying support points
        const pointField& points_;

        //- Reference to edgeList
        const edgeList& edges_;

        //- Subset of edges to work on
        labelList subset_;

        //- Edge bounding boxes (valid only if cacheBb_)
        treeBoundBoxList bbs_;

        //- Use subset of edges (subset_ active)
        const bool useSubset_;

        //- Whether to precalculate and store edge bounding box
        const bool cacheBb_;


    // Private Member Functions

        //- The object index corresponding to the shape index,
        //- returns -1 for an invalid range
        label mappedIndex(label i) const noexcept
        {
            return
            (
                (i >= 0 && i < size())
              ? (useSubset_ ? subset_.cdata()[i] : i)
              : -1
            );
        }

        //- Initialise all member data
        void update();


public:


    //- Forward to treeDataEdge findNearest operations
    class findNearestOp
    {
        const indexedOctree<treeDataEdge>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataEdge>& tree) noexcept
        :
            tree_(tree)
        {}

        void operator()
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& minIndex,
            point& nearestPoint
        ) const;

        void operator()
        (
            const labelUList& indices,
            const linePointRef& ln,

            treeBoundBox& tightest,
            label& minIndex,
            point& linePoint,
            point& nearestPoint
        ) const;
    };


    //- Forward to treeDataEdge findIntersect operations
    class findIntersectOp
    {
    public:

        findIntersectOp(const indexedOctree<treeDataEdge>&) noexcept {}

        //- Not implemented
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Declare name of the class and its debug switch
    ClassNameNoDebug("treeDataEdge");


    // Generated Methods

        //- Move construct
        treeDataEdge(treeDataEdge&&) = default;


    // Constructors (cachable versions)

        //- Construct from all edges.
        //  \note Holds references to edges and points!
        treeDataEdge
        (
            const bool cacheBb,
            const edgeList& edges,
            const pointField& points
        );

        //- Construct from range of edges.
        //  \note Holds references to edges and points!
        treeDataEdge
        (
            const bool cacheBb,
            const edgeList& edges,
            const pointField& points,
            const labelRange& range
        );

        //- Construct from selected edges.
        //  \note Holds references to edges and points!
        treeDataEdge
        (
            const bool cacheBb,
            const edgeList& edges,
            const pointField& points,
            const labelUList& edgeLabels
        );

        //- Construct from selected edges, transferring contents.
        //  \note Holds references to edges and points!
        treeDataEdge
        (
            const bool cacheBb,
            const edgeList& edges,
            const pointField& points,
            labelList&& edgeLabels
        );


    // Constructors (non-caching)

        //- Construct from all edges.
        //  \note Holds references to edges and points!
        treeDataEdge(const edgeList& edges, const pointField& points)
        :
            treeDataEdge(false, edges, points)
        {}

        //- Construct from range of edges.
        //  \note Holds references to edges and points!
        treeDataEdge
        (
            const edgeList& edges,
            const pointField& points,
            const labelRange& range
        )
        :
            treeDataEdge(false, edges, points, range)
        {}

        //- Construct from selected edges.
        //  \note Holds references to edges and points!
        treeDataEdge
        (
            const edgeList& edges,
            const pointField& points,
            const labelUList& edgeLabels
        )
        :
            treeDataEdge(false, edges, points, edgeLabels)
        {}

        //- Construct from selected edges, transferring contents.
        //  \note Holds references to edges and points!
        treeDataEdge
        (
            const edgeList& edges,
            const pointField& points,
            labelList&& edgeLabels
        )
        :
            treeDataEdge(false, edges, points, std::move(edgeLabels))
        {}


    // Static Functions

        //- Calculate and return bounding boxes for all edges
        static treeBoundBoxList boxes
        (
            const edgeList& edges,
            const pointField& points
        );

        //- Calculate and return bounding boxes for specified range of edges
        static treeBoundBoxList boxes
        (
            const edgeList& edges,
            const pointField& points,
            const labelRange& range
        );

        //- Calculate and return bounding boxes for specified edges
        static treeBoundBoxList boxes
        (
            const edgeList& edges,
            const pointField& points,
            const labelUList& edgeIds
        );

        //- Return bounding box of all edges
        static treeBoundBox bounds
        (
            const edgeList& edges,
            const pointField& points
        );

        //- Return bounding box of specified range of edges
        static treeBoundBox bounds
        (
            const edgeList& edges,
            const pointField& points,
            const labelRange& range
        );

        //- Return bounding box of specified edges
        static treeBoundBox bounds
        (
            const edgeList& edges,
            const pointField& points,
            const labelUList& edgeIds
        );


    // Member Functions

        //- Object dimension == 1 (line element)
        static constexpr int nDim() noexcept { return 1; }

        //- Return bounding box for the specified edge indices
        treeBoundBox bounds(const labelUList& indices) const;


    // Access

        //- The reference point field
        const pointField& points() const noexcept { return points_; }

        //- The original list of edges
        const edgeList& edges() const noexcept { return edges_; }

        //- The subset of edge ids to use
        const labelUList& objectIds() const noexcept { return subset_; }

        //- Use a subset of edges
        bool useSubset() const noexcept { return useSubset_; }

        //- Is the effective edge selection empty?
        bool empty() const noexcept
        {
            return useSubset_ ? subset_.empty() : edges_.empty();
        }

        //- The size of edge selection
        label size() const noexcept
        {
            return useSubset_ ? subset_.size() : edges_.size();
        }

        //- Map from shape index to original (non-subset) edge label
        label objectIndex(const label index) const
        {
            return useSubset_ && index >= 0 ? subset_[index] : index;
        }

        //- Edge at specified shape index
        const edge& operator[](const label index) const
        {
            return edges_[objectIndex(index)];
        }

        //- Geometric line for edge at specified shape index. Frequently used
        const linePointRef line(const label index) const
        {
            return edges_[objectIndex(index)].line(points_);
        }

        //- Representative point (edge centre) at shape index
        point centre(const label index) const
        {
            return edges_[objectIndex(index)].centre(points_);
        }

        //- Representative point cloud for contained shapes.
        //- One point per shape, corresponding to the edge centres.
        tmp<pointField> centres() const;


    // Search

        //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
        //  Only makes sense for closed surfaces.
        volumeType getVolumeType
        (
            const indexedOctree<treeDataEdge>&,
            const point&
        ) const;

        //- Does (bb of) shape at index searchBox
        bool overlaps
        (
            const label index,
            const treeBoundBox& searchBox
        ) const;

        //- Does (bb of) shape at index overlap bb
        bool overlaps
        (
            const label index,
            const point& centre,
            const scalar radiusSqr
        ) const;

        //- Calculates nearest (to sample) point in shape.
        //  Returns actual point and distance (squared)
        void findNearest
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& nearestIndex,
            point& nearestPoint
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
