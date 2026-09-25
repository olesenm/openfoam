/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::treeDataPoint

Description
    Holds (reference to) pointField. Encapsulation of data needed for
    octree searches.
    Used for searching for nearest point. No bounding boxes around points.
    Only overlaps and calcNearest are implemented, rest makes little sense.

    Optionally works on subset of points.

SourceFiles
    treeDataPoint.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_treeDataPoint_H
#define Foam_treeDataPoint_H

#include "pointField.H"
#include "line.H"
#include "treeBoundBox.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                        Class treeDataPoint Declaration
\*---------------------------------------------------------------------------*/

class treeDataPoint
{
    // Private Data

        //- Reference to the underlying point field
        const pointField& points_;

        //- Subset of points to work on (or empty)
        labelList subset_;

        //- Use subset of points (subset_ active)
        const bool useSubset_;


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

public:

    //- Forward to treeDataPoint findNearest operations
    class findNearestOp
    {
        const indexedOctree<treeDataPoint>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataPoint>& tree) noexcept
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


    //- Forward to treeDataPoint findIntersect operations (not possible)
    class findIntersectOp
    {
    public:

        findIntersectOp(const indexedOctree<treeDataPoint>& tree) noexcept {}

        //- Not implemented
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Declare name of the class
    ClassNameNoDebug("treeDataPoint");


    // Generated Methods

        //- Move construct
        treeDataPoint(treeDataPoint&&) = default;


    // Constructors (non-caching)

        //- Construct from pointField
        //  \note Holds reference to the points!
        explicit treeDataPoint(const pointField& points) noexcept;

        //- Construct from subset of pointField, copies point ids
        //  \note Holds reference to the points!
        treeDataPoint
        (
            const pointField& points,
            const labelUList& pointLabels,
            //! optional override of subset
            //! (eg, allow empty pointLabels to mean no subset)
            const bool withSubset = true
        );

        //- Construct from subset of pointField, moves point ids
        //  \note Holds reference to the points!
        treeDataPoint
        (
            const pointField& points,
            labelList&& pointLabels,
            //! optional override of subset
            //! (eg, allow empty pointLabels to mean no subset)
            const bool withSubset = true
        );


    // Member Functions

        //- Object dimension == 0 (point element)
        static constexpr int nDim() noexcept { return 0; }

        //- Return bounding box for the specified point indices
        treeBoundBox bounds(const labelUList& indices) const;


    // Access

        //- The original point field
        const pointField& points() const noexcept { return points_; }

        //- The subset of point ids being used
        const labelUList& objectIds() const noexcept { return subset_; }

        //- Use a subset of points
        bool useSubset() const noexcept { return useSubset_; }

        //- Is the effective point field empty?
        bool empty() const noexcept
        {
            return useSubset_ ? subset_.empty() : points_.empty();
        }

        //- The size of the effective point field
        label size() const noexcept
        {
            return useSubset_ ? subset_.size() : points_.size();
        }

        //- Map to the original (non-subset) point label
        label objectIndex(const label index) const
        {
            return useSubset_ && index >= 0 ? subset_[index] : index;
        }

        //- Point at specified shape index
        const point& operator[](const label index) const
        {
            return points_[objectIndex(index)];
        }

        //- Point at specified shape index
        const point& centre(const label index) const
        {
            return points_[objectIndex(index)];
        }

        //- Point cloud
        tmp<pointField> centres() const;


    // Search

        //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
        //  Only makes sense for closed surfaces.
        volumeType getVolumeType
        (
            const indexedOctree<treeDataPoint>& os,
            const point& sample
        ) const;

        //- Does (bb of) shape at index searchBox
        bool overlaps
        (
            const label index,
            const treeBoundBox& searchBox
        ) const;

        //- Does shape at index overlap the sphere
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


    // Housekeeping

        //- Map to the original (non-subset) point label
        ///FOAM_DEPRECATED_FOR(2022-10, "objectIndex()")
        label pointLabel(label index) const { return objectIndex(index); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
