/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
    Foam::boundBox

Description
    A bounding box defined in terms of min/max extrema points.

Note
    When a bounding box is created without any points, it creates an inverted
    bounding box. Points can be added later and the bounding box will grow to
    include them.

    \verbatim
    Octant vertex ordering               | Hex cell ordering
    (treeBoundBox)                       | (boundBox, blockMesh, hexCell)
                                         |
                    6 ---- 7             |         7 ---- 6
              f5    |\     :\    f3      |   f5    |\     :\    f3
               |    | 4 ---- 5    \      |    |    | 4 ---- 5    \
               |    2.|....3 |     \     |    |    3.|....2 |     \
               |     \|     \|     f2    |    |     \|     \|     f2
              f4      0 ---- 1           |   f4      0 ---- 1
      Y  Z                               |
       \ |          f0 ------ f1         |         f0 ------ f1
        \|                               |
         o--- X                          |
     \endverbatim

SeeAlso
    Foam::treeBoundBox

SourceFiles
    boundBoxI.H
    boundBox.cxx
    boundBox.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundBox_H
#define Foam_boundBox_H

#include "pointField.H"
#include "faceList.H"
#include "Pair.H"
#include "triangleFwd.H"

// Transitional features/misfeatures.
// Default is syncPar = true, an explicit parameter under "+strict"
#ifdef FOAM_COMPILE_STRICT
#  define Foam_boundBox_syncPar(name) name
#else
#  define Foam_boundBox_syncPar(name) name = true
#endif

// No default syncPar handling (an ideal future...)
// #undef  Foam_boundBox_syncPar
// #define Foam_boundBox_syncPar(name) name

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class boundBox;
class treeBoundBox;
class plane;
class Random;
template<class T> class MinMax;

Istream& operator>>(Istream& is, boundBox& bb);
Ostream& operator<<(Ostream& os, const boundBox& bb);

/*---------------------------------------------------------------------------*\
                          Class boundBox Declaration
\*---------------------------------------------------------------------------*/

class boundBox
{
    // Private Data

        //- Minimum and maximum points describing the bounding box
        point min_, max_;


protected:

    // Protected Member Functions

        //- Test for overlap of box and box (inclusive check)
        inline static bool box_box_overlaps
        (
            const point& minA,  // boxA (min)
            const point& maxA,  // boxA (max)
            const point& minB,  // boxB (min)
            const point& maxB   // boxB (max)
        );

        //- Test for overlap of box and bounding sphere (centre + sqr(radius))
        //  Note: ordering of corners is irrelevant
        inline static bool box_sphere_overlaps
        (
            const point& corner0,   // box corner
            const point& corner1,   // box corner
            const point& centre,    // sphere centre
            const scalar radiusSqr  // sqr(radius)
        );

        //- Test if box contains the bounding sphere (centre + sqr(radius))
        inline static bool box_sphere_contains
        (
            const point& min_corner,    // box corner (min)
            const point& max_corner,    // box corner (max)
            const point& centre,        // sphere centre
            const scalar radiusSqr      // sqr(radius)
        );


public:

    // Static Data Members

        //- Bits used for (x/y/z) directional encoding.
        enum directionBit : direction
        {
            XDIR = 1,  //!< 1: x-direction. Same as (1 << vector::X)
            YDIR = 2,  //!< 2: y-direction. Same as (1 << vector::Y)
            ZDIR = 4   //!< 4: z-direction. Same as (1 << vector::Z)
        };

        //- A large boundBox: min/max == -/+ ROOTVGREAT
        static const boundBox greatBox;

        //- A large inverted boundBox: min/max == +/- ROOTVGREAT
        static const boundBox invertedBox;

        //- The unit normal per face
        static const FixedList<vector, 6> faceNormals;


    // Properties

        //- Number of points for boundBox and HEX
        static constexpr label nPoints() noexcept { return 8; }

        //- Number of edges for boundBox and HEX
        static constexpr label nEdges() noexcept { return 12; }

        //- Number of faces for boundBox and HEX
        static constexpr label nFaces() noexcept { return 6; }


    // Static Methods

        //- The null boundBox is the same as an inverted box
        static const boundBox& null() noexcept { return invertedBox; }

        //- The boundBox faces as a hexCell, using hexCorner points.
        //- Same as hexCell::modelFaces()
        static const Foam::faceList& hexFaces();

        //- The boundBox faces using octCorner points.
        static const Foam::faceList& octFaces();

        //- The boundBox edges, using octCorner points.
        static const Foam::edgeList& octEdges();


    // Standard (Generated) Methods

        //- Default construct: an inverted bounding box
        inline boundBox();

        //- Copy construct
        boundBox(const boundBox&) = default;

        //- Copy assignment
        boundBox& operator=(const boundBox&) = default;


    // Constructors

        //- Copy construct with specified global reduction
        boundBox(const boundBox& bb, bool syncPar);

        //- Construct a bounding box containing a single initial point
        inline explicit boundBox(const point& p);

        //- Construct a 0/1 unit bounding box
        inline explicit boundBox(Foam::zero_one);

        //- Construct from bound box min/max points
        inline boundBox(const point& min, const point& max);

        ///TBD: Construct from bound box min/max points
        ///inline explicit boundBox(const MinMax<point>& bb);

        //- Construct from bound box min/max points
        inline explicit boundBox(const Pair<point>& bb);

        //- Construct with extent of two bounding boxes, no reduction
        inline boundBox(const boundBox& a, const boundBox& b);

        //- Construct as the bounding box of the given points
        //  Does parallel communication (syncPar = true)
        explicit boundBox
        (
            const UList<point>& points,
            //! Do parallel reduction
            bool Foam_boundBox_syncPar(syncPar)
        );

        //- Construct as the bounding box of the given temporary pointField.
        //  Does parallel communication (syncPar = true)
        explicit boundBox
        (
            const tmp<pointField>& tpoints,
            //! Do parallel reduction
            bool Foam_boundBox_syncPar(syncPar)
        );

        //- Construct bounding box as an indirect subset of the points.
        //  The indices could be from cell/face etc.
        //  Does parallel communication (syncPar = true)
        boundBox
        (
            const UList<point>& points,
            const labelUList& indices,
            //! Do parallel reduction
            bool Foam_boundBox_syncPar(syncPar)
        );

        //- Construct bounding box as an indirect subset of the points.
        //  The indices could be from edge/triFace etc.
        //  Does parallel communication (syncPar = true)
        template<unsigned N>
        boundBox
        (
            const UList<point>& points,
            const FixedList<label, N>& indices,
            //! Do parallel reduction
            bool Foam_boundBox_syncPar(syncPar)
        );

        //- Construct from Istream
        inline explicit boundBox(Istream& is);


    // Member Functions

        // Access

            //- Bounding box is inverted, contains no points.
            inline bool empty() const noexcept;

            //- Bounding box is non-inverted.
            inline bool good() const noexcept;

            //- Bounding box is non-inverted - same as good().
            bool valid() const { return good(); }

            //- Minimum describing the bounding box
            inline const point& min() const noexcept;

            //- Maximum describing the bounding box
            inline const point& max() const noexcept;

            //- Minimum describing the bounding box, non-const access
            inline point& min() noexcept;

            //- Maximum describing the bounding box, non-const access
            inline point& max() noexcept;

            //- The centre (midpoint) of the bounding box
            inline point centre() const;

            //- Component of centre (midpoint) of the bounding box
            inline scalar centre(direction cmpt) const;

            //- Linear interpolation between bounding box min/max
            inline point lerp(scalar t) const;

            //- Linear interpolation between bounding box min/max
            inline point lerp(const vector& t) const;

            //- The bounding box span (from minimum to maximum)
            inline vector span() const;

            //- Component of bounding box span (from minimum to maximum)
            inline scalar span(direction cmpt) const;

            //- The magnitude/length of the bounding box diagonal
            inline scalar mag() const;

            //- The magnitude/length squared of bounding box diagonal
            inline scalar magSqr() const;

            //- The volume of the bound box
            inline scalar volume() const;

            //- Smallest length/height/width dimension
            inline scalar minDim() const;

            //- Largest length/height/width dimension
            inline scalar maxDim() const;

            //- Average length/height/width dimension
            inline scalar avgDim() const;

            //- Direction (X/Y/Z) of the smallest span (for empty box: 0)
            inline direction minDir() const;

            //- Direction (X/Y/Z) of the largest span (for empty box: 0)
            inline direction maxDir() const;

            //- Count the number of positive, non-zero dimensions.
            //  \return
            //   - -1 : if any dimensions are negative
            //   -  0 : 0D (point)
            //   -  1 : 1D (line aligned with an axis)
            //   -  2 : 2D (plane aligned with an axis)
            //   -  3 : 3D (box)
            inline int nDim() const;


        // Geometry

            //- Return corner point [0..7] corresponding to a 'hex' cell
            template<direction I>
            inline point hexCorner() const;

            //- Fill in 8 corner points in an order corresponding
            //- to a 'hex' cell
            void hexCorners(UList<point>& output) const;

            //- Corner points in an order corresponding to a 'hex' cell
            pointField hexCorners() const;

            //- Return corner point of given octant [0..7]
            template<direction Octant>
            inline point octCorner() const;

            //- Fill in 8 corner points in octant ordering
            void octCorners(UList<point>& output) const;

            //- Corner points in octant coding
            pointField octCorners() const;

            //- Return corner point of given octant [0..7].
            inline point octCorner(direction octant) const;

            //- Return corner point of given octant [0..7] == octCorner
            point corner(direction octant) const { return octCorner(octant); }

            //- Corner points in an order corresponding to a 'hex' cell
            tmp<pointField> points() const
            {
                return tmp<pointField>::New(hexCorners());
            }

            //- Fill in 6 face midpoints (hex/oct order are identical)
            void faceCentres(UList<point>& output) const;

            //- Face midpoints (hex/oct order are identical)
            pointField faceCentres() const;

            //- Face centre of given face index (hex/oct order are identical)
            point faceCentre(direction facei) const;


        // Manipulate

            //- Reset to an inverted box
            inline void reset();

            //- Reset to a 0/1 unit bounding box
            inline void reset(Foam::zero_one);

            //- Reset min/max to be identical to the specified point
            inline void reset(const point& pt);

            //- Reset min/max to specified values
            inline void reset(const point& min, const point& max);

            //- Same as reset() - reset to an inverted box
            void clear() { reset(); }

            //- Extend to include the second box.
            inline void add(const boundBox& bb);

            //- Extend to include the second bound box.
            //- Here the point pair are treated semantically like
            //- a bound box min/max pair.
            inline void add(const Pair<point>& bb);

            //- Extend to include the point.
            inline void add(const point& pt);

            //- Extend to include two additional points.
            inline void add(const point& pt0, const point& pt1);

            //- Extend to include the points.
            inline void add(const UList<point>& points);

            //- Extend to include the points from the temporary point field.
            inline void add(const tmp<pointField>& tpoints);

            //- Extend to include the points.
            template<unsigned N>
            void add(const FixedList<point, N>& points);

            //- Extend to include a (subsetted) point field.
            //
            //  \tparam IndexContainer  A container with an iterator that
            //      dereferences to a label
            //      (eg, FixedList, labelUList, face, triFace, ...)
            template<class IndexContainer>
            void add
            (
                const UList<point>& points,
                const IndexContainer& indices
            );

            //- Expand box by adjusting min/max by specified amount
            //- in each dimension
            inline void grow(scalar delta);

            //- Expand box by adjusting min/max by specified amounts
            inline void grow(const vector& delta);

            //- Expand box by factor*mag(span) in all dimensions
            inline void inflate(scalar factor);

            //- Expand box so each component (axis) span grows by
            //- (1+epsilon), scaled about the centre.
            //- This makes this an anisotropic, shape-preserving expansion.
            inline void inflate_components(scalar epsilon);

            //- Expand box slightly by expanding all dimensions with
            //- factor*span*(random 0-1) and guarantees
            //- factor*mag(span) minimum width in any direction
            void inflate(Random& rndGen, scalar factor);

            //- As per two parameter version but with additional
            //- grow() by given amount in each dimension
            void inflate(Random& rndGen, scalar factor, scalar delta);

            //- Expand box slightly by expanding all dimensions with
            //- factor*span*(random 0-1) and guarantees
            //- factor*mag(span) minimum width in any direction
            void inflate_random(label seedValue, scalar factor);

            //- As per two parameter version but with additional
            //- grow() by given amount in each dimension
            void inflate_random(label seedValue, scalar factor, scalar delta);


        // Parallel

            //- Inplace parallel reduction of min/max values.
            //- Uses worldComm (default) or the specified communicator.
            void reduce(int communicator = -1);

            //- Perform a reduction on a copy and return the result.
            //- Uses worldComm (default) or the specified communicator.
            inline boundBox returnReduce(int communicator = -1) const;

            //- Parallel reduction of min/max values for multiple boxes.
            //- Uses worldComm (default) or the specified communicator.
            static void reduce
            (
                boundBox bbs[],
                label count,
                int communicator = -1
            );

            //- Parallel reduction of min/max values for multiple boxes.
            //- Uses worldComm (default) or the specified communicator.
            static void reduce
            (
                treeBoundBox bbs[],
                label count,
                int communicator = -1
            );

            //- Parallel broadcast of contents.
            //- Uses worldComm (default) or the specified communicator.
            void broadcast(int communicator = -1);

            //- Parallel broadcast of contents for multiple boxes.
            //- Uses worldComm (default) or the specified communicator.
            static void broadcast
            (
                boundBox bbs[],
                label count,
                int communicator = -1
            );

            //- Parallel broadcast of contents for multiple boxes.
            //- Uses worldComm (default) or the specified communicator.
            static void broadcast
            (
                treeBoundBox bbs[],
                label count,
                int communicator = -1
            );


        // Query

            //- Does plane intersect this bounding box.
            //  There is an intersection if the plane segments the corner points
            //  \note the results are unreliable when plane coincides almost
            //  exactly with a box face
            bool intersects(const plane& pln) const;

            //- Does triangle intersect this bounding box
            //- or is contained within this bounding box.
            //  \note results may be unreliable when it coincides almost
            //  exactly with a box face
            bool intersects(const triPointRef& tri) const;

            //- Does line segment intersect this bounding box?
            //  Forwards to treeBoundBox implementation.
            bool intersects(const point& from, const point& to) const;

            //- Does line segment intersect this bounding box?
            //  Forwards to treeBoundBox implementation.
            inline bool intersects(const linePointRef& ln) const;

            //- Overlaps/touches boundingBox?
            inline bool overlaps(const boundBox& bb) const;

            //- Overlaps bounding sphere (centre + sqr(radius))?
            inline bool overlaps(const point& centre, scalar radiusSqr) const;

            //- Contains point? (inside or on edge)
            inline bool contains(const point& pt) const noexcept;

            //- Contains point? (inside only)
            inline bool containsInside(const point& pt) const noexcept;

            //- Fully contains other boundingBox?
            inline bool contains(const boundBox& bb) const noexcept;

            //- Fully contains bounding sphere (centre + sqr(radius))
            inline bool contains(const point& centre, scalar radiusSqr) const;

            //- Contains both points? (eg, a line)
            inline bool contains_all(const point& p0, const point& p1) const;

            //- Contains either points (eg, line ends)
            inline bool contains_any(const point& p0, const point& p1) const;

            //- Contains both line ends?
            inline bool contains_all(const linePointRef& ln) const;

            //- Contains either of the line ends?
            inline bool contains_any(const linePointRef& ln) const;

            //- Same as contains_all(const linePointRef&)
            bool contains(const linePointRef& ln) const
            {
                return contains_all(ln);
            }

            //- Contains all three points? (eg, a triangle)
            inline bool contains_all
            (
                const point& p0,
                const point& p1,
                const point& p2
            ) const;

            //- Contains any of the three points (eg, from a triangle)
            inline bool contains_any
            (
                const point& p0,
                const point& p1,
                const point& p2
            ) const;

            //- Contains all points? (inside or on edge).
            //  \return True on empty range, or if all points are contained.
            bool contains_all(const UList<point>& points) const;

            //- Contains any of the points? (inside or on edge)
            bool contains_any(const UList<point>& points) const;

            //- Same as contains_all(const UList<point>&)
            bool contains(const UList<point>& points) const
            {
                return contains_all(points);
            }

            //- Contains all of the (subsetted) points? (inside or on edge)
            //  \return True on empty ranges, or if all points are contained.
            //
            //  \tparam IndexContainer  A container with an iterator that
            //      dereferences to a label
            //      (eg, FixedList, labelUList, face, triFace, ...)
            template<class IndexContainer>
            bool contains_all
            (
                const UList<point>& points,
                const IndexContainer& indices
            ) const;

            //- Contains any of the (subsetted) points? (inside or on edge)
            //
            //  \tparam IndexContainer  A container with an iterator that
            //      dereferences to a label
            //      (eg, FixedList, labelUList, face, triFace, ...)
            template<class IndexContainer>
            bool contains_any
            (
                const UList<point>& points,
                const IndexContainer& indices
            ) const;

            //- Same as contains_all()
            template<class IndexContainer>
            bool contains
            (
                const UList<point>& points,
                const IndexContainer& indices
            ) const
            {
                return contains_all(points, indices);
            }

            //- Return the nearest point on the boundBox to the supplied point.
            //  If point is inside the boundBox then the point is returned
            //  unchanged.
            point nearest(const point& p) const;

            //- The min/max distance squared from given point to box corners.
            //  The max will correspond to the radius squared of the bounding
            //  sphere (around the given point) that encloses the entire box.
            inline Pair<scalar> distSqr(const point& p) const;

            //- Minimum (absolute) distance from given point to a box corner
            inline scalar minDist(const point& p) const;

            //- Minimum distance squared from given point to a box corner
            inline scalar minDistSqr(const point& p) const;

            //- The radius of the bounding sphere
            //- (around the given point) that encloses the entire box.
            inline scalar maxDist(const point& p) const;

            //- The radius squared of the bounding sphere
            //- (around the given point) that encloses the entire box.
            inline scalar maxDistSqr(const point& p) const;

            //- Compare distance of two bounding boxes to the reference point
            //  return:
            //  -2 : all vertices of box \p a are closer than any of box \p b
            //  -1 : some vertex of box \p a is closer than any of box \p b
            //   0 : boxes appear to have an identical distance to the point
            //  +1 : some vertex of box \p b is closer than any of box \p a
            //  +2 : all vertices of box \p b are closer than any of \p a
            static int distanceCmp
            (
                const boundBox& a,
                const boundBox& b,
                const point& refPoint
            );


    // Member Operators

        //- Restrict min/max to union with other box.
        void operator&=(const boundBox& bb);

        //- Extend box to include the second box, as per the add() method.
        //  Can be used in a reduction operation.
        void operator+=(const boundBox& bb) { add(bb); }


    // Write

        //- Write boundBox in OBJ format (octant order)
        void writeOBJ
        (
            Ostream& os,
            //! [in,out] current vertex offset
            label& vertIndex,
            //! write as lines instead of faces
            const bool writeLinesOnly = false
        ) const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, boundBox& bb);
        friend Ostream& operator<<(Ostream& os, const boundBox& bb);


    // Housekeeping

        //- Same as contains_any()
        bool containsAny(const UList<point>& points) const
        {
            return contains_any(points);
        }

        //- Same as contains_any()
        template<class IdxContainer>
        bool containsAny(const UList<point>& pts, const IdxContainer& ids) const
        {
            return contains_any(pts, ids);
        }

        //- Deprecated(2022-10) - use 'operator&=' to avoid confusion with
        //- other intersects() methods.
        //  \deprecated(2022-10) - use 'operator&='
        FOAM_DEPRECATED_FOR(2022-10, "boundBox::operator&=(const boundBox&)")
        bool intersect(const boundBox& bb) { *this &= bb; return good(); }

        //- Identical to centre()
        point midpoint() const { return centre(); }
};

// Cleanup
#undef Foam_boundBox_syncPar

// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Handle boundBox content as two points
template<> struct contiguous_element<boundBox> { using type = point; };

//- Contiguous data for boundBox
template<> struct is_contiguous<boundBox> : is_contiguous<point> {};

//- Contiguous scalar data for boundBox
template<>
struct is_contiguous_scalar<boundBox> : is_contiguous_scalar<point> {};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const boundBox& a, const boundBox& b)
{
    return (a.min() == b.min() && a.max() == b.max());
}


inline bool operator!=(const boundBox& a, const boundBox& b)
{
    return !(a == b);
}


inline bool operator<(const boundBox& a, const boundBox& b)
{
    return
    (
        a.min() < b.min()
     || (!(b.min() < a.min()) && a.max() < b.max())
    );
}


inline bool operator<=(const boundBox& a, const boundBox& b)
{
    return !(b < a);
}


inline bool operator>(const boundBox& a, const boundBox& b)
{
    return (b < a);
}


inline bool operator>=(const boundBox& a, const boundBox& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "boundBoxI.H"

#ifdef NoRepository
    #include "boundBox.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
