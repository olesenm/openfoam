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

\*---------------------------------------------------------------------------*/

#include "boundBox.H"
#include "treeBoundBox.H"
#include "plane.H"
#include "hexCell.H"
#include "triangle.H"
#include "MinMax.H"
#include "Random.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Broadcast contents (single or multiple boxes)
template<class BoundBoxType, class IntType>
void broadcastBoxes_impl
(
    BoundBoxType bbs[],
    IntType count,
    int communicator
)
{
    using namespace Foam;

    if (communicator < 0) communicator = UPstream::worldComm;

    if (count <= 0 || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }

    // Note: UPstream_dataType<boundBox> will be "type_6double" etc
    // (for the min/max point pair), so no need to adjust the count here.
    UPstream::broadcast(bbs, count, communicator);
}


// Bundle multiple boxes together and reduce in two MPI calls
template<class BoundBoxType, class IntType>
void reduceBoxes_impl
(
    BoundBoxType bbs[],
    IntType count,
    int communicator
)
{
    using namespace Foam;

    if (communicator < 0) communicator = UPstream::worldComm;

    if (count <= 0 || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if (count == 1)
    {
        // Single element
        bbs[0].reduce(communicator);
        return;
    }

    // Bundle min/max together and perform in two reductions.
    List<point> work(count);

    // Min values: pack, reduce, unpack
    {
        constexpr auto bop = UPstream::opCodes::op_min;

        forAll(work, i) { work[i] = bbs[i].min(); }

        UPstream::mpiAllReduce<bop>(work.data(), work.size(), communicator);

        forAll(work, i) { bbs[i].min() = work[i]; }
    }

    // Max values: pack, reduce, unpack
    {
        constexpr auto bop = UPstream::opCodes::op_max;

        forAll(work, i) { work[i] = bbs[i].max(); }

        UPstream::mpiAllReduce<bop>(work.data(), work.size(), communicator);

        forAll(work, i) { bbs[i].max() = work[i]; }
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::boundBox Foam::boundBox::greatBox
(
    point::uniform(-ROOTVGREAT),
    point::uniform(ROOTVGREAT)
);

const Foam::boundBox Foam::boundBox::invertedBox
(
    point::uniform(ROOTVGREAT),
    point::uniform(-ROOTVGREAT)
);

const Foam::FixedList<Foam::vector, 6> Foam::boundBox::faceNormals
({
    vector(-1,  0,  0), // 0: x-min, left
    vector( 1,  0,  0), // 1: x-max, right
    vector( 0, -1,  0), // 2: y-min, bottom
    vector( 0,  1,  0), // 3: y-max, top
    vector( 0,  0, -1), // 4: z-min, back
    vector( 0,  0,  1)  // 5: z-max, front
});


// Faces, in octant point order
static int octModelFaces_[6][4]
{
    {0, 4, 6, 2},  // 0: x-min, left
    {1, 3, 7, 5},  // 1: x-max, right
    {0, 1, 5, 4},  // 2: y-min, bottom
    {2, 6, 7, 3},  // 3: y-max, top
    {0, 2, 3, 1},  // 4: z-min, back
    {4, 5, 7, 6}   // 5: z-max, front
};

// Edges, in octant point order
static int octModelEdges_[12][2]
{
    {0, 1}, // 0
    {1, 3},
    {2, 3}, // 2
    {0, 2},
    {4, 5}, // 4
    {5, 7},
    {6, 7}, // 6
    {4, 6},
    {0, 4}, // 8
    {1, 5},
    {3, 7}, // 10
    {2, 6}
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::faceList& Foam::boundBox::hexFaces()
{
    return hexCell::modelFaces();
}


const Foam::faceList& Foam::boundBox::octFaces()
{
    static std::unique_ptr<Foam::faceList> ptr(nullptr);

    if (!ptr)
    {
        ptr.reset(new Foam::faceList(boundBox::nFaces()));

        int facei = 0;
        for (auto& f : *ptr)
        {
            f.resize(4);
            f[0] = octModelFaces_[facei][0];
            f[1] = octModelFaces_[facei][1];
            f[2] = octModelFaces_[facei][2];
            f[3] = octModelFaces_[facei][3];
            ++facei;
        }
    }

    return *ptr;
}


const Foam::edgeList& Foam::boundBox::octEdges()
{
    static std::unique_ptr<Foam::edgeList> ptr(nullptr);

    if (!ptr)
    {
        ptr.reset(new Foam::edgeList(boundBox::nEdges()));

        int edgei = 0;
        for (auto& e : *ptr)
        {
            e[0] = octModelEdges_[edgei][0];
            e[1] = octModelEdges_[edgei][1];
            ++edgei;
        }
    }

    return *ptr;
}



int Foam::boundBox::distanceCmp
(
    const boundBox& a,
    const boundBox& b,
    const point& refPoint
)
{
    // Distances (squared) from point to nearest/furthest away vertex
    Pair<scalar> distA = a.distSqr(refPoint);
    Pair<scalar> distB = b.distSqr(refPoint);

    // Categorize
    if (distA.second() < distB.first())
    {
        // All vertices of 'A' are nearer to point than any vertex of 'B'
        return -2;
    }
    else if (distB.second() < distA.first())
    {
        // All vertices of 'B' are nearer to point than any vertex of 'A'
        return 2;
    }
    else
    {
        // Mixed: lexicographic compare (min, then max) as tie-breaker
        if (distA < distB) return -1;
        if (distB < distA) return 1;
        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBox::boundBox(const boundBox& bb, bool syncPar)
:
    boundBox(bb)
{
    if (syncPar)
    {
        this->reduce(UPstream::worldComm);
    }
}


Foam::boundBox::boundBox(const UList<point>& points, bool syncPar)
:
    boundBox()
{
    add(points);

    if (syncPar)
    {
        this->reduce(UPstream::worldComm);
    }
}


Foam::boundBox::boundBox(const tmp<pointField>& tpoints, bool syncPar)
:
    boundBox()
{
    add(tpoints);

    if (syncPar)
    {
        this->reduce(UPstream::worldComm);
    }
}


Foam::boundBox::boundBox
(
    const UList<point>& points,
    const labelUList& indices,
    bool syncPar
)
:
    boundBox()
{
    add(points, indices);

    if (syncPar)
    {
        this->reduce(UPstream::worldComm);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::boundBox::hexCorners(UList<point>& output) const
{
    if (output.size() < boundBox::nPoints())
    {
        FatalErrorInFunction
            << "Cannot copy 8 corner points into a list of size "
            << output.size() << nl
            << abort(FatalError);
    }

    output[0] = hexCorner<0>();
    output[1] = hexCorner<1>();
    output[2] = hexCorner<2>();
    output[3] = hexCorner<3>();
    output[4] = hexCorner<4>();
    output[5] = hexCorner<5>();
    output[6] = hexCorner<6>();
    output[7] = hexCorner<7>();
}


void Foam::boundBox::octCorners(UList<point>& output) const
{
    if (output.size() < boundBox::nPoints())
    {
        FatalErrorInFunction
            << "Cannot copy 8 corner points into a list of size "
            << output.size() << nl
            << abort(FatalError);
    }

    output[0] = octCorner<0>();
    output[1] = octCorner<1>();
    output[2] = octCorner<2>();
    output[3] = octCorner<3>();
    output[4] = octCorner<4>();
    output[5] = octCorner<5>();
    output[6] = octCorner<6>();
    output[7] = octCorner<7>();
}


Foam::pointField Foam::boundBox::hexCorners() const
{
    pointField pts(boundBox::nPoints());
    hexCorners(pts);
    return pts;
}


Foam::pointField Foam::boundBox::octCorners() const
{
    pointField pts(boundBox::nPoints());
    octCorners(pts);
    return pts;
}


void Foam::boundBox::faceCentres(UList<point>& output) const
{
    if (output.size() < boundBox::nFaces())
    {
        FatalErrorInFunction
            << "Cannot copy 6 faceCentre points into a list of size "
            << output.size() << nl
            << abort(FatalError);
    }

    for (direction facei = 0; facei < boundBox::nFaces(); ++facei)
    {
        output[facei] = faceCentre(facei);
    }
}


Foam::pointField Foam::boundBox::faceCentres() const
{
    pointField pts(boundBox::nFaces());
    faceCentres(pts);
    return pts;
}


Foam::point Foam::boundBox::faceCentre(direction facei) const
{
    point pt = boundBox::centre();

    switch (facei)
    {
        case 0: pt.x() = min().x(); break;  // 0: x-min, left
        case 1: pt.x() = max().x(); break;  // 1: x-max, right
        case 2: pt.y() = min().y(); break;  // 2: y-min, bottom
        case 3: pt.y() = max().y(); break;  // 3: y-max, top
        case 4: pt.z() = min().z(); break;  // 4: z-min, back
        case 5: pt.z() = max().z(); break;  // 5: z-max, front
        default:
        {
            FatalErrorInFunction
                << "Face:" << int(facei) << " should be [0..5]"
                << abort(FatalError);
        }
    }

    return pt;
}


void Foam::boundBox::reduce(int communicator)
{
    if (communicator < 0) communicator = UPstream::worldComm;

    const auto tag = UPstream::msgType();  // irrelevant (is_contiguous)
    Foam::reduce(min_, minOp<point>(), tag, communicator);
    Foam::reduce(max_, maxOp<point>(), tag, communicator);
}


void Foam::boundBox::reduce
(
    boundBox bbs[],
    label count,
    int communicator
)
{
    reduceBoxes_impl(bbs, count, communicator);
}


void Foam::boundBox::reduce
(
    treeBoundBox bbs[],
    label count,
    int communicator
)
{
    reduceBoxes_impl(bbs, count, communicator);
}


void Foam::boundBox::broadcast(int communicator)
{
    broadcastBoxes_impl(this, 1, communicator);
}


void Foam::boundBox::broadcast
(
    boundBox bbs[],
    label count,
    int communicator
)
{
    broadcastBoxes_impl(bbs, count, communicator);
}


void Foam::boundBox::broadcast
(
    treeBoundBox bbs[],
    label count,
    int communicator
)
{
    broadcastBoxes_impl(bbs, count, communicator);
}


bool Foam::boundBox::intersects(const plane& pln) const
{
    // Require a full 3D box
    if (nDim() != 3)
    {
        return false;
    }

    // Check as below(1) or above(2) - stop when it cuts both
    int side = 0;

    #undef  doLocalCode
    #define doLocalCode(Idx)                                                  \
    {                                                                         \
        side |= (pln.whichSide(hexCorner<Idx>()) == plane::BACK ? 1 : 2);     \
        if (side == 3) return true;  /* Both below and above: done */         \
    }

    // Each box corner
    doLocalCode(0);
    doLocalCode(1);
    doLocalCode(2);
    doLocalCode(3);
    doLocalCode(4);
    doLocalCode(5);
    doLocalCode(6);
    doLocalCode(7);

    #undef doLocalCode

    return false;
}


bool Foam::boundBox::intersects(const triPointRef& tri) const
{
    // Require a full 3D box
    if (nDim() != 3)
    {
        return false;
    }

    // Simplest check - if any points are inside
    if (contains(tri.a()) || contains(tri.b()) || contains(tri.c()))
    {
        return true;
    }


    // Extent of box points projected onto axis
    const auto project_box = []
    (
        const boundBox& bb,
        const vector& axis,
        scalarMinMax& extent
    ) -> void
    {
        extent.reset(axis & bb.hexCorner<0>());
        extent.add(axis & bb.hexCorner<1>());
        extent.add(axis & bb.hexCorner<2>());
        extent.add(axis & bb.hexCorner<3>());
        extent.add(axis & bb.hexCorner<4>());
        extent.add(axis & bb.hexCorner<5>());
        extent.add(axis & bb.hexCorner<6>());
        extent.add(axis & bb.hexCorner<7>());
    };


    // Use separating axis theorem to determine if triangle and
    // (axis-aligned) bounding box intersect.

    scalarMinMax tri_extent(0);
    scalarMinMax box_extent(0);
    const boundBox& bb = *this;

    // 1.
    // Test separating axis defined by the box normals
    // (project triangle points)
    // - do first (largely corresponds to normal bound box rejection test)
    //
    // No intersection if extent of projected triangle points are outside
    // of the box range

    {
        // vector::X
        tri_extent.reset(tri.a().x());
        tri_extent.add(tri.b().x());
        tri_extent.add(tri.c().x());

        box_extent.reset(bb.min().x(), bb.max().x());

        if (!tri_extent.overlaps(box_extent))
        {
            return false;
        }

        // vector::Y
        tri_extent.reset(tri.a().y());
        tri_extent.add(tri.b().y());
        tri_extent.add(tri.c().y());

        box_extent.reset(bb.min().y(), bb.max().y());

        if (!tri_extent.overlaps(box_extent))
        {
            return false;
        }

        // vector::Z
        tri_extent.reset(tri.a().z());
        tri_extent.add(tri.b().z());
        tri_extent.add(tri.c().z());

        box_extent.reset(bb.min().z(), bb.max().z());

        if (!tri_extent.overlaps(box_extent))
        {
            return false;
        }
    }


    // 2.
    // Test separating axis defined by the triangle normal
    // (project box points)
    // - can use area or unit normal since any scaling is applied to both
    //   sides of the comparison.
    // - by definition all triangle points lie in the plane defined by
    //   the normal. It doesn't matter which of the points we use to define
    //   the triangle offset (extent) when projected onto the triangle normal

    vector axis = tri.areaNormal();

    tri_extent.reset(axis & tri.a());
    project_box(bb, axis, box_extent);

    if (!tri_extent.overlaps(box_extent))
    {
        return false;
    }


    // 3.
    // Test separating axes defined by the triangle edges, which are the
    // cross product of the edge vectors and the box face normals

    for (const vector& edgeVec : { tri.vecA(), tri.vecB(), tri.vecC() })
    {
        for (direction faceDir = 0; faceDir < vector::nComponents; ++faceDir)
        {
            axis = Zero;
            axis[faceDir] = 1;

            axis = (edgeVec ^ axis);

            // project tri
            tri_extent.reset(axis & tri.a());
            tri_extent.add(axis & tri.b());
            tri_extent.add(axis & tri.c());

            project_box(bb, axis, box_extent);

            if (!tri_extent.overlaps(box_extent))
            {
                return false;
            }
        }
    }

    return true;
}


bool Foam::boundBox::intersects(const point& from, const point& to) const
{
    return treeBoundBox::intersects(*this, from, to);
}


bool Foam::boundBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (!contains(p))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (contains(p))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::boundBox::nearest(const point& p) const
{
    // Clip the point to the range of the bounding box
    return point
    (
        std::min(std::max(p.x(), min_.x()), max_.x()),
        std::min(std::max(p.y(), min_.y()), max_.y()),
        std::min(std::max(p.z(), min_.z()), max_.z())
    );
}


void Foam::boundBox::inflate(Random& rndGen, scalar factor)
{
    vector newSpan(span());

    // Make 3D
    const scalar minSpan = factor*mag();

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        newSpan[dir] = std::max(newSpan[dir], minSpan);
    }

    min_ -= cmptMultiply(factor*rndGen.sample01<vector>(), newSpan);
    max_ += cmptMultiply(factor*rndGen.sample01<vector>(), newSpan);

    // May still need extra handling for zero-sized boxes.
    // Eg, use inflate with the extra delta (ie, grow)
}


void Foam::boundBox::inflate
(
    Random& rndGen,
    scalar factor,
    scalar delta
)
{
    inflate(rndGen, factor);
    grow(delta);
}


void Foam::boundBox::inflate_random(label seedValue, scalar factor)
{
    Random rndGen(seedValue);
    inflate(rndGen, factor);
}


void Foam::boundBox::inflate_random
(
    label seedValue,
    scalar factor,
    scalar delta
)
{
    Random rndGen(seedValue);
    inflate(rndGen, factor, delta);
}


void Foam::boundBox::writeOBJ
(
    Ostream& os,
    label& vertIndex,
    const bool writeLinesOnly
) const
{
    // Annotate with '#box' which can be grep'd for later
    os  << "#box: " << min_ << ' ' << max_ << nl;

    for (direction oct = 0; oct < 8; ++oct)
    {
        point p = corner(oct);
        os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }

    if (writeLinesOnly)
    {
        for (const auto& e : boundBox::octEdges())
        {
            os  << "l " << e.first() + vertIndex + 1
                << ' ' << e.second() + vertIndex + 1 << nl;
        }
    }
    else
    {
        for (const auto& f : boundBox::octFaces())
        {
            os  << 'f';
            for (label fpi : f)
            {
                os  << ' ' << fpi + vertIndex + 1;
            }
            os  << nl;
        }
    }

    vertIndex += 8;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::boundBox::operator&=(const boundBox& bb)
{
    min_ = ::Foam::max(min_, bb.min_);
    max_ = ::Foam::min(max_, bb.max_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBox& bb)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundBox)
        );
    }
    else
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBox& bb)
{
    if (is.format() == IOstreamOption::BINARY)
    {
        Detail::readContiguous<boundBox>
        (
            is,
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBox)
        );
    }
    else
    {
        is >> bb.min_ >> bb.max_;
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
