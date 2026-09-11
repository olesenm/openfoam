/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::Vector2D

Description
    Templated 2D Vector derived from VectorSpace adding construction from
    2 components, element access using x() and y() member functions and
    the inner-product (dot-product).

SourceFiles
    Vector2DI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Vector2D_H
#define Foam_Vector2D_H

#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Cmpt> class Vector2D;

// Common Types

//! \class Foam::floatVector2D
//! \brief A Vector2D of values with float precision
typedef Vector2D<float> floatVector2D;

//! \class Foam::doubleVector2D
//! \brief A Vector2D of values with double precision
typedef Vector2D<double> doubleVector2D;

/*---------------------------------------------------------------------------*\
                           Class Vector2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Vector2D
:
    public VectorSpace<Vector2D<Cmpt>, Cmpt, 2>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Vector2D<label> labelType;


    // Member Constants

        //- Rank of Vector2D is 1
        static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { X, Y };


    // Generated Methods

        //- Default construct
        Vector2D() = default;

        //- Copy construct
        constexpr Vector2D(const Vector2D&) = default;

        //- Move construct
        constexpr Vector2D(Vector2D&&) = default;

        //- Copy assignment
        constexpr Vector2D& operator=(const Vector2D&) = default;

        //- Move assignment
        constexpr Vector2D& operator=(Vector2D&&) = default;


    // Constructors

        //- Construct initialized to zero
        inline constexpr Vector2D(Foam::zero);

        //- Copy construct from VectorSpace of the same rank
        inline Vector2D(const VectorSpace<Vector2D<Cmpt>, Cmpt, 2>& vs);

        //- Construct from two components
        inline constexpr Vector2D(const Cmpt& vx, const Cmpt& vy);

        //- Copy construct from 3D variant, discarding z component
        inline constexpr explicit Vector2D(const Vector<Cmpt>& vs);

        //- Construct from Istream
        inline explicit Vector2D(Istream& is);


    // Member Functions

    // Component Access

        //- Access to the vector x component
        constexpr const Cmpt& x() const noexcept { return this->v_[X]; }

        //- Access to the vector y component
        constexpr const Cmpt& y() const noexcept { return this->v_[Y]; }

        //- Access to the vector x component
        constexpr Cmpt& x() noexcept { return this->v_[X]; }

        //- Access to the vector y component
        constexpr Cmpt& y() noexcept { return this->v_[Y]; }

        // //- A 'fake' value for z component
        // static Cmpt z() { return Cmpt(Foam::zero{}); }


    // Vector Operations

        //- The length (L2-norm) of the vector
        inline scalar mag() const;

        //- The length (L2-norm) squared of the vector.
        inline constexpr scalar magSqr() const;

        //- The L2-norm distance from another vector.
        //- The mag() of the difference.
        inline scalar dist(const Vector2D<Cmpt>& v2) const;

        //- The L2-norm distance squared from another vector.
        //- The magSqr() of the difference.
        inline constexpr scalar distSqr(const Vector2D<Cmpt>& v2) const;

        //- Normalise the vector by its magnitude
        //  For small magnitudes (less than ROOTVSMALL) set to zero.
        //  Will not be particularly useful for a vector of labels
        inline Vector2D<Cmpt>& normalise(const scalar tol = ROOTVSMALL);

        //- Inplace removal of components that are collinear to the given
        //- unit vector.
        inline Vector2D<Cmpt>& removeCollinear(const Vector2D<Cmpt>& unitVec);

        //- Scalar-product of \c this with another Vector.
        inline constexpr Cmpt inner(const Vector2D<Cmpt>& v2) const;

        //- Perp dot product (dot product with perpendicular vector)
        inline constexpr Cmpt perp(const Vector2D<Cmpt>& b) const;

        //- Return true if vector is within tol
        inline bool isClose
        (
            const Vector2D<Cmpt>& b,
            const scalar tol = 1e-10
        ) const;


    // Comparison Operations

        //- Lexicographically compare \em a and \em b with order (x:y)
        static inline bool less_xy
        (
            const Vector2D<Cmpt>& a,
            const Vector2D<Cmpt>& b
        );

        //- Lexicographically compare \em a and \em b with order (y:x)
        static inline bool less_yx
        (
            const Vector2D<Cmpt>& a,
            const Vector2D<Cmpt>& b
        );


    // Member Operators

        //- Inherit VectorSpace += operations
        using Vector2D::vsType::operator+=;

        //- Inherit VectorSpace -= operations
        using Vector2D::vsType::operator-=;

        //- Add compatible 2D vector to this
        template<class Cmpt2>
        std::enable_if_t<std::is_convertible_v<Cmpt2, Cmpt>, void>
        inline operator+=(const Vector2D<Cmpt2>& b)
        {
            this->x() += b.x();
            this->y() += b.y();
        }

        //- Subtract compatible 2D vector from this
        template<class Cmpt2>
        std::enable_if_t<std::is_convertible_v<Cmpt2, Cmpt>, void>
        inline operator-=(const Vector2D<Cmpt2>& b)
        {
            this->x() -= b.x();
            this->y() -= b.y();
        }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Vector2D<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Vector2D<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Vector2D<Cmpt>> : is_contiguous_scalar<Cmpt> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Vector2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
