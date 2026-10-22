/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    Foam::VectorSpace

Description
    Templated vector space.

    Template arguments are the Form the vector space will be used to create,
    the type of the elements and the number of elements.

    The VectorSpace will normally contain trivial component types,
    but since there are certain places that use/misuse Vector to store
    a tuple of other things, we cannot categorically assert that the
    components are "trivially copyable".

SourceFiles
    VectorSpaceI.H
    VectorSpace.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_VectorSpace_H
#define Foam_VectorSpace_H

#include "contiguous.H"
#include "direction.H"
#include "scalar.H"
#include "word.H"
#include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class VectorSpace Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Cmpt, direction Ncmpts>
class VectorSpace
{
public:

    //- The components of this vector space
    Cmpt v_[Ncmpts];

    //- VectorSpace type
    typedef VectorSpace<Form, Cmpt, Ncmpts> vsType;

    //- Component type
    typedef Cmpt cmptType;

    //- Magnitude type
    typedef Cmpt magType;

    //- The type to represent the size of a VectorSpace
    typedef direction size_type;

    // #ifdef FULLDEBUG
    // static_assert(Ncmpts, "Number of components must be non-zero, positive");
    // #endif


    // Static Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Number of components in this vector space
        static constexpr direction nComponents = Ncmpts;


        // VectorSpace currently defaults to a column-vector
        // This will be removed when column-vector is introduced
        // as a specialization
        static constexpr direction mRows = Ncmpts;
        static constexpr direction nCols = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const Form zero;
        static const Form one;
        static const Form max;
        static const Form min;
        static const Form rootMax;
        static const Form rootMin;


    // Sub-Block Classes

        //- Const sub-block type
        template<class SubVector, direction BStart>
        class ConstBlock
        {
            const vsType& vs_;

        public:

            //- Number of components in this vector space
            static const direction nComponents = SubVector::nComponents;

            //- Construct for a given vector
            inline ConstBlock(const vsType& vs);

            //- [i] const element access operator
            inline const Cmpt& operator[](const direction i) const;

            //- (i, 0) const element access operator
            inline const Cmpt& operator()
            (
                const direction i,
                const direction
            ) const;
        };


    // Generated Methods

        //- Default construct
        VectorSpace() = default;

        //- Copy construct
        constexpr VectorSpace(const VectorSpace&) = default;

        //- Move construct
        constexpr VectorSpace(VectorSpace&&) = default;

        //- Copy assignment
        constexpr VectorSpace& operator=(const VectorSpace&) = default;

        //- Move assignment
        constexpr VectorSpace& operator=(VectorSpace&&) = default;


    // Constructors

        //- Construct initialized to zero
        inline constexpr VectorSpace(Foam::zero);

        //- Copy construct of a VectorSpace with the same size
        template<class Form2, class Cmpt2>
        inline constexpr explicit VectorSpace
        (
            const VectorSpace<Form2, Cmpt2, Ncmpts>&
        );

        //- Aggregate initialization using Ncmpts arguments
        template<class... Args>
        inline constexpr explicit VectorSpace
        (
            Foam::direct_construct_t,
            Args&&... vals
        );

        //- Construct from Istream
        explicit VectorSpace(Istream& is);


    // Member Functions

        //- The number of elements in the VectorSpace = Ncmpts.
        static constexpr direction size() noexcept { return Ncmpts; }

        //- The maximum size, corresponds to the VectorSpace = Ncmpts.
        static constexpr direction max_size() noexcept { return Ncmpts; }

        //- Access the vector-space component
        inline constexpr const Cmpt& component(const direction) const;

        //- Access the vector-space component
        inline constexpr Cmpt& component(const direction);

        //- Copy out the value of the vector-space component
        inline void component(Cmpt&, const direction) const;

        //- Copy in the value of the vector-space component
        inline void replace(const direction, const Cmpt&);

        //- Return const pointer to the first data element
        constexpr const Cmpt* cdata() const noexcept { return v_; }

        //- Return pointer to the first data element
        constexpr Cmpt* data() noexcept { return v_; }

        //- Assign all components to given value
        inline constexpr void fill(const Cmpt& s);

        //- Return a vector-space Form with uniform components.
        //
        //  \note only constexpr if Form itself has a constexpr default
        //        constructor, otherwise it remains a runtime call.
        inline static constexpr Form uniform(const Cmpt& s);

        template<class SubVector, direction BStart>
        inline const ConstBlock<SubVector, BStart> block() const;


    // Member Operators

        //- Access the vector-space component
        inline constexpr const Cmpt& operator[](const direction) const;

        //- Access the vector-space component
        inline constexpr Cmpt& operator[](const direction);

        inline void operator+=(const VectorSpace<Form, Cmpt, Ncmpts>&);
        inline void operator-=(const VectorSpace<Form, Cmpt, Ncmpts>&);

        //- Copy assign from a VectorSpace with the same size
        template<class Form2, class Cmpt2>
        inline constexpr void operator=
        (
            const VectorSpace<Form2, Cmpt2, Ncmpts>& vs
        );

        //- Assign all components to zero
        inline constexpr void operator=(Foam::zero);
        inline void operator*=(const scalar);
        inline void operator/=(const scalar);


    // Iterators

        //- Random access iterator for traversing VectorSpace
        typedef Cmpt* iterator;

        //- Random access iterator for traversing VectorSpace
        typedef const Cmpt* const_iterator;


    // Random access iterators (const and non-const)

        //- Return an iterator (pointer) to begin of VectorSpace
        constexpr iterator begin() noexcept { return v_; }

        //- Return const_iterator (const pointer) to begin of VectorSpace
        constexpr const_iterator begin() const noexcept { return v_; }

        //- Return const_iterator (const pointer) to begin of VectorSpace
        constexpr const_iterator cbegin() const noexcept { return v_; }

        //- Return an iterator (pointer) to end of VectorSpace
        constexpr iterator end() noexcept { return (v_ + Ncmpts); }

        //- Return const_iterator (const pointer) to end of VectorSpace
        constexpr const_iterator end() const noexcept { return (v_ + Ncmpts); }

        //- Return const_iterator (const pointer) to end of VectorSpace
        constexpr const_iterator cend() const noexcept { return (v_ + Ncmpts); }
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- A word representation of a VectorSpace
template<class Form, class Cmpt, direction Ncmpts>
word name(const VectorSpace<Form, Cmpt, Ncmpts>& vs);


// * * * * * * * * * * * * * * * IOstream Operators * * * * * * * * * * * * //

template<class Form, class Cmpt, direction Ncmpts>
Istream& operator>>
(
    Istream& is,
    VectorSpace<Form, Cmpt, Ncmpts>& vs
);

template<class Form, class Cmpt, direction Ncmpts>
Ostream& operator<<
(
    Ostream& os,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VectorSpaceI.H"

#ifdef NoRepository
    #include "VectorSpace.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
