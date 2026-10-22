/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::MatrixSpace

Description
    Templated matrix space.

    Template arguments are the Form the matrix space will be used to create,
    the type of the elements and the number of rows and columns of the matrix.

SourceFiles
    MatrixSpaceI.H

See also
    Foam::VectorSpace

\*---------------------------------------------------------------------------*/

#ifndef Foam_MatrixSpace_H
#define Foam_MatrixSpace_H

#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MatrixSpace Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Cmpt, direction Mrows, direction Ncols>
class MatrixSpace
:
    public VectorSpace<Form, Cmpt, Mrows*Ncols>
{
public:

    // Typedefs

        //- MatrixSpace type
        typedef MatrixSpace<Form, Cmpt, Mrows, Ncols> msType;


    // Member Constants

        //- The number of rows
        static constexpr direction mRows = Mrows;

        //- The number of columns
        static constexpr direction nCols = Ncols;


    // Static Member Functions

        //- The number of rows
        static constexpr direction m() noexcept { return Mrows; }

        //- The number of columns
        static constexpr direction n() noexcept { return Ncols; }

        //- An identity matrix for square matrix-spaces
        inline static msType identity();


    // Sub-Block Classes

        //- Const sub-block type
        template<class SubTensor, direction BRowStart, direction BColStart>
        class ConstBlock
        {
            //- Reference to parent matrix
            const msType& matrix_;

        public:

            static const direction mRows = SubTensor::mRows;
            static const direction nCols = SubTensor::nCols;

            //- Return the number of rows in the block
            static direction m() noexcept { return mRows; }

            //- Return the number of columns in the block
            static direction n() noexcept { return nCols; }

            //- Construct for the given matrix
            inline ConstBlock(const msType& matrix);

            //- Construct and return the sub-tensor corresponding to this block
            inline SubTensor operator()() const;

            //- (i, j) const element access operator
            inline const Cmpt& operator()
            (
                const direction i,
                const direction j
            ) const;
        };


        //- Sub-block type
        template
        <
            class SubTensor,
            direction BRowStart,
            direction BColStart
        >
        class Block
        {
            //- Reference to parent matrix
            msType& matrix_;

        public:

            static const direction mRows = SubTensor::mRows;
            static const direction nCols = SubTensor::nCols;

            //- The number of rows in the block
            static direction m() noexcept { return mRows; }

            //- The number of columns in the block
            static direction n() noexcept { return nCols; }

            //- Construct for the given matrix
            inline Block(msType& matrix);

            //- Assignment to a matrix
            template<class Form2>
            inline void operator=
            (
                const MatrixSpace
                <
                    Form2,
                    Cmpt,
                    SubTensor::mRows,
                    SubTensor::nCols
                >& matrix
            );

            //- Assignment to a column vector
            template<class VSForm>
            inline void operator=
            (
                const VectorSpace<VSForm, Cmpt, SubTensor::mRows>& v
            );

            //- Construct and return the sub-tensor corresponding to this block
            inline SubTensor operator()() const;

            //- (i, j) const element access operator
            inline const Cmpt& operator()
            (
                const direction i,
                const direction j
            ) const;

            //- (i, j) element access operator
            inline Cmpt& operator()(const direction i, const direction j);
        };


    // Generated Methods

        //- Default construct
        MatrixSpace() = default;

        //- Copy construct
        constexpr MatrixSpace(const MatrixSpace&) = default;

        //- Move construct
        constexpr MatrixSpace(MatrixSpace&&) = default;

        //- Copy assignment
        constexpr MatrixSpace& operator=(const MatrixSpace&) = default;

        //- Move assignment
        constexpr MatrixSpace& operator=(MatrixSpace&&) = default;


    // Constructors

        //- Construct initialized to zero
        inline constexpr MatrixSpace(Foam::zero);

        //- Construct as copy of a VectorSpace with the same size
        template<class Form2, class Cmpt2>
        inline constexpr explicit MatrixSpace
        (
            const VectorSpace<Form2, Cmpt2, Mrows*Ncols>&
        );

        //- Aggregate initialization, forwarded to VectorSpace
        template<class... Args>
        inline constexpr explicit MatrixSpace
        (
            Foam::direct_construct_t,
            Args&&... vals
        );

        //- Construct from a block of another matrix space
        template
        <
            template<class, direction, direction> class Block2,
            direction BRowStart,
            direction BColStart
        >
        inline MatrixSpace
        (
            const Block2<Form, BRowStart, BColStart>& block
        );

        //- Construct from Istream
        explicit MatrixSpace(Istream& is);


    // Member Functions

        //- Access the matrix-space element: compile-time check of indices
        template<direction Row, direction Col>
        inline constexpr const Cmpt& elmt() const noexcept;

        //- Access the matrix-space element: compile-time check of indices
        template<direction Row, direction Col>
        inline constexpr Cmpt& elmt() noexcept;


        // Const element access functions for a 3x3
        // Compile-time errors are generated for inappropriate use

            constexpr const Cmpt& xx() const noexcept { return elmt<0,0>(); }
            constexpr const Cmpt& xy() const noexcept { return elmt<0,1>(); }
            constexpr const Cmpt& xz() const noexcept { return elmt<0,2>(); }
            constexpr const Cmpt& yx() const noexcept { return elmt<1,0>(); }
            constexpr const Cmpt& yy() const noexcept { return elmt<1,1>(); }
            constexpr const Cmpt& yz() const noexcept { return elmt<1,2>(); }
            constexpr const Cmpt& zx() const noexcept { return elmt<2,0>(); }
            constexpr const Cmpt& zy() const noexcept { return elmt<2,1>(); }
            constexpr const Cmpt& zz() const noexcept { return elmt<2,2>(); }


        // Element access functions for a 3x3
        // Compile-time errors are generated for inappropriate use

            constexpr Cmpt& xx() noexcept { return elmt<0,0>(); }
            constexpr Cmpt& xy() noexcept { return elmt<0,1>(); }
            constexpr Cmpt& xz() noexcept { return elmt<0,2>(); }
            constexpr Cmpt& yx() noexcept { return elmt<1,0>(); }
            constexpr Cmpt& yy() noexcept { return elmt<1,1>(); }
            constexpr Cmpt& yz() noexcept { return elmt<1,2>(); }
            constexpr Cmpt& zx() noexcept { return elmt<2,0>(); }
            constexpr Cmpt& zy() noexcept { return elmt<2,1>(); }
            constexpr Cmpt& zz() noexcept { return elmt<2,2>(); }


        //- Return the transpose of the matrix
        inline typename typeOfTranspose<Cmpt, Form>::type T() const;

        //- Return a const sub-block corresponding to the specified type
        //  starting at the specified row and column
        template<class SubTensor, direction BRowStart, direction BColStart>
        inline ConstBlock<SubTensor, BRowStart, BColStart> block() const;

        //- Return a sub-block corresponding to the specified type
        //  starting at the specified row and column
        template<class SubTensor, direction BRowStart, direction BColStart>
        inline Block<SubTensor, BRowStart, BColStart> block();

        //- (i, j) const element access operator
        inline const Cmpt& operator()
        (
            const direction i,
            const direction j
        ) const;

        //- (i, j) element access operator
        inline Cmpt& operator()(const direction i, const direction j);


    // Member Operators

        //- Assignment to zero
        inline void operator=(Foam::zero);

        //- Assignment to a block of another matrix space
        template
        <
            template<class, direction, direction> class Block2,
            direction BRowStart,
            direction BColStart
        >
        inline void operator=
        (
            const Block2<Form, BRowStart, BColStart>& block
        );

        //- Inner product with a compatible square matrix
        template<class Form2>
        inline void operator&=
        (
            const MatrixSpace<Form, Cmpt, Ncols, Ncols>& matrix
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MatrixSpaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
