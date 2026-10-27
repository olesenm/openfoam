/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::lduMatrix

Description
    lduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.

    Addressing arrays must be supplied for the upper and lower triangles.

    It might be better if this class were organised as a hierarchy starting
    from an empty matrix, then deriving diagonal, symmetric and asymmetric
    matrices.

SourceFiles
    lduMatrixATmul.C
    lduMatrix.C
    lduMatrixTemplates.C
    lduMatrixOperations.C
    lduMatrixSolver.C
    lduMatrixPreconditioner.C
    lduMatrixTests.C
    lduMatrixUpdateMatrixInterfaces.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_lduMatrix_H
#define Foam_lduMatrix_H

#include "lduMesh.H"
#include "primitiveFieldsFwd.H"
#include "DynamicField.H"
#include "SubField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "solverPerformance.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "InfoProxy.H"
#include "Enum.H"
#include "profilingTrigger.H"
#include <functional>  // For reference_wrapper

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class lduMatrix;

Ostream& operator<<(Ostream&, const lduMatrix&);
Ostream& operator<<(Ostream&, const InfoProxy<lduMatrix>&);


/*---------------------------------------------------------------------------*\
                           Class lduMatrix Declaration
\*---------------------------------------------------------------------------*/

class lduMatrix
{
    // Private Data

        //- LDU mesh reference
        //const lduMesh& lduMesh_;
        std::reference_wrapper<const lduMesh> lduMesh_;

        //- Diagonal coefficients
        std::unique_ptr<DynamicField<scalar>> diagPtr_;

        //- Off-diagonal coefficients (not including interfaces)
        std::unique_ptr<DynamicField<scalar>> lowerPtr_;

        //- Off-diagonal coefficients (not including interfaces)
        std::unique_ptr<DynamicField<scalar>> upperPtr_;

        //- Off-diagonal coefficients (not including interfaces) in CSR ordering
        mutable std::unique_ptr<DynamicField<scalar>> lowerCSRPtr_;

        //- Work space
        mutable std::unique_ptr<solveScalarField> workPtr_;


public:

    // Public Types

        //- Enumerated matrix normalisation types
        enum class normTypes : char
        {
            NO_NORM,            //!< "none" norm (returns 1)
            DEFAULT_NORM,       //!< "default" norm (== L1_scaled)
            L1_SCALED_NORM,     //!< "L1_scaled" norm
        };

        //- Names for the normTypes
        static const Enum<normTypes> normTypesNames_;

        //- Default maximum number of iterations for solvers (1000)
        static constexpr const label defaultMaxIter = 1000;

        //- Default (absolute) tolerance (1e-6)
        static const scalar defaultTolerance;


    // -----------------------------------------------------------------------
    //- Abstract base-class for lduMatrix solvers
    class solver
    {
    protected:

        // Protected Data

            word fieldName_;
            const lduMatrix& matrix_;

            //- Flat storage for boundary interface coefficients (may be owned or ref)
            scalarField interfaceBouCoeffsData_;

            //- Flat storage for internal interface coefficients (may be owned or ref)
            scalarField interfaceIntCoeffsData_;

            //- Per-interface offsets into flat storage (size nInterfaces+1)
            labelList coeffsOffsets_;

            lduInterfaceFieldPtrsList interfaces_;

            //- Dictionary of solution controls
            dictionary controlDict_;

            //- Verbosity level for solver output statements
            int log_;

            //- Minimum number of iterations in the solver
            label minIter_;

            //- Maximum number of iterations in the solver
            label maxIter_;

            //- The normalisation type
            lduMatrix::normTypes normType_;

            //- Final convergence tolerance
            scalar tolerance_;

            //- Convergence tolerance relative to the initial
            scalar relTol_;

            //- Profiling instrumentation
            profilingTrigger profiling_;


        // Protected Member Functions

            //- Read the control parameters from controlDict_
            virtual void readControls();


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                symMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const scalarField& interfaceBouCoeffsData,
                    const scalarField& interfaceIntCoeffsData,
                    const labelUList& coeffsOffsets,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& solverControls
                ),
                (
                    fieldName,
                    matrix,
                    interfaceBouCoeffsData,
                    interfaceIntCoeffsData,
                    coeffsOffsets,
                    interfaces,
                    solverControls
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                asymMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const scalarField& interfaceBouCoeffsData,
                    const scalarField& interfaceIntCoeffsData,
                    const labelUList& coeffsOffsets,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& solverControls
                ),
                (
                    fieldName,
                    matrix,
                    interfaceBouCoeffsData,
                    interfaceIntCoeffsData,
                    coeffsOffsets,
                    interfaces,
                    solverControls
                )
            );


        // Constructors

            //- Construct solver for given field name, matrix, flat coeffs etc
            solver
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const scalarField& interfaceBouCoeffsData,
                const scalarField& interfaceIntCoeffsData,
                const labelUList& coeffsOffsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& solverControls
            );

        public:

        // Selectors

            //- Return a new solver with flat coefficient storage
            static autoPtr<solver> New
            (
                const word& solverName,
                const word& fieldName,
                const lduMatrix& matrix,
                const scalarField& interfaceBouCoeffsData,
                const scalarField& interfaceIntCoeffsData,
                const labelUList& coeffsOffsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& solverControls
            );

            //- Return a new solver from dictionary, flat coefficient storage
            static autoPtr<solver> New
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const scalarField& interfaceBouCoeffsData,
                const scalarField& interfaceIntCoeffsData,
                const labelUList& coeffsOffsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& solverControls
            );



        //- Destructor
        virtual ~solver() = default;


        // Member Functions

            const word& fieldName() const noexcept
            {
                return fieldName_;
            }

            const lduMatrix& matrix() const noexcept
            {
                return matrix_;
            }

            //- Return flat storage for boundary coefficients
            const scalarField& interfaceBouCoeffsData() const noexcept
            {
                return interfaceBouCoeffsData_;
            }

            //- Return flat storage for internal coefficients
            const scalarField& interfaceIntCoeffsData() const noexcept
            {
                return interfaceIntCoeffsData_;
            }

            //- Return per-interface offsets into flat coefficient storage
            const labelList& coeffsOffsets() const noexcept
            {
                return coeffsOffsets_;
            }

            //- Return SubField view of boundary coefficients for interface i
            SubField<scalar> interfaceBouCoeffs(const label inti) const
            {
                return SubField<scalar>
                (
                    interfaceBouCoeffsData_,
                    coeffsOffsets_[inti+1] - coeffsOffsets_[inti],
                    coeffsOffsets_[inti]
                );
            }

            //- Return SubField view of internal coefficients for interface i
            SubField<scalar> interfaceIntCoeffs(const label inti) const
            {
                return SubField<scalar>
                (
                    interfaceIntCoeffsData_,
                    coeffsOffsets_[inti+1] - coeffsOffsets_[inti],
                    coeffsOffsets_[inti]
                );
            }

            const lduInterfaceFieldPtrsList& interfaces() const noexcept
            {
                return interfaces_;
            }


            //- Read and reset the solver parameters from the given stream
            virtual void read(const dictionary&);

            //- Solve with given field and rhs
            virtual solverPerformance solve
            (
                scalarField& psi,
                const scalarField& source,
                const direction cmpt=0
            ) const = 0;

            //- Solve with given field and rhs (in solveScalar precision).
            //  Default is to call solve routine
            virtual solverPerformance scalarSolve
            (
                solveScalarField& psi,
                const solveScalarField& source,
                const direction cmpt=0
            ) const;

            //- Return the matrix norm using the specified norm method
            solveScalarField::cmptType normFactor
            (
                const solveScalarField& psi,
                const solveScalarField& source,
                const solveScalarField& Apsi,
                solveScalarField& tmpField,
                const lduMatrix::normTypes normType
            ) const;

            //- Return the matrix norm used to normalise the residual for the
            //- stopping criterion
            solveScalarField::cmptType normFactor
            (
                const solveScalarField& psi,
                const solveScalarField& source,
                const solveScalarField& Apsi,
                solveScalarField& tmpField
            ) const
            {
                return this->normFactor(psi, source, Apsi, tmpField, normType_);
            }
    };


    // -----------------------------------------------------------------------
    //- Abstract base-class for lduMatrix smoothers
    class smoother
    {
    protected:

        // Protected Data

            word fieldName_;
            const lduMatrix& matrix_;
            const scalarField& interfaceBouCoeffsData_;
            const scalarField& interfaceIntCoeffsData_;
            const labelUList& coeffsOffsets_;
            const lduInterfaceFieldPtrsList& interfaces_;


    public:

        //- Find the smoother name (directly or from a sub-dictionary)
        static word getName(const dictionary&);

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                symMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const scalarField& interfaceBouCoeffsData,
                    const scalarField& interfaceIntCoeffsData,
                    const labelUList& coeffsOffsets,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& solverControls
                ),
                (
                    fieldName,
                    matrix,
                    interfaceBouCoeffsData,
                    interfaceIntCoeffsData,
                    coeffsOffsets,
                    interfaces,
                    solverControls
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                asymMatrix,
                (
                    const word& fieldName,
                    const lduMatrix& matrix,
                    const scalarField& interfaceBouCoeffsData,
                    const scalarField& interfaceIntCoeffsData,
                    const labelUList& coeffsOffsets,
                    const lduInterfaceFieldPtrsList& interfaces,
                    const dictionary& solverControls
                ),
                (
                    fieldName,
                    matrix,
                    interfaceBouCoeffsData,
                    interfaceIntCoeffsData,
                    coeffsOffsets,
                    interfaces,
                    solverControls
                )
            );


        // Constructors

            //- Construct for given field name, matrix etc (flat storage)
            smoother
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const scalarField& interfaceBouCoeffsData,
                const scalarField& interfaceIntCoeffsData,
                const labelUList& coeffsOffsets,
                const lduInterfaceFieldPtrsList& interfaces
            );


        // Selectors

            //- Return a new smoother
            static autoPtr<smoother> New
            (
                const word& fieldName,
                const lduMatrix& matrix,
                const scalarField& interfaceBouCoeffsData,
                const scalarField& interfaceIntCoeffsData,
                const labelUList& coeffsOffsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const dictionary& solverControls
            );


        //- Destructor
        virtual ~smoother() = default;


        // Member Functions

            const word& fieldName() const noexcept
            {
                return fieldName_;
            }

            const lduMatrix& matrix() const noexcept
            {
                return matrix_;
            }

            //- Return flat storage for boundary coefficients
            const scalarField& interfaceBouCoeffsData() const noexcept
            {
                return interfaceBouCoeffsData_;
            }

            //- Return flat storage for internal coefficients
            const scalarField& interfaceIntCoeffsData() const noexcept
            {
                return interfaceIntCoeffsData_;
            }

            //- Return per-interface offsets into flat coefficient storage
            const labelList& coeffsOffsets() const noexcept
            {
                return coeffsOffsets_;
            }

            //- Convenience accessor (returns flat data)
            const scalarField& interfaceBouCoeffs() const noexcept
            {
                return interfaceBouCoeffsData_;
            }

            //- Convenience accessor (returns flat data)
            const scalarField& interfaceIntCoeffs() const noexcept
            {
                return interfaceIntCoeffsData_;
            }

            const lduInterfaceFieldPtrsList& interfaces() const noexcept
            {
                return interfaces_;
            }


            //- Smooth the solution for a given number of sweeps
            virtual void smooth
            (
                solveScalarField& psi,
                const scalarField& source,
                const direction cmpt,
                const label nSweeps
            ) const = 0;

            //- Smooth the solution for a given number of sweeps
            virtual void scalarSmooth
            (
                solveScalarField& psi,
                const solveScalarField& source,
                const direction cmpt,
                const label nSweeps
            ) const = 0;
    };


    // -----------------------------------------------------------------------
    //- Abstract base-class for lduMatrix preconditioners
    class preconditioner
    {
    protected:

        // Protected Data

            //- Reference to the base-solver this preconditioner is used with
            const solver& solver_;


    public:

        //- Find the preconditioner name (directly or from a sub-dictionary)
        static word getName(const dictionary&);

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                symMatrix,
                (
                    const solver& sol,
                    const dictionary& solverControls
                ),
                (sol, solverControls)
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                asymMatrix,
                (
                    const solver& sol,
                    const dictionary& solverControls
                ),
                (sol, solverControls)
            );


        // Constructors

            //- Construct for given solver
            explicit preconditioner(const solver& sol)
            :
                solver_(sol)
            {}


        // Selectors

            //- Return a new preconditioner
            static autoPtr<preconditioner> New
            (
                const solver& sol,
                const dictionary& solverControls
            );


        //- Destructor
        virtual ~preconditioner() = default;


        // Member Functions

            //- Read and reset the preconditioner parameters
            //- from the given stream
            virtual void read(const dictionary&)
            {}

            //- Return wA the preconditioned form of residual rA
            virtual void precondition
            (
                solveScalarField& wA,
                const solveScalarField& rA,
                const direction cmpt=0
            ) const = 0;

            //- Return wT the transpose-matrix preconditioned form of
            //- residual rT.
            //  This is only required for preconditioning asymmetric matrices.
            virtual void preconditionT
            (
                solveScalarField& wT,
                const solveScalarField& rT,
                const direction cmpt=0
            ) const
            {
                NotImplemented;
            }

            //- Signal end of solver
            virtual void setFinished(const solverPerformance& perf) const
            {}
    };


    // -----------------------------------------------------------------------

    // Static Data

        // Declare name of the class and its debug switch
        ClassName("lduMatrix");


    // Constructors

        //- Construct (without coefficients) for an LDU addressed mesh.
        //  Not yet 'explicit' (legacy code may rely on implicit construct)
        lduMatrix(const lduMesh& mesh);

        //- Copy construct
        lduMatrix(const lduMatrix&);

        //- Move construct
        lduMatrix(lduMatrix&&);

        //- Construct as copy or re-use as specified.
        lduMatrix(lduMatrix&, bool reuse);

        //- Construct given an LDU addressed mesh and an Istream
        //- from which the coefficients are read
        lduMatrix(const lduMesh& mesh, Istream& is);


    //- Destructor
    ~lduMatrix() = default;


    // Member Functions

    // Addressing

        //- Return the LDU mesh from which the addressing is obtained
        const lduMesh& mesh() const noexcept
        {
            return lduMesh_;
        }

        //- Set the LDU mesh containing the addressing
        void setLduMesh(const lduMesh& m)
        {
            lduMesh_ = m;
        }

        //- Return the LDU addressing
        const lduAddressing& lduAddr() const
        {
            return mesh().lduAddr();
        }

        //- Return the patch evaluation schedule
        const lduSchedule& patchSchedule() const
        {
            return mesh().lduAddr().patchSchedule();
        }


    // Coefficients

        const scalarField& diag() const;
        const scalarField& upper() const;
        const scalarField& lower() const;

        scalarField& diag();
        scalarField& upper();
        scalarField& lower();

        // Size with externally provided sizes
        // (for constructing with 'fake' mesh in GAMG)

        scalarField& diag(label size);
        scalarField& upper(label nCoeffs);
        scalarField& lower(label nCoeffs);


    // CSR

        bool hasLowerCSR() const noexcept { return bool(lowerCSRPtr_); }

        const scalarField& lowerCSR() const;

        scalarField& lowerCSR();

        //- Work array
        const solveScalarField& work() const;

        //- Work array
        solveScalarField& work(label size) const;


    // Characteristics

        //- The matrix type (empty, diagonal, symmetric, ...)
        word matrixTypeName() const;

        bool hasDiag() const noexcept { return bool(diagPtr_); }
        bool hasUpper() const noexcept { return bool(upperPtr_); }
        bool hasLower() const noexcept { return bool(lowerPtr_); }

        //- Matrix has diagonal only
        bool diagonal() const noexcept
        {
            return (diagPtr_ && !lowerPtr_ && !upperPtr_);
        }

        //- Matrix is symmetric
        bool symmetric() const noexcept
        {
            return (diagPtr_ && !lowerPtr_ && upperPtr_);
        }

        //- Matrix is asymmetric (ie, full)
        bool asymmetric() const noexcept
        {
            return (diagPtr_ && lowerPtr_ && upperPtr_);
        }


        // Operations

            void sumDiag();
            void negSumDiag();

            void sumMagOffDiag(scalarField& sumOff) const;

            //- Matrix multiplication with flat field + offset array
            void Amul
            (
                solveScalarField&,
                const tmp<solveScalarField>&,
                const solveScalarField& interfaceBouCoeffsData,
                const labelUList& offsets,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;

            //- Matrix transpose multiplication with flat field + offset array
            void Tmul
            (
                solveScalarField&,
                const tmp<solveScalarField>&,
                const solveScalarField& interfaceIntCoeffsData,
                const labelUList& offsets,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;


            //- Sum the coefficients on each row of the matrix
            void sumA
            (
                solveScalarField&,
                const solveScalarField& interfaceBouCoeffsData,
                const labelUList& offsets,
                const lduInterfaceFieldPtrsList&
            ) const;


            void residual
            (
                solveScalarField& rA,
                const solveScalarField& psi,
                const scalarField& source,
                const solveScalarField& interfaceBouCoeffsData,
                const labelUList& offsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;

            //- Residual with flat field + offset array
            tmp<solveScalarField> residual
            (
                const solveScalarField& psi,
                const scalarField& source,
                const solveScalarField& interfaceBouCoeffsData,
                const labelUList& offsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;


            //- Initialise interfaced interfaces for matrix operations
            //- using flat field + offset array
            void initMatrixInterfaces
            (
                const bool add,
                const solveScalarField& interfaceCoeffsData,
                const labelUList& offsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const solveScalarField& psiif,
                solveScalarField& result,
                const direction cmpt
            ) const;

            //- Update interfaced interfaces for matrix operations
            //- using flat field + offset array
            void updateMatrixInterfaces
            (
                const bool add,
                const solveScalarField& interfaceCoeffsData,
                const labelUList& offsets,
                const lduInterfaceFieldPtrsList& interfaces,
                const solveScalarField& psiif,
                solveScalarField& result,
                const direction cmpt,
                const int startRequest  //!< starting request (non-blocking)
            ) const;

            //- Set the residual field using an IOField on the object registry
            //- if it exists
            void setResidualField
            (
                const scalarField& residual,
                const word& fieldName,
                const bool initial
            ) const;

            template<class Type>
            tmp<Field<Type>> H(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type>> H(const tmp<Field<Type>>&) const;

            tmp<scalarField> H1() const;

            template<class Type>
            tmp<Field<Type>> faceH(const Field<Type>&) const;

            template<class Type>
            tmp<Field<Type>> faceH(const tmp<Field<Type>>&) const;


    // Info

        //- Return info proxy,
        //- used to print matrix information to a stream
        InfoProxy<lduMatrix> info() const noexcept { return *this; }


    // Member Operators

        //- Copy assignment
        void operator=(const lduMatrix&);

        //- Move assignment
        void operator=(lduMatrix&&);

        void negate();

        void operator+=(const lduMatrix&);
        void operator-=(const lduMatrix&);

        void operator*=(const scalarField&);
        void operator*=(scalar);


    // Ostream Operators

        friend Ostream& operator<<(Ostream&, const lduMatrix&);
        friend Ostream& operator<<(Ostream&, const InfoProxy<lduMatrix>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lduMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
