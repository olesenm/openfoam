/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 M. Janssens
    Copyright (C) 2025 OpenCFD Ltd.
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

InNamespace
    Foam::Expression

Description
    Expression templates for fvMatrix.

    TBD:
    - useImplicit_ etc.
    - lduMatrix::lowerCSR ?

SourceFiles
    fvMatrixExpression.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMatrixExpression_H
#define Foam_fvMatrixExpression_H

#include "ListExpression.H"
#include "GeometricFieldExpression.H"
#include "fvMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class fvMatrix;

namespace Expression
{

/*---------------------------------------------------------------------------*\
                     Class fvMatrixExpression Declaration
\*---------------------------------------------------------------------------*/

//- Expression of fvMatrix
//
// Bit different:
// - immediately stores state (hasDiag, hasUpper etc.)
// - delays execution of building coefficients
// - big problem: individual components can be built differently so
//   different template argument ... Needs to be templated on all outputs
//   separately!
template
<
    typename E,
    typename DiagExpr,
    typename UpperExpr,
    typename LowerExpr,
    typename FaceFluxExpr,
    typename SourceExpr
>
class fvMatrixExpression
{
protected:

    const bool hasDiag_;
    const bool hasUpper_;
    const bool hasLower_;
    const bool hasFaceFluxCorrection_;
    dimensionSet dimensions_;

public:
    static constexpr bool is_leaf = false;

    // Constructors

        fvMatrixExpression
        (
            const bool hasDiag,
            const bool hasUpper,
            const bool hasLower,
            const bool hasFaceFluxCorrection,
            const dimensionSet& dimensions
        )
        :
            hasDiag_(hasDiag),
            hasUpper_(hasUpper),
            hasLower_(hasLower),
            hasFaceFluxCorrection_(hasFaceFluxCorrection),
            dimensions_(dimensions)
        {
            #ifdef DEBUG_ET
            Pout<< "fvMatrixExpression :" << nl
                << "    E              : " << error::demangle<E>() << nl;
            if (hasDiag_)
            {
                Pout<< "    DiagExpr       : "
                    << error::demangle<DiagExpr>() << nl;
            }
            {
                Pout<< "    SourceExpr     : "
                    << error::demangle<SourceExpr>() << nl;
            }
            if (hasUpper_)
            {
                Pout<< "    UpperExpr      : "
                    << error::demangle<UpperExpr>() << nl;
            }
            if (hasLower_)
            {
                Pout<< "    LowerExpr      : "
                    << error::demangle<LowerExpr>() << nl;
            }
            if (hasFaceFluxCorrection_)
            {
                Pout<< "    FaceFluxExpr   : "
                    << error::demangle<FaceFluxExpr>() << nl;
            }
            #endif
        }


    // lduMatrix part

        DiagExpr diag() const
        {
            return static_cast<const E&>(*this).diag();
        }

        DiagExpr diag()
        {
            return static_cast<E&>(*this).diag();
        }

        UpperExpr upper() const
        {
            return static_cast<const E&>(*this).upper();
        }

        UpperExpr upper()
        {
            return static_cast<E&>(*this).upper();
        }

        LowerExpr lower() const
        {
            return static_cast<const E&>(*this).lower();
        }

        LowerExpr lower()
        {
            return static_cast<E&>(*this).lower();
        }

        bool hasDiag() const
        {
             return hasDiag_;
        }
        bool hasUpper() const
        {
             return hasUpper_;
        }
        bool hasLower() const
        {
             return hasLower_;
        }
        bool diagonal() const
        {
             return (hasDiag() && !hasLower() && !hasUpper());
        }
        bool symmetric() const
        {
             return (hasDiag() && !hasLower() && hasUpper());
        }
        bool asymmetric() const
        {
             return (hasDiag() && hasLower() && hasUpper());
        }
        bool hasFaceFluxCorrection() const
        {
             return hasFaceFluxCorrection_;
        }


    // fvMatrix

        SourceExpr source()
        {
            return static_cast<E&>(*this).source();
        }

        SourceExpr source() const
        {
            return static_cast<const E&>(*this).source();
        }

        UpperExpr internalCoeffs(const label i)
        {
            return static_cast<E&>(*this).internalCoeffs(i);
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return static_cast<const E&>(*this).internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i)
        {
            return static_cast<E&>(*this).boundaryCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return static_cast<const E&>(*this).boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return static_cast<const E&>(*this).faceFluxCorrection();
        }

        FaceFluxExpr faceFluxCorrection()
        {
            return static_cast<E&>(*this).faceFluxCorrection();
        }

        const auto& dimensions() const noexcept
        {
            return dimensions_;
        }

        auto& dimensions() noexcept
        {
            return dimensions_;
        }


    // Other

        //- Helper to evaluate (=construct) an fvMatrix
        template<class Matrix>
        Matrix& evaluate(Matrix& m) const
        {
            {
                // Assign to source.
                source().evaluate(m.source());
            }
            if (hasDiag())
            {
                // Assign to diag. Also sets hasDiag.
                diag().evaluate(m.diag());
            }
            if (hasUpper())
            {
                upper().evaluate(m.upper());
            }
            if (hasLower())
            {
                lower().evaluate(m.lower());
            }
            if (hasFaceFluxCorrection())
            {
                faceFluxCorrection().evaluate(*m.faceFluxCorrectionPtr());
            }

            // Do boundary
            if (hasUpper() || hasLower())
            {
                auto& intCoeffs = m.internalCoeffs();
                auto& bouCoeffs = m.boundaryCoeffs();
                const label n = m.psi().boundaryField().size();

                typedef
                    typename std::remove_reference_t<decltype(m.source())>::
                    value_type value_type;

#if defined(FOAM_OFFLOAD)
                constexpr bool foam_offload = true;
#else
                constexpr bool foam_offload = false;
#endif
                if constexpr (foam_offload)
                {
                    //- Parallel via ListsRefWrap. See GeometricFieldExpression.
                    //  Flat storage: internalCoeffs()/boundaryCoeffs() are
                    //  single contiguous fields, so evaluate the concatenation
                    //  of all per-patch expressions into them in one pass.
                    std::vector<UpperExpr> intExprs;
                    std::vector<UpperExpr> bouExprs;
                    intExprs.reserve(n);
                    bouExprs.reserve(n);
                    for (label i = 0; i < n; ++i)
                    {
                        intExprs.emplace_back(internalCoeffs(i));
                        bouExprs.emplace_back(boundaryCoeffs(i));
                    }

                    // Internal : evaluate all patches at once
                    ListsRefWrap<UpperExpr> intSrc(intExprs);
                    ListRefWrap<value_type>(intCoeffs, intSrc);

                    // Boundary : evaluate all patches at once
                    ListsRefWrap<UpperExpr> bouSrc(bouExprs);
                    ListRefWrap<value_type>(bouCoeffs, bouSrc);
                }
                else
                {
                    //- Sequential
                    for (label i = 0; i < n; ++i)
                    {
                        const auto intExpr = internalCoeffs(i);
                        if (intExpr.size())
                        {
                            auto intSlice = m.internalCoeffs(i);
                            intExpr.fill(intSlice);
                            const auto bouExpr = boundaryCoeffs(i);
                            auto bouSlice = m.boundaryCoeffs(i);
                            bouExpr.fill(bouSlice);
                        }
                    }
                }
            }

            const_cast<dimensionSet&>(m.dimensions()) = dimensions();

            return m;
        }
};


/*---------------------------------------------------------------------------*\
                       Class fvMatrixRefWrap Declaration
\*---------------------------------------------------------------------------*/

//- Expression wrap of non-const reference to fvMatrix
template<class Matrix>
class fvMatrixRefWrap
:
    public fvMatrixExpression
    <
        fvMatrixRefWrap<Matrix>,
        ListRefWrap<typename Matrix::psiFieldType::value_type>, // diag
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // upper
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // lower
        GeometricFieldRefWrap<typename Matrix::faceFluxFieldType>,// fflux
        ListRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef Matrix this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef ListRefWrap<value_type> DiagExpr;
    typedef ListConstRefWrap<value_type> UpperExpr;
    typedef ListConstRefWrap<value_type> LowerExpr;
    typedef GeometricFieldRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListRefWrap<value_type> SourceExpr;


private:

    this_type& elems_;


public:

    //- Copy construct
    fvMatrixRefWrap(this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}

    //- Move construct
    fvMatrixRefWrap(this_type&& elems)
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}

    // Construct from Matrix, forcing its evaluation.
    template<typename E>
    fvMatrixRefWrap
    (
        this_type& elems,
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            expr.hasDiag(),
            expr.hasUpper(),
            expr.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {
        expr.evaluate(elems_);
    }

    //- Assignment
    template<typename E>
    void operator=
    (
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    {
        expr.evaluate(elems_);
    }

    //- Evaluate and return as GeoField. Rename to evaluate to make it clear
    //- it takes time? Or leave as indexing for convenience?
    template<typename E>
    this_type& evaluate
    (
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    {
        return expr.evaluate(elems_);
    }


    // lduMatrix part

        DiagExpr diag() const
        {
            return elems_.diag();
        }

        DiagExpr diag()
        {
            return elems_.diag();
        }

        UpperExpr upper() const
        {
            return elems_.upper();
        }

        UpperExpr upper()
        {
            return elems_.upper();
        }

        LowerExpr lower() const
        {
            return elems_.lower();
        }

        LowerExpr lower()
        {
            return elems_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_.source();
        }

        SourceExpr source()
        {
            return elems_.source();
        }

        UpperExpr internalCoeffs(const label patchi)
        {
            return elems_.internalCoeffs(patchi);
        }

        UpperExpr internalCoeffs(const label patchi) const
        {
            return elems_.internalCoeffs(patchi);
        }

        UpperExpr boundaryCoeffs(const label patchi)
        {
            return elems_.boundaryCoeffs(patchi);
        }

        UpperExpr boundaryCoeffs(const label patchi) const
        {
            return elems_.boundaryCoeffs(patchi);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return *const_cast<this_type&>
            (
                elems_
            ).faceFluxCorrectionPtr().get();
        }

        FaceFluxExpr faceFluxCorrection()
        {
            return *elems_.faceFluxCorrectionPtr().get();
        }
};


/*---------------------------------------------------------------------------*\
                    Class fvMatrixConstRefWrap Declaration
\*---------------------------------------------------------------------------*/

//- Expression wrap of const reference to fvMatrix
template<class Matrix>
class fvMatrixConstRefWrap
:
    public fvMatrixExpression
    <
        fvMatrixConstRefWrap<Matrix>,
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // diag
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // upper
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // lower
        GeometricFieldConstRefWrap<typename Matrix::faceFluxFieldType>, // fflux
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef Matrix this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef ListConstRefWrap<value_type> DiagExpr;
    typedef ListConstRefWrap<value_type> UpperExpr;
    typedef ListConstRefWrap<value_type> LowerExpr;
    typedef GeometricFieldConstRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListConstRefWrap<value_type> SourceExpr;


private:

    const this_type& elems_;


public:

    // Construct from components
    fvMatrixConstRefWrap(const this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixConstRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return elems_.diag();
        }

        UpperExpr upper() const
        {
            return elems_.upper();
        }

        LowerExpr lower() const
        {
            return elems_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_.source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return elems_.internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return elems_.boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            typename Matrix::faceFluxFieldType& ffc =
                *const_cast<this_type&>(elems_).faceFluxCorrectionPtr();
            return ffc;
        }
};


/*---------------------------------------------------------------------------*\
                    Class fvMatrixConstTmpWrap Declaration
\*---------------------------------------------------------------------------*/

//- Expression wrap of const reference to fvMatrix
template<class Matrix>
class fvMatrixConstTmpWrap
:
    public fvMatrixExpression
    <
        fvMatrixConstTmpWrap<Matrix>,
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // diag
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // upper
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // lower
        GeometricFieldConstRefWrap<typename Matrix::faceFluxFieldType>, // fflux
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef tmp<Matrix> this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef ListConstRefWrap<value_type> DiagExpr;
    typedef ListConstRefWrap<value_type> UpperExpr;
    typedef ListConstRefWrap<value_type> LowerExpr;
    typedef GeometricFieldConstRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListConstRefWrap<value_type> SourceExpr;


private:

    const this_type elems_;


public:

    // Construct from components
    fvMatrixConstTmpWrap(const this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixConstTmpWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems().hasDiag(),
            elems().hasUpper(),
            elems().hasLower(),
            elems().hasFaceFluxCorrection(),
            elems().dimensions()
        ),
        elems_(elems)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return elems_().diag();
        }

        UpperExpr upper() const
        {
            return elems_().upper();
        }

        LowerExpr lower() const
        {
            return elems_().lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_().source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return elems_().internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return elems_().boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return *const_cast<Matrix&>
            (
                elems_()
            ).faceFluxCorrectionPtr();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Expressions on fvMatrices
// ~~~~~~~~~~~~~~~~~~~~~~~~~

// operator+
template<typename E1, typename E2>
class fvm_add
:
    public fvMatrixExpression
    <
        fvm_add<E1, E2>,
        List_add                    // DiagExpr
        <
            typename E1::DiagExpr,
            typename E2::DiagExpr
        >,
        List_add                    // UpperExpr
        <
            typename E1::UpperExpr,
            typename E2::UpperExpr
        >,
        List_add                    // LowerExpr
        <
            typename E1::LowerExpr,
            typename E2::LowerExpr
        >,
        GF_add                    // FaceFluxExpr
        <
            typename E1::FaceFluxExpr,
            typename E2::FaceFluxExpr
        >,
        List_add                    // SourceExpr
        <
            typename E1::SourceExpr,
            typename E2::SourceExpr
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_add
    <
        typename E1::DiagExpr,
        typename E2::DiagExpr
    > DiagExpr;
    typedef List_add
    <
        typename E1::UpperExpr,
        typename E2::UpperExpr
    > UpperExpr;
    typedef List_add
    <
        typename E1::LowerExpr,
        typename E2::LowerExpr
    > LowerExpr;
    typedef GF_add
    <
        typename E1::FaceFluxExpr,
        typename E2::FaceFluxExpr
    > FaceFluxExpr;
    typedef List_add
    <
        typename E1::SourceExpr,
        typename E2::SourceExpr
    > SourceExpr;

    fvm_add(const E1& u, const E2& v)
    :
        fvMatrixExpression
        <
            fvm_add<E1, E2>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag() || v.hasDiag(),
            u.hasUpper() || v.hasUpper(),
            u.hasLower() || v.hasLower(),
            u.hasFaceFluxCorrection() || v.hasFaceFluxCorrection(),
            u.dimensions()
        ),
        u_(u),
        v_(v)
    {
        //Pout<< "fvm_add :" << nl
        //    << "    u_ : " << error::demangle(typeid(u_)) << nl
        //    << "    v_ : " << error::demangle(typeid(v_)) << endl;
    }


    // lduMatrix part

        DiagExpr diag() const
        {
            //return u_.diag() + v_.diag();
            return DiagExpr(u_.diag(), v_.diag());
        }

        UpperExpr upper() const
        {
            return UpperExpr(u_.upper(), v_.upper());
        }

        LowerExpr lower() const
        {
            return LowerExpr(u_.lower(), v_.lower());
        }

    // fvMatrix

        SourceExpr source() const
        {
            return SourceExpr(u_.source(), v_.source());
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return UpperExpr(u_.boundaryCoeffs(i), v_.boundaryCoeffs(i));
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return UpperExpr(u_.internalCoeffs(i), v_.internalCoeffs(i));
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return FaceFluxExpr(u_.faceFluxCorrection(), v_.faceFluxCorrection());
        }
};
template<typename E1, typename E2>
fvm_add<E1, E2>
operator+
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    return fvm_add<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}


// operator-, operator==
template<typename E1, typename E2>
class fvm_subtract
:
    public fvMatrixExpression
    <
        fvm_subtract<E1, E2>,
        List_subtract                   // diag
        <
            typename E1::DiagExpr,
            typename E2::DiagExpr
        >,
        List_subtract                   // upper
        <
            typename E1::UpperExpr,
            typename E2::UpperExpr
        >,
        List_subtract                   // lower
        <
            typename E1::LowerExpr,
            typename E2::LowerExpr
        >,
        GF_subtract                   // faceflux
        <
            typename E1::FaceFluxExpr,
            typename E2::FaceFluxExpr
        >,
        List_subtract                   // source
        <
            typename E1::SourceExpr,
            typename E2::SourceExpr
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_subtract
    <
        typename E1::DiagExpr,
        typename E2::DiagExpr
    > DiagExpr;
    typedef List_subtract
    <
        typename E1::UpperExpr,
        typename E2::UpperExpr
    > UpperExpr;
    typedef List_subtract
    <
        typename E1::LowerExpr,
        typename E2::LowerExpr
    > LowerExpr;
    typedef GF_subtract
    <
        typename E1::FaceFluxExpr,
        typename E2::FaceFluxExpr
    > FaceFluxExpr;
    typedef List_subtract
    <
        typename E1::SourceExpr,
        typename E2::SourceExpr
    > SourceExpr;

    fvm_subtract(const E1& u, const E2& v)
    :
        fvMatrixExpression
        <
            fvm_subtract<E1, E2>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag() || v.hasDiag(),
            u.hasUpper() || v.hasUpper(),
            u.hasLower() || v.hasLower(),
            u.hasFaceFluxCorrection() || v.hasFaceFluxCorrection(),
            u.dimensions()
        ),
        u_(u),
        v_(v)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return u_.diag() - v_.diag();
        }

        UpperExpr upper() const
        {
            return u_.upper() - v_.upper();
        }

        LowerExpr lower() const
        {
            return u_.lower() - v_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return u_.source() - v_.source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return u_.internalCoeffs(i) - v_.internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return u_.boundaryCoeffs(i) - v_.boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return u_.faceFluxCorrection() - v_.faceFluxCorrection();
        }
};
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator-
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    return fvm_subtract<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator==
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    // workaround: Su,Sp etc return dimless
    if (v.dimensions() == dimless)
    {
        const_cast<dimensionSet&>(v.dimensions()).reset(u.dimensions());
    }
    else if (u.dimensions() == dimless)
    {
        const_cast<dimensionSet&>(u.dimensions()).reset(v.dimensions());
    }
    return fvm_subtract<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}


// Negation : operator-
template<typename E1>
class fvm_negate
:
    public fvMatrixExpression
    <
        fvm_negate<E1>,
        List_negate<typename E1::DiagExpr>,
        List_negate<typename E1::UpperExpr>,
        List_negate<typename E1::LowerExpr>,
        GF_negate<typename E1::FaceFluxExpr>,
        List_negate<typename E1::SourceExpr>
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_negate<typename E1::DiagExpr> DiagExpr;
    typedef List_negate<typename E1::UpperExpr> UpperExpr;
    typedef List_negate<typename E1::LowerExpr> LowerExpr;
    typedef GF_negate<typename E1::FaceFluxExpr> FaceFluxExpr;
    typedef List_negate<typename E1::SourceExpr> SourceExpr;

    fvm_negate(const E1& u)
    :
        fvMatrixExpression
        <
            fvm_negate<E1>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag(),
            u.hasUpper(),
            u.hasLower(),
            u.hasFaceFluxCorrection(),
            u.dimensions()
        ),
        u_(u)
    {}


    // lduMatrix part

        // Note: conversion constructor from diag to List_negate<> takes
        // care of negation
        DiagExpr diag() const
        {
            return u_.diag();
        }

        UpperExpr upper() const
        {
            return u_.upper();
        }

        LowerExpr lower() const
        {
            return u_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return u_.source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return u_.internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return u_.boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return u_.faceFluxCorrection();
        }
};
template<typename E1>
auto operator-
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u
)
{
    return fvm_negate<E1>(static_cast<const E1&>(u));
}


// Source manipulation
// ~~~~~~~~~~~~~~~~~~~

//- fvMatrix, internal parts and dimensions only. No upper/lower/boundaryCoeffs
//  etc. Used for source manipulation.
template<class E1, class E2>
class fvMatrixInternal
:
    public fvMatrixExpression
    <
        fvMatrixInternal<E1, E2>,
        E1, // Diag
        UniformListWrap<typename E1::value_type>, // Upper, Not used
        UniformListWrap<typename E1::value_type>, // Lower, Not used
        GeometricFieldConstRefWrap
        <
            GeometricField<typename E1::value_type, fvsPatchField, surfaceMesh>
        >,  // FaceFlux, Not used
        E2  // Source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef E1 DiagExpr;

    // typedef ListConstRefWrap<value_type> UpperExpr;
    typedef UniformListWrap<value_type> UpperExpr;
    // typedef ListConstRefWrap<value_type> LowerExpr;
    typedef UniformListWrap<value_type> LowerExpr;

    typedef GeometricField<value_type, fvsPatchField, surfaceMesh>
            faceFluxFieldType;
    typedef GeometricFieldConstRefWrap<faceFluxFieldType> FaceFluxExpr;
    typedef E2 SourceExpr;


private:

    //- Underlying matrix - used for sizes only
    const fvMatrix<value_type>& m_;

    const E1 diag_;
    const SourceExpr source_;


public:

    // Construct from component expressions
    fvMatrixInternal
    (
        const fvMatrix<value_type>& m,
        const DiagExpr& diag,
        const SourceExpr& source
    )
    :
        fvMatrixExpression
        <
            fvMatrixInternal<DiagExpr, SourceExpr>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            FaceFluxExpr,    // face flux not used
            SourceExpr
        >
        (
            true,       //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            m.dimensions()  //elems.dimensions()
        ),
        m_(m),
        diag_(diag),
        source_(source)
    {}

    // Construct from component expressions
    fvMatrixInternal
    (
        const fvMatrix<value_type>& m,
        const DiagExpr&& diag,
        const SourceExpr&& source
    )
    :
        fvMatrixExpression
        <
            fvMatrixInternal<DiagExpr, SourceExpr>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            FaceFluxExpr,       // faceflux not yet used
            SourceExpr
        >
        (
            true,       //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            m.dimensions()  //elems.dimensions()
        ),
        m_(m),
        diag_(diag),
        source_(source)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return diag_;
        }

        UpperExpr upper() const
        {
            // return List<value_type>::null();
            return UpperExpr(m_.upper().size(), value_type(0));
        }

        LowerExpr lower() const
        {
            // return List<value_type>::null();
            return LowerExpr(m_.lower().size(), value_type(0));
        }


    // fvMatrix

        SourceExpr source() const
        {
            return source_;
        }

        UpperExpr internalCoeffs(const label i) const
        {
            // return List<value_type>::null();
            return UpperExpr
            (
                m_.internalCoeffs(i).size(),
                value_type(0)
            );
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            // return List<value_type>::null();
            return UpperExpr
            (
                m_.boundaryCoeffs(i).size(),
                value_type(0)
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            //return faceFluxFieldType::null();
            return
            (
                *const_cast<fvMatrix<value_type>&>(m_).
                    faceFluxCorrectionPtr().get()
            );
        }
};

//- Sources using matrix - provides dimensions
template<typename E1>
auto SuSp
(
    const E1& expression,
    const fvMatrix<typename E1::value_type>& m
)
{
    // Equivalent of fvm::SuSp

    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename E1::value_type> expr;

    // Wrap constants
    const constant constantZero(m.psi().size(), 0.0);
    const constant constantMinOne(m.psi().size(), -1.0);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(m.psi().mesh().V());
    return fvMatrixInternal
    (
        m,  // m.dimensions(),
        V*max(expression, constantZero),        // diag
        constantMinOne                          // source
       *V
       *min(expression, constantZero)
       *expr(m.psi().internalField())
    );
}
// Provide ListExpression for diag; dimensions from matrix
template<typename E1>
auto Sp
(
    const E1& expression,
    const fvMatrix<typename E1::value_type>& m
)
{
    // Equivalent of fvm::Sp

    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;

    // Wrap zero
    const constant constantZero(m.psi().size(), 0.0);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(m.psi().mesh().V());

    return fvMatrixInternal
    (
        m,              // m.dimensions(),
        V*expression,   // diag
        constantZero    // source
    );
}
template<typename E1>
auto Su
(
    const E1& expression,
    const fvMatrix<typename E1::value_type>& m
)
{
    // Equivalent of fvm::Su

    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;

    // Wrap zero
    const constant constantZero(m.psi().size(), 0.0);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(m.psi().mesh().V());

    return fvMatrixInternal
    (
        m,              // m.dimensions(),
        constantZero,   // Still provides an expression for diag. TBD.
        -V*expression
    );
}


// // Alternative in-place fvMatrix sources. Can be removed once ones above
// // work.
// template<class Expr>
// void rhs
// (
//     fvMatrix<typename Expr::value_type>& m,
//     const Expr& expression
// )
// {
//     //- Wrap of List as an expression
//     typedef ListConstRefWrap<typename Expr::value_type> expr;
//     //- Evaluator of an expression
//     typedef ListRefWrap<typename Expr::value_type> evaluator;

//     // Wrap mesh volume
//     const ListConstRefWrap<scalar> V(m.psi().mesh().V());

//     // Add expression to source
//     auto& s = m.source();
//     evaluator(s, expr(s) + V*expression);
// }
// template<class Expr>
// void rhs
// (
//     fvMatrix<typename Expr::value_type>& m,
//     const Expr& expression
// )
// {
//     Su(m, expression);
// }
// template<class Expr, class Expr2>
// void Sp
// (
//     fvMatrix<typename Expr::value_type>& m,
//     const Expr2& mult,
//     const Expr& expression
// )
// {
//     //- Wrap of List as an expression
//     typedef ListConstRefWrap<typename Expr::value_type> expr;
//     //- Evaluator of an expression
//     typedef ListRefWrap<typename Expr::value_type> evaluator;

//     // Wrap mesh volume
//     const ListConstRefWrap<scalar> V(m.psi().mesh().V());

//     // Add expression onto diag
//     auto& d = m.diag();
//     evaluator(d, expr(d) - mult*V*expression);
// }
// template<class Expr, class Expr2>
// void SuSp
// (
//     fvMatrix<typename Expr::value_type>& m,
//     const Expr2& mult,
//     const Expr& expression
// )
// {
//     //- Wrap of constant as a list expression
//     typedef UniformListWrap<scalar> constant;
//     //- Wrap of List as an expression
//     typedef ListConstRefWrap<typename Expr::value_type> expr;
//     //- Evaluator of an expression
//     typedef ListRefWrap<typename Expr::value_type> evaluator;

//     const constant constantZero(expression.size(), 0.0);

//     // Wrap mesh volume
//     const ListConstRefWrap<scalar> V(m.psi().mesh().V());

//     // Linearise expression
//     auto& diag = m.diag();
//     auto& source = m.source();

//     // diag() += domain*max(susp.field(), scalar(0));
//     evaluator
//     (
//         diag,
//         expr(diag) - mult*V*max(expression, constantZero)
//     );
//     // source() -= domain*min(susp.field(), scalar(0))*fld.primitiveField();
//     evaluator
//     (
//         source,
//         expr(source)
//       + (
//             mult
//            *V
//            *min(expression, constantZero)
//            *expr(m.psi().internalField())
//         )
//     );
// }


// //- Variant of fvMatrixInternal: return -1 size.
// template<class E1, class E2>
// class fvMatrixInternal2
// :
//     public fvMatrixExpression
//     <
//         fvMatrixInternal2<E1, E2>,
//         E1, // Diag
//         UniformListWrap<typename E1::value_type>, // Upper, Not used
//         UniformListWrap<typename E1::value_type>, // Lower, Not used
//         UniformGeometricFieldWrap2<typename E1::value_type>,  // FFlux, Not used
//         E2  // Source
//     >
// {
// public:

//     static constexpr bool is_leaf = false;  //true;

//     //- Type to return for internal field
//     typedef typename E1::value_type value_type;

//     //- Type to return for containers
//     typedef E1 DiagExpr;

//     typedef UniformListWrap<value_type> UpperExpr;
//     typedef UniformListWrap<value_type> LowerExpr;
//     typedef UniformGeometricFieldWrap2<typename E1::value_type> FaceFluxExpr;
//     typedef E2 SourceExpr;


// private:

//     const E1 diag_;
//     const SourceExpr source_;


// public:

//     // Construct from component expressions
//     fvMatrixInternal2
//     (
//         const GeometricField<value_type, fvPatchField, volMesh>& psi,
//         const DiagExpr& diag,
//         const SourceExpr& source
//     )
//     :
//         fvMatrixExpression
//         <
//             fvMatrixInternal2<DiagExpr, SourceExpr>,
//             DiagExpr,
//             UpperExpr,       // upper not yet used
//             LowerExpr,       // lower not yet used
//             FaceFluxExpr,    // face flux not used
//             SourceExpr
//         >
//         (
//             true,       //elems.hasDiag(),
//             false,      //elems.hasUpper(),
//             false,      //elems.hasLower(),
//             false,      //elems.hasFaceFluxCorrection(),
//             dimless     //elems.dimensions()
//         ),
//         diag_(diag),
//         source_(source)
//     {}


//     // lduMatrix part

//         DiagExpr diag() const
//         {
//             return diag_;
//         }

//         UpperExpr upper() const
//         {
//             return UpperExpr(-1, value_type(0));
//         }

//         LowerExpr lower() const
//         {
//             return LowerExpr(-1, value_type(0));
//         }


//     // fvMatrix

//         SourceExpr source() const
//         {
//             return source_;
//         }

//         UpperExpr internalCoeffs(const label i) const
//         {
//             return UpperExpr(-1, value_type(0));
//         }

//         UpperExpr boundaryCoeffs(const label i) const
//         {
//             return UpperExpr(-1, value_type(0));
//         }

//         FaceFluxExpr faceFluxCorrection() const
//         {
//             return FaceFluxExpr(value_type(0));
//         }
// };
// template<typename E1>
// auto SuSp
// (
//     const E1& expression,
//     const GeometricField<typename E1::value_type, fvPatchField, volMesh>& fld
// )
// {
//     //- Wrap of constant as a list expression
//     typedef Expression::UniformListWrap<scalar> constant;
//     //- Wrap of List as an expression
//     typedef ListConstRefWrap<typename E1::value_type> expr;

//     // Wrap constants
//     const constant constantZero(fld.size(), 0.0);
//     const constant constantMinOne(fld.size(), -1.0);

//     // Wrap mesh volume
//     const ListConstRefWrap<scalar> V(fld.mesh().V());

//     return fvMatrixInternal2
//     (
//         fld,
//         V*max(expression, constantZero),        // diag
//         V*min(expression, constantZero)*constantMinOne*expr(fld.internalField())
//     );
// }
// template<typename E1>
// auto Sp
// (
//     const E1& expression,
//     const GeometricField<typename E1::value_type, fvPatchField, volMesh>& fld
// )
// {
//     // Wrap zero
//     const Expression::UniformListWrap<scalar> constantZero(fld.size(), 0.0);

//     // Wrap mesh volume
//     const ListConstRefWrap<scalar> V(fld.mesh().V());

//     return fvMatrixInternal2(fld, V*expression, constantZero);
// }
// template<typename E1>
// auto Su
// (
//     const E1& expression,
//     const GeometricField<typename E1::value_type, fvPatchField, volMesh>& fld
// )
// {
//     // Wrap zero
//     const Expression::UniformListWrap<scalar> constantZero(fld.size(), 0.0);

//     // Wrap mesh volume
//     const ListConstRefWrap<scalar> V(fld.mesh().V());

//     return fvMatrixInternal2(fld, constantZero, -V*expression);
// }


// Some discretisation using expressions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<typename E1, typename E2>
class interpolate
:
    public GeometricFieldExpression
    <
        interpolate<E1, E2>,

        //- Type to return for internal field
        List_interpolate
        <
            IndirectConstWrap<typename E1::IntExpr>,
            IndirectConstWrap<typename E1::IntExpr>,
            typename E2::IntExpr
        >,

        //- Type to return for uncoupled patch field
        typename E1::UncoupledPatchExpr,

        //- Type to return for coupled patch field
        List_interpolate
        <
            IndirectConstWrap<typename E1::IntExpr>,
            typename E1::CoupledPatchExpr,
            typename E2::CoupledPatchExpr
        >,

        //- Type of actual values
        typename E1::value_type
    >
{
    //- Internal values
    typename std::conditional
    <
        E1::is_leaf,
        const E1&,
        const E1
    >::type cellVals_;

    //- Interpolation weights
    typename std::conditional
    <
        E2::is_leaf,
        const E2&,
        const E2
    >::type faceWeights_;

    //- Addressing container
    const fvMesh& mesh_;

public:

    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef List_interpolate
    <
        IndirectConstWrap<typename E1::IntExpr>,
        IndirectConstWrap<typename E1::IntExpr>,
        typename E2::IntExpr
    > IntExpr;

    //- Type to return for uncoupled patch field
    typedef typename E1::UncoupledPatchExpr UncoupledPatchExpr;

    //- Type to return for coupled patch field
    typedef List_interpolate
    <
        IndirectConstWrap<typename E1::IntExpr>,
        typename E1::CoupledPatchExpr,
        typename E2::CoupledPatchExpr
    > CoupledPatchExpr;

    //- Type of actual values
    typedef typename E1::value_type value_type;


    // Construct from components
    interpolate
    (
        const E1& cellVals,
        const E2& faceWeights,
        const fvMesh& mesh
    )
    :
        GeometricFieldExpression
        <
            interpolate<E1, E2>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            cellVals.dimensions()*faceWeights.dimensions(),
            faceWeights.oriented()
        ),
        cellVals_(cellVals),
        faceWeights_(faceWeights),
        mesh_(mesh)
    {}

    value_type operator[](const label i) const
    {
        const auto ownVal = cellVals_[mesh_.owner()[i]];
        const auto neiVal = cellVals_[mesh_.neighbour()[i]];
        return faceWeights_[i]*(ownVal-neiVal) + neiVal;
    }

    label size() const
    {
        return faceWeights_.size(); // number of internal faces
    }

    IntExpr internalField() const
    {
        return IntExpr
        (
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.owner()
            ),
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.neighbour()
            ),
            faceWeights_.internalField()
        );
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        return UncoupledPatchExpr(cellVals_.patchField(i));
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        return CoupledPatchExpr
        (
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.boundaryMesh()[i].faceCells()
            ),
            cellVals_.access
            (
                [&](const auto& fld, const label i)
                {
                    return ListConstTmpWrap<Field<value_type>>
                    (
                        fld.boundaryField()[i].patchNeighbourField()
                    );
                },
                i
            ),
            faceWeights_.coupledPatchField(i)
        );
    }
};


template<typename E1>
class linearInterpolate
:
    public interpolate
    <
        E1,
        GeometricFieldConstRefWrap<surfaceScalarField>
    >
{
public:

    // Construct from components
    linearInterpolate
    (
        const E1& cellVals,
        const fvMesh& mesh
    )
    :
        interpolate
        <
            E1,
            GeometricFieldConstRefWrap<surfaceScalarField>
        >
        (
            cellVals,
            mesh.surfaceInterpolation::weights().expr(),
            mesh
        )
    {}
};


template<class Type, class E1, class E2>
void fvmLaplacianUncorrected
(
    fvMatrix<Type>& fvm,
    const E1& gammaMagSf,
    const E2& deltaCoeffs
)
{
    typedef Expression::ListConstRefWrap<Type> expr;
    typedef Expression::ListRefWrap<Type> evaluator;

    const auto& vf = fvm.psi();

    // Set dimensions. Or should we just check?
    const_cast<dimensionSet&>(fvm.dimensions()).reset
    (
        deltaCoeffs.dimensions()
      * gammaMagSf.dimensions()
      * vf.dimensions()
    );

    Expression::ListRefWrap<scalar>
    (
        fvm.upper(),
        deltaCoeffs.internalField()*gammaMagSf.internalField()
    );
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const auto& pvf = vf.boundaryField()[patchi];
        auto intCoeffs = fvm.internalCoeffs(patchi);
        auto bouCoeffs = fvm.boundaryCoeffs(patchi);

        if (pvf.coupled())
        {
            // Evaluate deltaCoeffs
            scalarField pDeltaCoeffs;
            ListRefWrap<scalar>(pDeltaCoeffs, deltaCoeffs.patchField(patchi));

            //intCoeffs = pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            evaluator
            (
                intCoeffs,
                expr(pvf.gradientInternalCoeffs(pDeltaCoeffs)())
               *gammaMagSf.coupledPatchField(patchi)
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            evaluator
            (
                bouCoeffs,
               -expr(pvf.gradientBoundaryCoeffs(pDeltaCoeffs)())
               *gammaMagSf.coupledPatchField(patchi)
            );
        }
        else
        {
            //intCoeffs = pGamma*pvf.gradientInternalCoeffs();
            evaluator
            (
                intCoeffs,
                expr(pvf.gradientInternalCoeffs()())
               *gammaMagSf.patchField(patchi)
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs();
            evaluator
            (
                bouCoeffs,
               -expr(pvf.gradientBoundaryCoeffs()())
               *gammaMagSf.patchField(patchi)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
