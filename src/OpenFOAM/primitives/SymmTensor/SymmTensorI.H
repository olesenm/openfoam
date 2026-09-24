/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "Tensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor(Foam::zero)
:
    SymmTensor::vsType(Foam::zero{})
{}


template<class Cmpt>
template<class Cmpt2>
inline Foam::SymmTensor<Cmpt>::SymmTensor
(
    const VectorSpace<SymmTensor<Cmpt2>, Cmpt2, 6>& vs
)
:
    SymmTensor::vsType(vs)
{}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = Zero; this->v_[XZ] = Zero;
    this->v_[YY] = st.ii(); this->v_[YZ] = Zero;
    this->v_[ZZ] = st.ii();
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor
(
    const Vector<Cmpt>& x,
    const Vector<Cmpt>& y,
    const Vector<Cmpt>& z,
    const bool transposed  /* ignored */
)
{
    this->rows(x, y, z);
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor
(
    const Cmpt txx, const Cmpt txy, const Cmpt txz,
                    const Cmpt tyy, const Cmpt tyz,
                                    const Cmpt tzz
)
{
    this->v_[XX] = txx; this->v_[XY] = txy; this->v_[XZ] = txz;
                        this->v_[YY] = tyy; this->v_[YZ] = tyz;
                                            this->v_[ZZ] = tzz;
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor(Istream& is)
:
    SymmTensor::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SymmTensor<Cmpt>::x() const
{
    return Vector<Cmpt>(this->v_[XX], this->v_[XY], this->v_[XZ]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SymmTensor<Cmpt>::y() const
{
    return Vector<Cmpt>(this->v_[XY], this->v_[YY], this->v_[YZ]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SymmTensor<Cmpt>::z() const
{
    return Vector<Cmpt>(this->v_[XZ], this->v_[YZ], this->v_[ZZ]);
}


template<class Cmpt>
template<Foam::direction Idx>
inline Foam::Vector<Cmpt> Foam::SymmTensor<Cmpt>::row() const
{
    if (Idx == 0) return x();
    else if (Idx == 1) return y();
    else if (Idx == 2) return z();

    static_assert(Idx < 3, "Invalid row access");
    return Zero;
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SymmTensor<Cmpt>::row(const direction r) const
{
    switch (r)
    {
        case 0: return x(); break;
        case 1: return y(); break;
        case 2: return z(); break;
        default:
            FatalErrorInFunction
                << "Invalid row access " << r << abort(FatalError);
    }

    return Zero;
}


template<class Cmpt>
template<Foam::direction Idx>
inline void Foam::SymmTensor<Cmpt>::row(const Vector<Cmpt>& v)
{
    if (Idx == 0)
    {
        this->v_[XX] = v.x(); this->v_[XY] = v.y(); this->v_[XZ] = v.z();
    }
    else if (Idx == 1)
    {
        this->v_[XY] = v.x(); this->v_[YY] = v.y(); this->v_[YZ] = v.z();
    }
    else if (Idx == 2)
    {
        this->v_[XZ] = v.x(); this->v_[YZ] = v.y(); this->v_[ZZ] = v.z();
    }

    static_assert(Idx < 3, "Invalid row access");
}


template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::rows
(
    const Vector<Cmpt>& x,
    const Vector<Cmpt>& y,
    const Vector<Cmpt>& z
)
{
    this->v_[XX] = x.x(); this->v_[XY] = x.y(); this->v_[XZ] = x.z();
    this->v_[YY] = y.y(); this->v_[YZ] = y.z();
    this->v_[ZZ] = z.z();
}


template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::row
(
    const direction r,
    const Vector<Cmpt>& v
)
{
    switch (r)
    {
        case 0: row<0>(v); break;
        case 1: row<1>(v); break;
        case 2: row<2>(v); break;
        default:
            FatalErrorInFunction
                << "Invalid row access " << r << abort(FatalError);
    }
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SymmTensor<Cmpt>::diag() const
{
    return Vector<Cmpt>(this->v_[XX], this->v_[YY], this->v_[ZZ]);
}


template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::diag(const Vector<Cmpt>& v)
{
    this->v_[XX] = v.x(); this->v_[YY] = v.y(); this->v_[ZZ] = v.z();
}


template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::addDiag(const Vector<Cmpt>& v)
{
    this->v_[XX] += v.x(); this->v_[YY] += v.y(); this->v_[ZZ] += v.z();
}


template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::subtractDiag(const Vector<Cmpt>& v)
{
    this->v_[XX] -= v.x(); this->v_[YY] -= v.y(); this->v_[ZZ] -= v.z();
}


template<class Cmpt>
inline Foam::scalar Foam::SymmTensor<Cmpt>::diagSqr() const
{
    return
    (
        Foam::magSqr(this->xx())
      + Foam::magSqr(this->yy())
      + Foam::magSqr(this->zz())
    );
}


template<class Cmpt>
inline bool Foam::SymmTensor<Cmpt>::is_identity(const scalar tol) const
{
    return
    (
        Foam::mag(xx() - pTraits<Cmpt>::one) < tol
     && Foam::mag(yy() - pTraits<Cmpt>::one) < tol
     && Foam::mag(zz() - pTraits<Cmpt>::one) < tol
     && Foam::mag(xy()) < tol
     && Foam::mag(xz()) < tol
     && Foam::mag(yz()) < tol
    );
}


template<class Cmpt>
inline Cmpt Foam::SymmTensor<Cmpt>::det() const
{
    return
    (
        xx()*yy()*zz() + xy()*yz()*xz()
      + xz()*xy()*yz() - xx()*yz()*yz()
      - xy()*xy()*zz() - xz()*yy()*xz()
    );
}


template<class Cmpt>
inline Cmpt Foam::SymmTensor<Cmpt>::det2D(const direction excludeCmpt) const
{
    switch (excludeCmpt)
    {
        case 0:  // Eliminate x
        {
            return (yy()*zz() - yz()*yz());
        }

        case 1:  // Eliminate y
        {
            return (xx()*zz() - xz()*xz());
        }
    }

    // Fall-through: Eliminate z
    return (xx()*yy() - xy()*xy());
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt> Foam::SymmTensor<Cmpt>::adjunct() const
{
    // symmetric: cof() == adjunct()
    return SymmTensor<Cmpt>
    (
        yy()*zz() - yz()*yz(), xz()*yz() - xy()*zz(), xy()*yz() - xz()*yy(),
                               xx()*zz() - xz()*xz(), xy()*xz() - xx()*yz(),
                                                      xx()*yy() - xy()*xy()
    );
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt> Foam::SymmTensor<Cmpt>::cof() const
{
    // symmetric: cof() == adjunct()
    return this->adjunct();
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt> Foam::SymmTensor<Cmpt>::adjunct2D
(
    const direction excludeCmpt
) const
{
    switch (excludeCmpt)
    {
        case 0:  // Eliminate x
        {
            return SymmTensor<Cmpt>
            (
                Zero, Zero,  Zero,
                      zz(), -yz(),
                             yy()
            );
        }

        case 1:  // Eliminate y
        {
            return SymmTensor<Cmpt>
            (
                zz(), Zero, -xz(),
                      Zero,  Zero,
                             xx()
            );
        }
    }

    // Fall-through: Eliminate z
    return SymmTensor<Cmpt>
    (
        yy(), -xy(), Zero,
               xx(), Zero,
                     Zero
    );
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt> Foam::SymmTensor<Cmpt>::inv2D
(
    const direction excludeCmpt
) const
{
    const Cmpt detval = this->det2D(excludeCmpt);

    return this->adjunct2D(excludeCmpt)/detval;
}


// Invert without much error handling
template<class Cmpt>
inline Foam::SymmTensor<Cmpt> Foam::SymmTensor<Cmpt>::inv() const
{
    const Cmpt detval = this->det();

    #ifdef FULLDEBUG
    if (mag(detval) < VSMALL)
    {
        FatalErrorInFunction
            << "Tensor not properly invertible, determinant:"
            << detval << " tensor:" << *this << nl
            << abort(FatalError);
    }
    #endif

    return this->adjunct()/detval;
}


// Invert with some error handling
template<class Cmpt>
inline Foam::SymmTensor<Cmpt> Foam::SymmTensor<Cmpt>::safeInv() const
{
    {
        // Attempt to identify and handle 2-D cases.
        // - use diagSqr instead of magSqr for fewer operations
        const scalar magSqr_xx = Foam::magSqr(xx());
        const scalar magSqr_yy = Foam::magSqr(yy());
        const scalar magSqr_zz = Foam::magSqr(zz());

        // SMALL: 1e-15 (double), 1e-6 (float), but 1e-6 may be adequate

        const scalar threshold = SMALL * (magSqr_xx + magSqr_yy + magSqr_zz);

        const bool small_xx = (magSqr_xx < threshold);
        const bool small_yy = (magSqr_yy < threshold);
        const bool small_zz = (magSqr_zz < threshold);

        if (small_xx || small_yy || small_zz)
        {
            SymmTensor<Cmpt> work(*this);

            if (small_xx) { work.xx() += 1.0; }
            if (small_yy) { work.yy() += 1.0; }
            if (small_zz) { work.zz() += 1.0; }

            const Cmpt detval = work.det();

            if (mag(detval) < ROOTVSMALL)
            {
                // Appears to be nearly zero - leave untouched?
                return SymmTensor<Cmpt>(Zero);
            }

            work = work.adjunct()/detval;

            if (small_xx) { work.xx() -= 1.0; }
            if (small_yy) { work.yy() -= 1.0; }
            if (small_zz) { work.zz() -= 1.0; }

            return work;
        }
    }

    const Cmpt detval = this->det();

    if (mag(detval) < ROOTVSMALL)
    {
        // Appears to be nearly zero - leave untouched?
        return SymmTensor<Cmpt>(Zero);
    }

    return this->adjunct()/detval;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::operator=(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = Zero;    this->v_[XZ] = Zero;
                            this->v_[YY] = st.ii(); this->v_[YZ] = Zero;
                                                    this->v_[ZZ] = st.ii();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the trace of a SymmTensor
template<class Cmpt>
inline Cmpt tr(const SymmTensor<Cmpt>& st)
{
    return st.xx() + st.yy() + st.zz();
}


//- Return the spherical part of a SymmTensor
template<class Cmpt>
inline SphericalTensor<Cmpt> sph(const SymmTensor<Cmpt>& st)
{
    return SphericalTensor<Cmpt>
    (
        (1.0/3.0)*tr(st)
    );
}


//- Return the symmetric part of a SymmTensor, i.e. itself
template<class Cmpt>
inline const SymmTensor<Cmpt>& symm(const SymmTensor<Cmpt>& st)
{
    return st;
}


//- Return twice the symmetric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> twoSymm(const SymmTensor<Cmpt>& st)
{
    return 2*st;
}


//- Return the deviatoric part of the symmetric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> devSymm(const SymmTensor<Cmpt>& st)
{
    return dev(st);
}


//- Return the deviatoric part of twice the symmetric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> devTwoSymm(const SymmTensor<Cmpt>& st)
{
    return dev2(st);
}


//- Return the deviatoric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> dev(const SymmTensor<Cmpt>& st)
{
    return st - sph(st);
}


//- Return the two-third deviatoric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> dev2(const SymmTensor<Cmpt>& st)
{
    return st - 2*sph(st);
}


//- Return the determinant of a SymmTensor
template<class Cmpt>
inline Cmpt det(const SymmTensor<Cmpt>& st)
{
    return st.det();
}


//- Return the cofactor of a SymmTensor as a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> cof(const SymmTensor<Cmpt>& st)
{
    return st.cof();
}


//- Return the inverse of a SymmTensor, using given determinant value
template<class Cmpt>
inline SymmTensor<Cmpt> inv(const SymmTensor<Cmpt>& st, const Cmpt detval)
{
    #ifdef FULLDEBUG
    if (mag(detval) < VSMALL)
    {
        FatalErrorInFunction
            << "SymmTensor not properly invertible, determinant:"
            << detval << " tensor:" << st << nl
            << abort(FatalError);
    }
    #endif

    return st.adjunct()/detval;
}


//- Return the inverse of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> inv(const SymmTensor<Cmpt>& st)
{
    return st.inv();
}


//- Return the 1st invariant of a SymmTensor
template<class Cmpt>
inline Cmpt invariantI(const SymmTensor<Cmpt>& st)
{
    return tr(st);
}


//- Return the 2nd invariant of a SymmTensor
template<class Cmpt>
inline Cmpt invariantII(const SymmTensor<Cmpt>& st)
{
    return
    (
        st.xx()*st.yy() + st.yy()*st.zz() + st.xx()*st.zz()
      - st.xy()*st.xy() - st.yz()*st.yz() - st.xz()*st.xz()
    );
}


//- Return the 3rd invariant of a SymmTensor
template<class Cmpt>
inline Cmpt invariantIII(const SymmTensor<Cmpt>& st)
{
    return det(st);
}


//- Return the inner-product of a SymmTensor with itself
template<class Cmpt>
inline SymmTensor<Cmpt>
innerSqr(const SymmTensor<Cmpt>& st)
{
    return SymmTensor<Cmpt>
    (
        st.xx()*st.xx() + st.xy()*st.xy() + st.xz()*st.xz(),
        st.xx()*st.xy() + st.xy()*st.yy() + st.xz()*st.yz(),
        st.xx()*st.xz() + st.xy()*st.yz() + st.xz()*st.zz(),

        st.xy()*st.xy() + st.yy()*st.yy() + st.yz()*st.yz(),
        st.xy()*st.xz() + st.yy()*st.yz() + st.yz()*st.zz(),

        st.xz()*st.xz() + st.yz()*st.yz() + st.zz()*st.zz()
    );
}


//- Return the square of Frobenius norm of a SymmTensor
template<class Cmpt>
inline Foam::scalar magSqr(const SymmTensor<Cmpt>& st)
{
    return
    (
        magSqr(st.xx()) + 2*magSqr(st.xy()) + 2*magSqr(st.xz())
                        +   magSqr(st.yy()) + 2*magSqr(st.yz())
                                            +   magSqr(st.zz())
    );
}


//- Return the square of a Vector as a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> sqr(const Vector<Cmpt>& v)
{
    return SymmTensor<Cmpt>
    (
        v.x()*v.x(), v.x()*v.y(), v.x()*v.z(),
                     v.y()*v.y(), v.y()*v.z(),
                                  v.z()*v.z()
    );
}


//- Linear interpolation of symmetric tensors a and b by factor t
template<class Cmpt>
inline SymmTensor<Cmpt> lerp
(
    const SymmTensor<Cmpt>& a,
    const SymmTensor<Cmpt>& b,
    const scalar t
)
{
    const scalar onet = (1-t);

    return SymmTensor<Cmpt>
    (
        onet*a.xx() + t*b.xx(),
        onet*a.xy() + t*b.xy(),
        onet*a.xz() + t*b.xz(),
        onet*a.yy() + t*b.yy(),
        onet*a.yz() + t*b.yz(),
        onet*a.zz() + t*b.zz()
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Sum of a SphericalTensor and a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator+(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return SymmTensor<Cmpt>
    (
        spt1.ii() + st2.xx(), st2.xy(),             st2.xz(),
                              spt1.ii() + st2.yy(), st2.yz(),
                                                    spt1.ii() + st2.zz()
    );
}


//- Sum of a SymmTensor and a SphericalTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator+(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return SymmTensor<Cmpt>
    (
        st1.xx() + spt2.ii(), st1.xy(),             st1.xz(),
                              st1.yy() + spt2.ii(), st1.yz(),
                                                    st1.zz() + spt2.ii()
    );
}


//- Subtract a SymmTensor from a SphericalTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator-(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return SymmTensor<Cmpt>
    (
        spt1.ii() - st2.xx(), -st2.xy(),             -st2.xz(),
                               spt1.ii() - st2.yy(), -st2.yz(),
                                                      spt1.ii() - st2.zz()
    );
}


//- Subtract a SphericalTensor from a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator-(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return SymmTensor<Cmpt>
    (
        st1.xx() - spt2.ii(), st1.xy(),             st1.xz(),
                              st1.yy() - spt2.ii(), st1.yz(),
                                                    st1.zz() - spt2.ii()
    );
}


//- Return the Hodge dual of a SymmTensor as a Vector
template<class Cmpt>
inline Vector<Cmpt> operator*(const SymmTensor<Cmpt>& st)
{
    return Vector<Cmpt>(st.yz(), -st.xz(), st.xy());
}


//- Division of a SymmTensor by a Cmpt
template<class Cmpt>
inline SymmTensor<Cmpt>
operator/(const SymmTensor<Cmpt>& st, const Cmpt s)
{
    return SymmTensor<Cmpt>
    (
        st.xx()/s, st.xy()/s, st.xz()/s,
                   st.yy()/s, st.yz()/s,
                              st.zz()/s
    );
}


//- Inner-product of a SymmTensor and a SymmTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const SymmTensor<Cmpt>& st1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        st1.xx()*st2.xx() + st1.xy()*st2.xy() + st1.xz()*st2.xz(),
        st1.xx()*st2.xy() + st1.xy()*st2.yy() + st1.xz()*st2.yz(),
        st1.xx()*st2.xz() + st1.xy()*st2.yz() + st1.xz()*st2.zz(),

        st1.xy()*st2.xx() + st1.yy()*st2.xy() + st1.yz()*st2.xz(),
        st1.xy()*st2.xy() + st1.yy()*st2.yy() + st1.yz()*st2.yz(),
        st1.xy()*st2.xz() + st1.yy()*st2.yz() + st1.yz()*st2.zz(),

        st1.xz()*st2.xx() + st1.yz()*st2.xy() + st1.zz()*st2.xz(),
        st1.xz()*st2.xy() + st1.yz()*st2.yy() + st1.zz()*st2.yz(),
        st1.xz()*st2.xz() + st1.yz()*st2.yz() + st1.zz()*st2.zz()
    );
}


//- Inner-product of a SphericalTensor and a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return SymmTensor<Cmpt>
    (
        spt1.ii()*st2.xx(), spt1.ii()*st2.xy(), spt1.ii()*st2.xz(),
                            spt1.ii()*st2.yy(), spt1.ii()*st2.yz(),
                                                spt1.ii()*st2.zz()
    );
}


//- Inner-product of a SymmTensor and a SphericalTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator&(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return SymmTensor<Cmpt>
    (
        st1.xx()*spt2.ii(), st1.xy()*spt2.ii(), st1.xz()*spt2.ii(),
                            st1.yy()*spt2.ii(), st1.yz()*spt2.ii(),
                                                st1.zz()*spt2.ii()
    );
}


//- Inner-product of a SymmTensor and a Vector
template<class Cmpt>
inline Vector<Cmpt>
operator&(const SymmTensor<Cmpt>& st, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        st.xx()*v.x() + st.xy()*v.y() + st.xz()*v.z(),
        st.xy()*v.x() + st.yy()*v.y() + st.yz()*v.z(),
        st.xz()*v.x() + st.yz()*v.y() + st.zz()*v.z()
    );
}


//- Inner-product of a Vector and a SymmTensor
template<class Cmpt>
inline Vector<Cmpt>
operator&(const Vector<Cmpt>& v, const SymmTensor<Cmpt>& st)
{
    return Vector<Cmpt>
    (
        v.x()*st.xx() + v.y()*st.xy() + v.z()*st.xz(),
        v.x()*st.xy() + v.y()*st.yy() + v.z()*st.yz(),
        v.x()*st.xz() + v.y()*st.yz() + v.z()*st.zz()
    );
}


//- Double-inner-product of a SymmTensor and a SymmTensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor<Cmpt>& st1, const SymmTensor<Cmpt>& st2)
{
    return
    (
        st1.xx()*st2.xx() + 2*st1.xy()*st2.xy() + 2*st1.xz()*st2.xz()
                          +   st1.yy()*st2.yy() + 2*st1.yz()*st2.yz()
                                                +   st1.zz()*st2.zz()
    );
}


//- Double-inner-product of a SphericalTensor and a SymmTensor
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return (spt1.ii()*st2.xx() + spt1.ii()*st2.yy() + spt1.ii()*st2.zz());
}


//- Double-inner-product of a SymmTensor and a SphericalTensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return (st1.xx()*spt2.ii() + st1.yy()*spt2.ii() + st1.zz()*spt2.ii());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class outerProduct<SymmTensor<Cmpt>, Cmpt>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, Vector<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<SphericalTensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class typeOfSum<SymmTensor<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

#ifdef FOAM_OFFLOAD

// GPU/stdpar-safe: no SVD
inline SymmTensor<scalar> pinv(const SymmTensor<scalar>& st)
{
    const scalar detval = st.det();

    // If invertible enough -> exact inverse
    if (mag(detval) >= ROOTVSMALL)
    {
        return inv(st, detval);
    }

    // Near-singular fallback: diagonal regularization (Tikhonov-like)
    // This is not Moore–Penrose, but is GPU-safe and robust.
    const scalar eps = ROOTVSMALL;  // you can tune to SMALL if needed

    SymmTensor<scalar> reg(st);
    reg.xx() += eps;
    reg.yy() += eps;
    reg.zz() += eps;

    const scalar det2 = reg.det();
    if (mag(det2) < ROOTVSMALL)
    {
        return SymmTensor<scalar>(Zero);
    }

    return inv(reg, det2);
}

inline void sort3Ascending(scalar& a, scalar& b, scalar& c)
{
    if (a > b) { scalar t=a; a=b; b=t; }
    if (b > c) { scalar t=b; b=c; c=t; }
    if (a > b) { scalar t=a; a=b; b=t; }
}

inline Vector<scalar> eigenValues(const SymmTensor<scalar>& T)
{
    const scalar a11 = T.xx();
    const scalar a22 = T.yy();
    const scalar a33 = T.zz();
    const scalar a12 = T.xy();
    const scalar a13 = T.xz();
    const scalar a23 = T.yz();

    const scalar off2 = sqr(a12) + sqr(a13) + sqr(a23);
    if (off2 < ROOTSMALL)
    {
        scalar l1=a11, l2=a22, l3=a33;
        sort3Ascending(l1, l2, l3);
        return Vector<scalar>(l1, l2, l3);
    }

    // Analytic eigensolver for symmetric 3x3
    const scalar m = (a11 + a22 + a33) / scalar(3);

    const scalar b11 = a11 - m;
    const scalar b22 = a22 - m;
    const scalar b33 = a33 - m;

    const scalar p2 =
        (b11*b11 + b22*b22 + b33*b33
      + scalar(2)*(a12*a12 + a13*a13 + a23*a23)) / scalar(6);

    if (p2 < VSMALL)
    {
        return Vector<scalar>(m, m, m);
    }

    const scalar p = Foam::sqrt(p2);

    const scalar detB =
          b11*(b22*b33 - a23*a23)
        - a12*(a12*b33 - a13*a23)
        + a13*(a12*a23 - a13*b22);

    scalar r = detB / (scalar(2) * p*p*p);
    if (r < scalar(-1)) r = scalar(-1);
    if (r > scalar( 1)) r = scalar( 1);

    const scalar phi = Foam::acos(r) / scalar(3);

    const scalar twoPiOver3  = scalar(2.0943951023931954923);
    const scalar fourPiOver3 = scalar(4.1887902047863909846);

    scalar l1 = m + scalar(2)*p*Foam::cos(phi);
    scalar l2 = m + scalar(2)*p*Foam::cos(phi + twoPiOver3);
    scalar l3 = m + scalar(2)*p*Foam::cos(phi + fourPiOver3);

    sort3Ascending(l1, l2, l3);
    return Vector<scalar>(l1, l2, l3);
}

inline Vector<scalar> safeNormalize(const Vector<scalar>& v, const Vector<scalar>& fallback)
{
    const scalar mv = mag(v);
    return (mv > SMALL) ? (v/mv) : fallback;
}

inline Vector<scalar> eigenVector
(
    const SymmTensor<scalar>& T,
    const scalar eVal,
    const Vector<scalar>& standardBasis1,
    const Vector<scalar>& standardBasis2
)
{
    // A = T - eVal*I
    const scalar Axx = T.xx() - eVal;
    const scalar Ayy = T.yy() - eVal;
    const scalar Azz = T.zz() - eVal;

    const scalar Axy = T.xy();
    const scalar Axz = T.xz();
    const scalar Ayz = T.yz();

    const scalar sd0 = Ayy*Azz - Ayz*Ayz;
    const scalar sd1 = Azz*Axx - Axz*Axz;
    const scalar sd2 = Axx*Ayy - Axy*Axy;

    const scalar magSd0 = mag(sd0);
    const scalar magSd1 = mag(sd1);
    const scalar magSd2 = mag(sd2);

    const Vector<scalar> fallback =
        safeNormalize(standardBasis1 ^ standardBasis2, standardBasis1);

    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        const Vector<scalar> eVec
        (
            1,
            (Ayz*Axz - Azz*Axy)/sd0,
            (Ayz*Axy - Ayy*Axz)/sd0
        );
        return safeNormalize(eVec, fallback);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        const Vector<scalar> eVec
        (
            (Axz*Ayz - Azz*Axy)/sd1,
            1,
            (Axz*Axy - Axx*Ayz)/sd1
        );
        return safeNormalize(eVec, fallback);
    }
    else if (magSd2 > SMALL)
    {
        const Vector<scalar> eVec
        (
            (Axy*Ayz - Ayy*Axz)/sd2,
            (Axy*Axz - Axx*Ayz)/sd2,
            1
        );
        return safeNormalize(eVec, fallback);
    }

    // repeated eigenvalue fallback
    return fallback;
}

inline Tensor<scalar> eigenVectors(const SymmTensor<scalar>& T, const Vector<scalar>& eVals)
{
    Vector<scalar> Ux(1,0,0), Uy(0,1,0), Uz(0,0,1);

    Ux = eigenVector(T, eVals.x(), Uy, Uz);
    Uy = eigenVector(T, eVals.y(), Uz, Ux);

    // orthogonalize Uy against Ux for stability
    Uy -= (Uy & Ux)*Ux;
    Uy = safeNormalize(Uy, safeNormalize(Ux ^ Uz, Vector<scalar>(0,1,0)));

    Uz = Ux ^ Uy;
    Uz = safeNormalize(Uz, Vector<scalar>(0,0,1));

    return Tensor<scalar>(Ux, Uy, Uz);
}

inline Tensor<scalar> eigenVectors(const SymmTensor<scalar>& T)
{
    const Vector<scalar> eVals(eigenValues(T));
    return eigenVectors(T, eVals);
}

#endif // FOAM_OFFLOAD

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
