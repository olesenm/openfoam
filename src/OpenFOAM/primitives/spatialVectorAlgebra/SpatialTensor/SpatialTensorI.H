/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "Identity.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline constexpr Foam::SpatialTensor<Cmpt>::SpatialTensor(Foam::zero)
:
    SpatialTensor::msType(Foam::zero{})
{}


template<class Cmpt>
inline Foam::SpatialTensor<Cmpt>::SpatialTensor
(
    const typename SpatialTensor::msType& ms
)
:
    SpatialTensor::msType(ms)
{}


template<class Cmpt>
inline Foam::SpatialTensor<Cmpt>::SpatialTensor
(
    const Tensor<Cmpt>& t00, const Tensor<Cmpt>& t01,
    const Tensor<Cmpt>& t10, const Tensor<Cmpt>& t11
)
{
    // Block (0, 0)
    this->v_[0] = t00.xx();   this->v_[1] = t00.xy();   this->v_[2] = t00.xz();
    this->v_[6] = t00.yx();   this->v_[7] = t00.yy();   this->v_[8] = t00.yz();
    this->v_[12] = t00.zx();  this->v_[13] = t00.zy();  this->v_[14] = t00.zz();

    // Block (0, 1)
    this->v_[3] = t01.xx();   this->v_[4] = t01.xy();   this->v_[5] = t01.xz();
    this->v_[9] = t01.yx();   this->v_[10] = t01.yy();  this->v_[11] = t01.yz();
    this->v_[15] = t01.zx();  this->v_[16] = t01.zy();  this->v_[17] = t01.zz();

    // Block (1, 0)
    this->v_[18] = t10.xx();  this->v_[19] = t10.xy();  this->v_[20] = t10.xz();
    this->v_[24] = t10.yx();  this->v_[25] = t10.yy();  this->v_[26] = t10.yz();
    this->v_[30] = t10.zx();  this->v_[31] = t10.zy();  this->v_[32] = t10.zz();

    // Block (1, 1)
    this->v_[21] = t11.xx();  this->v_[22] = t11.xy();  this->v_[23] = t11.xz();
    this->v_[27] = t11.yx();  this->v_[28] = t11.yy();  this->v_[29] = t11.yz();
    this->v_[33] = t11.zx();  this->v_[34] = t11.zy();  this->v_[35] = t11.zz();
}


template<class Cmpt>
inline Foam::SpatialTensor<Cmpt>::SpatialTensor
(
    const Cmpt& t00, const Cmpt& t01, const Cmpt& t02,
    const Cmpt& t03, const Cmpt& t04, const Cmpt& t05,

    const Cmpt& t10, const Cmpt& t11, const Cmpt& t12,
    const Cmpt& t13, const Cmpt& t14, const Cmpt& t15,

    const Cmpt& t20, const Cmpt& t21, const Cmpt& t22,
    const Cmpt& t23, const Cmpt& t24, const Cmpt& t25,

    const Cmpt& t30, const Cmpt& t31, const Cmpt& t32,
    const Cmpt& t33, const Cmpt& t34, const Cmpt& t35,

    const Cmpt& t40, const Cmpt& t41, const Cmpt& t42,
    const Cmpt& t43, const Cmpt& t44, const Cmpt& t45,

    const Cmpt& t50, const Cmpt& t51, const Cmpt& t52,
    const Cmpt& t53, const Cmpt& t54, const Cmpt& t55
)
:
    // aggregate initialization
    SpatialTensor::msType
    (
        Foam::direct_construct_t{},
        t00, t01, t02,
        t03, t04, t05,
        t10, t11, t12,
        t13, t14, t15,
        t20, t21, t22,
        t23, t24, t25,
        t30, t31, t32,
        t33, t34, t35,
        t40, t41, t42,
        t43, t44, t45,
        t50, t51, t52,
        t53, t54, t55
    )
{}


template<class Cmpt>
inline Foam::SpatialTensor<Cmpt>::SpatialTensor(Istream& is)
:
    SpatialTensor::msType(is)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Return the cross-product tensor
template<class Cmpt>
inline Foam::SpatialTensor<Cmpt> operator^
(
    const SpatialVector<Cmpt>& v,
    const Identity<Cmpt>&
)
{
    return SpatialTensor<Cmpt>
    (
        0,       -v.wz(),   v.wy(),   0,        0,        0,
        v.wz(),   0,       -v.wx(),   0,        0,        0,
       -v.wy(),   v.wx(),   0,        0,        0,        0,
        0,       -v.lz(),   v.ly(),   0,       -v.wz(),   v.wy(),
        v.lz(),   0,       -v.lx(),   v.wz(),   0,       -v.wx(),
       -v.ly(),   v.lx(),   0,       -v.wy(),   v.wx(),   0
    );
}


//- Return the dual cross-product tensor
template<class Cmpt>
inline Foam::SpatialTensor<Cmpt> operator^
(
    const SpatialVector<Cmpt>& f,
    const typename Identity<Cmpt>::dual&
)
{
    return SpatialTensor<Cmpt>
    (
        0,       -f.wz(),   f.wy(),   0,       -f.lz(),   f.ly(),
        f.wz(),   0,       -f.wx(),   f.lz(),   0,       -f.lx(),
       -f.wy(),   f.wx(),   0,       -f.ly(),   f.lx(),   0,
        0,        0,        0,        0,       -f.wz(),   f.wy(),
        0,        0,        0,        f.wz(),   0,       -f.wx(),
        0,        0,        0,       -f.wy(),   f.wx(),   0
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
