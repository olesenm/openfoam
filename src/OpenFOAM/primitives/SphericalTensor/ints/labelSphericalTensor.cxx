/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "labelSphericalTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::SphericalTensor<Foam::label>::vsType::typeName
(
    "labelSphericalTensor"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  defineTraits
#define defineTraits(Type, Prefix)                                            \
                                                                              \
template<>                                                                    \
const char* const                                                             \
Foam::SphericalTensor<Type>::vsType::componentNames[] = { "ii" };             \
                                                                              \
template<>                                                                    \
const Foam::SphericalTensor<Type>                                             \
Foam::SphericalTensor<Type>::vsType::zero                                     \
(                                                                             \
    SphericalTensor<Type>::uniform(0)                                         \
);                                                                            \
                                                                              \
template<>                                                                    \
const Foam::SphericalTensor<Type>                                             \
Foam::SphericalTensor<Type>::vsType::one                                      \
(                                                                             \
    SphericalTensor<Type>::uniform(1)                                         \
);                                                                            \
                                                                              \
template<>                                                                    \
const Foam::SphericalTensor<Type>                                             \
Foam::SphericalTensor<Type>::vsType::max                                      \
(                                                                             \
    SphericalTensor<Type>::uniform(Prefix##Max)                               \
);                                                                            \
                                                                              \
template<>                                                                    \
const Foam::SphericalTensor<Type>                                             \
Foam::SphericalTensor<Type>::vsType::min                                      \
(                                                                             \
    SphericalTensor<Type>::uniform(-Prefix##Max)                              \
);                                                                            \
                                                                              \
template<>                                                                    \
const Foam::SphericalTensor<Type>                                             \
Foam::SphericalTensor<Type>::vsType::rootMax                                  \
(                                                                             \
    SphericalTensor<Type>::uniform(std::sqrt(double(Prefix##Max)))            \
);                                                                            \
                                                                              \
template<>                                                                    \
const Foam::SphericalTensor<Type>                                             \
Foam::SphericalTensor<Type>::vsType::rootMin                                  \
(                                                                             \
    SphericalTensor<Type>::uniform(-std::sqrt(double(Prefix##Max)))           \
);                                                                            \
                                                                              \
template<>                                                                    \
const Foam::SphericalTensor<Type>                                             \
Foam::SphericalTensor<Type>::I(1);


defineTraits(Foam::label, label);

#undef defineTraits


// ************************************************************************* //
