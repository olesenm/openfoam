/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "patchFieldProbe.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchFieldProbe::sample(const VolumeField<Type>& vField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    const auto& pbm = thisMesh_.boundaryMesh();
    const auto& bField = vField.boundaryField();

    forAll(*this, probei)
    {
        if (label facei = faceIds_[probei]; facei >= 0)
        {
            auto patchFacePair = pbm.whichPatchFace(facei);
            values[probei] = bField[patchFacePair];
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchFieldProbe::sample(const SurfaceField<Type>& sField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    const auto& pbm = thisMesh_.boundaryMesh();
    const auto& bField = sField.boundaryField();

    forAll(*this, probei)
    {
        if (label facei = faceIds_[probei]; facei >= 0)
        {
            auto patchFacePair = pbm.whichPatchFace(facei);
            values[probei] = bField[patchFacePair];
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchFieldProbe::sample(const word& fieldName) const
{
    return sample(thisMesh_.lookupObject<VolumeField<Type>>(fieldName));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchFieldProbe::sampleSurfaceField(const word& fieldName) const
{
    return sample(thisMesh_.lookupObject<SurfaceField<Type>>(fieldName));
}


// ************************************************************************* //
