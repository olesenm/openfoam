/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "fieldExtents.H"
#include "volFields.H"
#include "treeBoundBox.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldExtents::calcFieldExtents
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const auto* fieldPtr = obr_.cfindObject<VolFieldType>(fieldName);

    if (!fieldPtr)
    {
        return false;
    }
    const auto& field = *fieldPtr;


    // If the field value (or mag) is above the threshold, add its
    // corresponding position into the *local* bounding box.
    // The reference point offset is applied later.

    const scalar threshold = threshold_;

    auto filtered_bounds =
        [threshold](const UList<Type>& input, const auto& centres)
        {
            boundBox bb;

            const label len = std::min(input.size(), centres.size());
            for (label i = 0; i < len; ++i)
            {
                if constexpr (std::is_same_v<scalar, Type>)
                {
                    if (input[i] > threshold)
                    {
                        bb.add(centres[i]);
                    }
                }
                else
                {
                    if (Foam::mag(input[i]) > threshold)
                    {
                        bb.add(centres[i]);
                    }
                }
            }

            return bb;
        };

    Log << "field: " << fieldName << nl;

    writeCurrentTime(file());


    // Note: calculate bounding boxes as separate min/max limits
    List<point> bb_min(patchIDs_.size()+1);
    List<point> bb_max(patchIDs_.size()+1);

    // Internal field. Its boundBox is stored after the patch indices
    if (auto idx = patchIDs_.size(); internalField_)
    {
        const auto& input = field.primitiveField();
        const auto& centres = mesh_.C().primitiveField();

        const boundBox bb = filtered_bounds(input, centres);
        bb_min[idx] = bb.min();
        bb_max[idx] = bb.max();
    }
    else
    {
        bb_min[idx] = point::zero;
        bb_max[idx] = point::zero;
    }

    // Patches
    {
        const auto& fieldBoundary = field.boundaryField();
        const auto& CfBoundary = mesh_.C().boundaryField();

        forAll(patchIDs_, idx)
        {
            const auto patchi = patchIDs_[idx];
            const Field<Type>& input = fieldBoundary[patchi];
            const auto& centres = CfBoundary[patchi];

            const boundBox bb = filtered_bounds(input, centres);
            bb_min[idx] = bb.min();
            bb_max[idx] = bb.max();
        }
    }


    // Reduce all boundBox min/max values
    UPstream::mpiAllReduce_min(bb_min.data(), bb_min.size());
    UPstream::mpiAllReduce_max(bb_max.data(), bb_max.size());

    // Internal field
    if (auto idx = patchIDs_.size(); internalField_)
    {
        boundBox bb(bb_min[idx], bb_max[idx]);
        if (bb.empty())
        {
            bb.reset(point::zero);
        }
        else
        {
            // Apply reference position offset
            bb.min() -= this->C0_;
            bb.max() -= this->C0_;
        }

        Log << "    internal field: " << bb << nl;
        file() << bb;
        this->setResult(fieldName + "_internal_min", bb.min());
        this->setResult(fieldName + "_internal_max", bb.max());
    }

    // Patches
    forAll(patchIDs_, idx)
    {
        const label patchi = patchIDs_[idx];
        const word& patchName = mesh_.boundaryMesh()[patchi].name();

        boundBox bb(bb_min[idx], bb_max[idx]);
        if (bb.empty())
        {
            bb.reset(point::zero);
        }
        else
        {
            // Apply reference position offset
            bb.min() -= this->C0_;
            bb.max() -= this->C0_;
        }

        Log << "    patch " << patchName << ": " << bb << nl;
        file() << bb;
        this->setResult(fieldName + "_" + patchName + "_min", bb.min());
        this->setResult(fieldName + "_" + patchName + "_max", bb.max());
    }

    Log << endl;
    file() << endl;

    return true;
}


// ************************************************************************* //
