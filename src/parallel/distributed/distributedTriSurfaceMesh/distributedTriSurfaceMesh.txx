/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "distributedTriSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::distributedTriSurfaceMesh::distributeFields
(
    const mapDistribute& map
)
{
    typedef DimensionedField<Type, triSurfaceGeoMesh> fieldType;

    for (fieldType& field : objectRegistry::sorted<fieldType>())
    {
        const label oldSize = field.size();

        map.distribute(field);

        DebugInfo
            << "Mapped " << field.typeName << ' ' << field.name()
            << " from size " << oldSize << " to size " << field.size() << endl;
    }
}


// ************************************************************************* //
