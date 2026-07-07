/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "processorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(processorField, 0);
    addToRunTimeSelectionTable(functionObject, processorField, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::processorField::processorField
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    procFieldPtr_(nullptr)
{
    read(dict);

    createProcessorIDField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::processorField::createProcessorIDField()
{
    if (procFieldPtr_)
    {
        return;
    }

    procFieldPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "processorID",
                mesh_.time().timeName(),
                mesh_,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, UPstream::myProcNo())
        )
    );
}


bool Foam::functionObjects::processorField::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::processorField::execute()
{
    createProcessorIDField();

    *procFieldPtr_ == dimensionedScalar(dimless, UPstream::myProcNo());

    return true;
}


bool Foam::functionObjects::processorField::write()
{
    if (!procFieldPtr_)
    {
        execute();
    }

    procFieldPtr_().write();

    return true;
}


void Foam::functionObjects::processorField::updateMesh(const mapPolyMesh&)
{
    procFieldPtr_.clear();
}


// ************************************************************************* //
