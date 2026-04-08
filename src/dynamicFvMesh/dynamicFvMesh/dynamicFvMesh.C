/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "dynamicFvMesh.H"
#include "dynamicLoadBalancingModel.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicFvMesh, 0);
    defineRunTimeSelectionTable(dynamicFvMesh, IOobject);
    defineRunTimeSelectionTable(dynamicFvMesh, doInit);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicFvMesh::readDict()
{
    IOobject dictHeader
    (
        "dynamicMeshDict",
        thisDb().time().constant(),
        thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (dictHeader.typeHeaderOk<IOdictionary>(false, false))
    {
        IOdictionary dict(dictHeader);
        timeControl_.read(dict);

        if (!timeControl_.always())
        {
            // Feedback about the trigger mechanism
            Info<< "Controlled mesh update triggered on "
                << timeControl_.type() << " interval "
                << timeControl_.interval() << nl;
        }
    }
    else
    {
        // Ensure it is indeed pass-through
        timeControl_.clear();
    }
}


void Foam::dynamicFvMesh::readDynamicLoadBalancingModel()
{
    if (!dlbModelInitialised_)
    {
        autoPtr<dynamicLoadBalancingModel> modelPtr =
            dynamicLoadBalancingModel::New(*this);

        dlbModelPtr_.reset(modelPtr.ptr());
        dlbModelInitialised_ = true;
    }
}


bool Foam::dynamicFvMesh::executeDynamicLoadBalancingModel()
{
    readDynamicLoadBalancingModel();

    if (!dlbModelPtr_)
    {
        return false;
    }

    return dlbModelPtr_->execute();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicFvMesh::dynamicFvMesh(const IOobject& io, const bool doInit)
:
    fvMesh(io, doInit),
    timeControl_(io.time(), "update"),   // assume has no side effects
    dlbModelPtr_(nullptr),
    dlbModelInitialised_(false)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        fvMesh::init(doInit);
    }

    readDict();

    // Assume something changed
    return true;
}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    const Foam::zero,
    bool syncPar
)
:
    fvMesh(io, Foam::zero{}, syncPar),
    timeControl_(io.time(), "update"),
    dlbModelPtr_(nullptr),
    dlbModelInitialised_(false)
{
    readDict();
}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    timeControl_(io.time(), "update"),
    dlbModelPtr_(nullptr),
    dlbModelInitialised_(false)
{
    readDict();
}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(cells),
        syncPar
    ),
    timeControl_(io.time(), "update"),
    dlbModelPtr_(nullptr),
    dlbModelInitialised_(false)
{
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::dynamicFvMesh::~dynamicFvMesh()
{}  // dynamicLoadBalancingModel was forward declared


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::dynamicFvMesh::controlledUpdate()
{
    bool meshUpdated = false;

    if (timeControl_.execute())
    {
        if (!timeControl_.always())
        {
            // Feedback that update has been triggered
            Info<< "Mesh update triggered based on "
                << timeControl_.type() << nl;
        }

        addProfiling(mesh, "mesh.update()");
        meshUpdated = this->update();
    }

    const bool dlbModelUpdated = executeDynamicLoadBalancingModel();

    return meshUpdated || dlbModelUpdated;
}


// ************************************************************************* //
