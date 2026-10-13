/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "dynamicLoadBalancingModel.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicLoadBalancingModel>
Foam::dynamicLoadBalancingModel::New(fvMesh& mesh)
{
    IOobject dictHeader
    (
        "dynamicLoadBalanceDict",
        mesh.time().constant(),
        mesh,
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    if (!dictHeader.typeHeaderOk<IOdictionary>(false, false))
    {
        return nullptr;
    }

    dictionary dict(IOdictionary::readContents(dictHeader));

    if (dict.empty() || !dict.getOrDefault<bool>("active", true))
    {
        return nullptr;
    }

    dictionary modelDict(dict);
    const word modelType
    (
        modelDict.getOrDefault<word>("type", "dynamicLoadBalancer")
    );
    modelDict.set("type", modelType);

    mesh.time().libs().open("libs", modelDict, true);

    if (modelType == "dynamicLoadBalancer")
    {
        mesh.time().libs().open
        (
            {dlLibraryTable::fullname("dynamicLoadBalancing")},
            true
        );
    }

    Info<< "Selecting dynamic-load-balancing model " << modelType << endl;

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "dynamicLoadBalancingModel table is empty" << exit(FatalError);
    }

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            modelDict,
            "type",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<dynamicLoadBalancingModel>(ctorPtr(mesh, modelDict));
}


// ************************************************************************* //
