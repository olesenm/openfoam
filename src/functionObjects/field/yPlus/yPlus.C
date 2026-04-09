/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "yPlus.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "wallFunctionYPlus.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(yPlus, 0);
    addToRunTimeSelectionTable(functionObject, yPlus, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::yPlus::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "y+ ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::yPlus::yPlus
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    UName_("U"),
    useWallFunction_(true),
    writeFields_(true)  // May change in the future
{
    read(dict);

    writeFileHeader(file());

    volScalarField* yPlusPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );

    regIOobject::store(yPlusPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::yPlus::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        useWallFunction_ = true;
        writeFields_ = true;   // May change in the future

        dict.readIfPresent("useWallFunction", useWallFunction_);
        dict.readIfPresent("writeFields", writeFields_);
        UName_ = dict.getOrDefault<word>("field", "U");

        return true;
    }

    return false;
}


bool Foam::functionObjects::yPlus::execute()
{
    auto& yPlus = lookupObjectRef<volScalarField>(scopedName(typeName));

    if (foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        volScalarField::Boundary& yPlusBf = yPlus.boundaryFieldRef();

        const turbulenceModel& model =
            lookupObject<turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        const fvPatchList& patches = mesh_.boundary();

        forAll(patches, patchi)
        {
            const fvPatch& patch = patches[patchi];

            if (isA<wallFvPatch>(patch))
            {
                yPlusBf[patchi] =
                    wallFunctionYPlus(model, patchi, useWallFunction_, UName_);
            }
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: yPlus will not be calculated" << endl;

        if (postProcess)
        {
            WarningInFunction
                << "Please try to use the solver option -postProcess, e.g.:"
                << " <solver> -postProcess -func yPlus" << endl;
        }

        return false;
    }

    return true;
}


bool Foam::functionObjects::yPlus::write()
{
    const auto& yPlus = obr_.lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << ' ' << name() << " write:" << nl;

    if (writeFields_)
    {
        Log << "    writing field " << yPlus.name() << endl;
        yPlus.write();
    }

    const volScalarField::Boundary& yPlusBf = yPlus.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& yPlusp = yPlusBf[patchi];

            auto limits = gMinMax(yPlusp);
            auto avg = gAverage(yPlusp);

            if (UPstream::master())
            {
                writeCurrentTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << limits.min()
                    << token::TAB << limits.max()
                    << token::TAB << avg
                    << endl;
            }

            Log << "    patch " << patch.name()
                << " y+ : min = " << limits.min()
                << ", max = " << limits.max()
                << ", average = " << avg << endl;
        }
    }

    return true;
}


// ************************************************************************* //
