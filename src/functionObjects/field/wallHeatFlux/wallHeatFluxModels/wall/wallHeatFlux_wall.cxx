/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

#include "wallHeatFlux_wall.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatFluxModels
{
    defineTypeNameAndDebug(wall, 0);
    addToRunTimeSelectionTable
    (
        wallHeatFluxModel,
        wall,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::wallHeatFluxModels::wall::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os  << endl;

    writtenHeader_ = true;
}


void Foam::wallHeatFluxModels::wall::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    volScalarField& wallHeatFlux
)
{
    volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

    const auto& heBf = he.boundaryField();

    const auto& alphaBf = alpha.boundaryField();

    for (const label patchi : patchIDs_)
    {
        wallHeatFluxBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
    }


    if (const auto* qrPtr = mesh().cfindObject<volScalarField>(qrName()))
    {
        const auto& radHeatFluxBf = qrPtr->boundaryField();

        for (const label patchi : patchIDs_)
        {
            wallHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatFluxModels::wall::wall
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& name,
    const word objName,
    functionObjects::stateFunctionObject& state
)
:
    wallHeatFluxModel(dict, mesh, name, objName, state)
{
    auto* wallHeatFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                objName,
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    mesh.objectRegistry::store(wallHeatFluxPtr);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::wallHeatFluxModels::wall::read(const dictionary& dict)
{
    if (!wallHeatFluxModel::read(dict))
    {
        return false;
    }


    qrName_ = dict.getOrDefault<word>("qr", "qr");


    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    wordRes patchNames;
    labelHashSet patchSet;
    if (dict.readIfPresent("patches", patchNames) && !patchNames.empty())
    {
        patchSet = pbm.patchSet(patchNames);
    }

    labelHashSet allWalls(pbm.findPatchIDs<wallPolyPatch>());

    Info<< state().type() << " " << state().name() << ":" << nl;

    if (patchSet.empty())
    {
        patchIDs_ = allWalls.sortedToc();

        Info<< "    processing all (" << patchIDs_.size()
            << ") wall patches" << nl << endl;
    }
    else
    {
        allWalls &= patchSet;
        patchSet -= allWalls;
        patchIDs_ = allWalls.sortedToc();

        if (!patchSet.empty())
        {
            WarningInFunction
                << "Requested wall heat-flux on ("
                << patchSet.size() << ") non-wall patches:" << nl;

            for (const label patchi : patchSet.sortedToc())
            {
                Info<< "        " << pbm[patchi].name() << nl;
            }
            Info<< nl;
        }

        Info<< "    processing (" << patchIDs_.size()
            << ") wall patches:" << nl;

        for (const label patchi : patchIDs_)
        {
            Info<< "        " << pbm[patchi].name() << nl;
        }
        Info<< endl;
    }


    if (writeFile::canResetFile())
    {
        writeFile::resetFile(objName());
    }

    if (writeFile::canWriteHeader())
    {
        writeFileHeader(file());
    }


    return true;
}


bool Foam::wallHeatFluxModels::wall::execute()
{
    auto& wallHeatFlux = mesh().lookupObjectRef<volScalarField>(objName());

    if
    (
        const auto* ptr
      = mesh().cfindObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
        ptr
    )
    {
        const auto& turbModel = *ptr;

        calcHeatFlux
        (
            turbModel.alphaEff()(),
            turbModel.transport().he(),
            wallHeatFlux
        );
    }
    else if
    (
        const auto* ptr
      = mesh().cfindObject<fluidThermo>(fluidThermo::dictName);
        ptr
    )
    {
        const auto& thermo = *ptr;

        calcHeatFlux
        (
            thermo.alpha(),
            thermo.he(),
            wallHeatFlux
        );
    }
    else if
    (
        const auto* ptr
      = mesh().cfindObject<solidThermo>(solidThermo::dictName);
        ptr
    )
    {
        const auto& thermo = *ptr;

        calcHeatFlux(thermo.alpha(), thermo.he(), wallHeatFlux);
    }
    else if
    (
        const auto* ptr
      = mesh().cfindObject<basicThermo>("phaseProperties");
        ptr
    )
    {
        const auto& thermo = *ptr;

        calcHeatFlux(thermo.kappaEff()(), thermo.T(), wallHeatFlux);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    const fvPatchList& patches = mesh().boundary();

    const surfaceScalarField::Boundary& magSf = mesh().magSf().boundaryField();

    for (const label patchi : patchIDs_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFlux.boundaryField()[patchi];

        auto limitsHfp = gMinMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << limitsHfp.min()
                << token::TAB << limitsHfp.max()
                << token::TAB << integralHfp
                << endl;
        }

        if (state().log)
        {
            Info<< "    min/max/integ(" << pp.name() << ") = "
                << limitsHfp.min() << ", "
                << limitsHfp.max() << ", "
                << integralHfp << endl;
        }

        state().setResult("min(" + pp.name() + ")", limitsHfp.min());
        state().setResult("max(" + pp.name() + ")", limitsHfp.max());
        state().setResult("int(" + pp.name() + ")", integralHfp);
    }


    return true;
}


bool Foam::wallHeatFluxModels::wall::write()
{
    const auto& wallHeatFlux =
        mesh().lookupObject<volScalarField>(objName());

    if (state().log)
    {
        Info<< "    writing field " << wallHeatFlux.name() << endl;
    }

    wallHeatFlux.write();

    return true;
}


// ************************************************************************* //
