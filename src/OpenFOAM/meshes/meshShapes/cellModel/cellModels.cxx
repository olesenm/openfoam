/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "cellModel.H"
#include "etcFiles.H"
#include "IFstream.H"
#include "HashSet.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum<Foam::cellModel::modelType> Foam::cellModel::modelNames
({
    { modelType::UNKNOWN, "unknown" },
    { modelType::HEX, "hex" },
    { modelType::WEDGE, "wedge" },
    { modelType::PRISM, "prism" },
    { modelType::PYR, "pyr" },
    { modelType::TET, "tet" },
    { modelType::TETWEDGE, "tetWedge" },
    { modelType::SPLITHEX, "splitHex" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// List of predefined models
Foam::PtrList<Foam::cellModel> models_;

// Lookup of model pointers (in models_) by index
Foam::List<const Foam::cellModel*> modelPtrs_;


// Construct from central "etc/cellModels" file.
void constructModels(bool verbose = false)
{
    using namespace Foam;

    if (models_.size())
    {
        FatalErrorInFunction
            << "attempt to re-construct cellModels when it already exists"
            << exit(FatalError);
    }

    IFstream is(findEtcFile("cellModels", true));  // Mandatory file

    PtrList<cellModel> newPtrs(is);
    models_.swap(newPtrs);

    // Build two lookups: by index, by name
    // Since there are relatively few models, use straight lookup for the index
    // and a linear (non-hashed) search for the name.
    // Lookup by name is less likely than lookup by enum anyhow.

    label maxIndex = -1;
    forAll(models_, i)
    {
        if (maxIndex < models_[i].index())
        {
            maxIndex = models_[i].index();
        }
    }

    if (verbose)
    {
        InfoErr
            << "loaded " << models_.size()
            << " cell models from etc/controlDict"
            << ", max index=" << maxIndex << endl;
    }

    modelPtrs_.clear();
    modelPtrs_.resize(maxIndex+1, nullptr);

    wordHashSet used;
    used.reserve(models_.size());

    forAll(models_, i)
    {
        const cellModel* p = models_.get(i);
        if (!p) continue;

        const label modelIndex = p->index();
        const word& modelName = p->name();

        if (!used.insert(modelName))
        {
            FatalErrorInFunction
                << "more than one model share the name "
                << modelName
                << exit(FatalError);
        }
        else
        {
            if (modelPtrs_[modelIndex])
            {
                FatalErrorInFunction
                    << "more than one model share the index "
                    << modelIndex
                    << exit(FatalError);
            }

            modelPtrs_[modelIndex] = p;
        }
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::cellModel* Foam::cellModel::ptr(const modelType model)
{
    return ptr(label(model));
}


const Foam::cellModel* Foam::cellModel::ptr(const word& modelName)
{
    if (models_.empty())
    {
        constructModels();
    }

    const label n = models_.size();
    for (label i = 0; i < n; ++i)
    {
        const cellModel* p = models_.get(i);

        if (p && p->name() == modelName)
        {
            return p;
        }
    }

    return nullptr;
}


const Foam::cellModel* Foam::cellModel::ptr(const label modelIndex)
{
    if (models_.empty())
    {
        constructModels();
    }

    const cellModel* p(nullptr);

    if (modelIndex >= 0 && modelIndex < modelPtrs_.size())
    {
        p = modelPtrs_[modelIndex];
    }

    return p;
}


const Foam::cellModel& Foam::cellModel::ref(const modelType model)
{
    const cellModel* p = ptr(model);

    if (!p)
    {
        FatalErrorInFunction
            << "No such cellModel: " << modelNames[model]
            << exit(FatalError);
    }

    return *p;
}


const Foam::cellModel& Foam::cellModel::ref(const word& modelName)
{
    const cellModel* p = ptr(modelName);

    if (!p)
    {
        FatalErrorInFunction
            << "No such cellModel: " << modelName
            << exit(FatalError);
    }

    return *p;
}


const Foam::cellModel& Foam::cellModel::ref(const label modelIndex)
{
    const cellModel* p = ptr(modelIndex);

    if (!p)
    {
        FatalErrorInFunction
            << "No such cellModel: " << modelIndex
            << exit(FatalError);
    }

    return *p;
}


// ************************************************************************* //
