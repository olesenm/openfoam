/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "randomDecomp.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(randomDecomp);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        randomDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

//- Random distribution on the 0-nCells interval
Foam::labelList randomMap
(
    const Foam::label numBins,
    const Foam::label nCells,
    const Foam::label agglomerate = 0
)
{
    using namespace Foam;

    // The upper sample limit
    const label upper = (numBins - 1);

    if (upper < 1)
    {
        return labelList(nCells, Foam::zero{});
    }

    labelList distribution(nCells);

    Random rndGen(0);

    if (agglomerate > 1)
    {
        label cached = 0;
        label repeat = 0;

        for (label& val : distribution)
        {
            cached = rndGen.position(0, upper);

            if (!repeat)
            {
                cached = rndGen.position<label>(0, upper);
                repeat = agglomerate;
            }
            --repeat;

            val = cached;
        }
    }
    else
    {
        for (label& val : distribution)
        {
            val = rndGen.position<label>(0, upper);
        }
    }

    return distribution;
}

} // End anonyous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomDecomp::randomDecomp(label numDomains)
:
    decompositionMethod(numDomains),
    agglom_(0)
{}


Foam::randomDecomp::randomDecomp
(
    const dictionary& decompDict,
    const word& regionName,
    int select
)
:
    decompositionMethod(decompDict, regionName),
    agglom_(0)
{
    const dictionary& coeffs = findCoeffsDict(typeName + "Coeffs", select);

    // No sanity check needed here (done in randomMap routine)
    coeffs.readIfPresent("agglom", agglom_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::randomDecomp::decompose
(
    const pointField& points,
    const scalarField&
) const
{
    return randomMap(nDomains_, points.size(), agglom_);
}


Foam::labelList Foam::randomDecomp::decompose
(
    const polyMesh& mesh,
    const pointField&,
    const scalarField&
) const
{
    return randomMap(nDomains_, mesh.nCells(), agglom_);
}


Foam::labelList Foam::randomDecomp::decompose
(
    const CompactListList<label>& globalCellCells,
    const pointField&,
    const scalarField&
) const
{
    return randomMap(nDomains_, globalCellCells.size(), agglom_);
}


Foam::labelList Foam::randomDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField&,
    const scalarField&
) const
{
    return randomMap(nDomains_, globalCellCells.size(), agglom_);
}


// ************************************************************************* //
