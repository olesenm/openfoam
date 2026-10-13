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

#include "cellNumber.H"
#include "Pstream.H"
#include "profilingPstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{
    defineTypeNameAndDebug(cellNumber, 0);
    addToRunTimeSelectionTable
    (
        performanceMetric,
        cellNumber,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::performanceMetrics::cellNumber::cellNumber
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    performanceMetric(mesh, dict),
    prevUsage_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::performanceMetrics::cellNumber::read(const dictionary& dict)
{
    if (!performanceMetric::read(dict))
    {
        return false;
    }

    return true;
}


Foam::scalar Foam::performanceMetrics::cellNumber::sampleLocalUsage() const
{
    // Raw usage is the local number of cells
    scalar usage = scalar(mesh().nCells());

    // Smooth usage to avoid large jumps. prevUsage_ < 0 is the sentinel
    // for "no prior sample yet" (first call, or just after
    // resetTimeCounters()).
    if (smoothing() > SMALL && smoothing() <= 1 && prevUsage_ >= 0)
    {
        usage = smoothing()*usage + (1 - smoothing())*prevUsage_;
    }
    prevUsage_ = usage;

    // Clamp to at least one cell
    usage = Foam::max(usage, scalar(1));
    return usage;
}


bool Foam::performanceMetrics::cellNumber::resetTimeCounters()
{
    profilingPstream::reset();

    // Drop the smoothing history so a new "case epoch" starts clean.
    prevUsage_ = -1;

    return true;
}


Foam::tmp<Foam::scalarField>
Foam::performanceMetrics::cellNumber::processorWeights() const
{
    return tmp<scalarField>::New(UPstream::nProcs(), scalar(1));
}


// ************************************************************************* //
