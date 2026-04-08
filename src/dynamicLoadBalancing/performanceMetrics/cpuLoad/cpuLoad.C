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

#include "cpuLoad.H"
#include "Pstream.H"
#include "profilingPstream.H"
#include <algorithm>
#include <limits>
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{
    defineTypeNameAndDebug(cpuLoad, 0);
    addToRunTimeSelectionTable
    (
        performanceMetric,
        cpuLoad,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::performanceMetrics::cpuLoad::cpuLoad
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    performanceMetric(mesh, dict),
    clock_(),
    cpuTime_(),
    prevUsage_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::performanceMetrics::cpuLoad::read(const dictionary& dict)
{
    if (!performanceMetric::read(dict))
    {
        return false;
    }

    return true;
}


Foam::scalar Foam::performanceMetrics::cpuLoad::sampleLocalUsage() const
{
    // Measure wall-clock time increment since the previous restart
    scalar dLocalWall = clock_.clockTimeIncrement();

    // A gate where all processors wait the slowest one
    scalar dWall = returnReduce(dLocalWall, maxOp<scalar>());


    // Measure the simplistic MPI profiling times after all the processors
    // have reached this point thanks to the returnReduce above
    scalar dMpi = 0;
    const profilingPstream::timingList& mpiTimes = profilingPstream::times();

    for (unsigned i = 0; i < profilingPstream::nCategories; ++i)
    {
        dMpi += mpiTimes[i];
    }


    // Measure CPU time increment since the previous restart; ideally, the
    // dCpu includes the MPI profiling times as well
    const scalar dCpuAndMpi  = cpuTime_.cpuTimeIncrement();

    // Estimate the CPU time increment excluding the MPI times
    scalar dCpu = Foam::max(dCpuAndMpi - dMpi, 0.0);

    // Avoid problems with zero wall-clock time increments
    if (dWall <= std::numeric_limits<scalar>::epsilon())
    {
        dWall = 1.0;
        dCpu  = 0.0;
    }


    scalar usage = dCpu/dWall;


    // Smooth usage to avoid large jumps. prevUsage_ < 0 is the sentinel
    // for "no prior sample yet" (first call, or just after
    // resetTimeCounters()).
    if (smoothing() > SMALL && smoothing() <= 1 && prevUsage_ >= 0)
    {
        usage = smoothing()*usage + (1 - smoothing())*prevUsage_;
    }
    prevUsage_ = usage;

    DebugInFunction
        << "Clock time increment: " << dLocalWall << nl
        << "Global wall time increment: " << dWall << nl
        << "CPU time & MPI time increment: " << dCpuAndMpi << nl
        << "MPI time increment: " << dMpi << nl
        << "CPU time increment: " << dCpu << nl
        << "Local usage: " << usage << nl;

    // Clamp to minimum
    usage = Foam::max(usage, scalar(minUsage()));
    return usage;
}


bool Foam::performanceMetrics::cpuLoad::resetTimeCounters()
{
    // The 'clock' is reset internally in clockTimeIncrement()

    cpuTime_.resetCpuTimeIncrement();

    profilingPstream::reset();

    // Drop the smoothing history so a new "case epoch" starts clean.
    prevUsage_ = -1;

    return true;
}


// ************************************************************************* //
