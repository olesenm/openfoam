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

Class
    Foam::performanceMetrics::cpuLoad

Description
    The \c cpuLoad performance metric measures processor imbalance from
    the ratio of CPU time to wall-clock time per rank.

Usage
    Example fragment in \c constant/dynamicLoadBalanceDict:
    \verbatim
    performanceMetric  cpuLoad;
    \endverbatim

    where the entries mean:
    \table
      Property          | Description                | Type | Reqd | Deflt
      performanceMetric | Model name: cpuLoad        | word | yes  | -
    \endtable

Note
 - Using \c cpuLoad and the \c parProfiling function object will cause
   issues in time measurements due to the singleton nature of the
   \c profilingPstream class.

SourceFiles
    cpuLoad.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_performanceMetrics_cpuLoad_H
#define Foam_performanceMetrics_cpuLoad_H

#include "performanceMetric.H"
#include "cpuTime.H"
#include "clock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{

/*---------------------------------------------------------------------------*\
                            Class cpuLoad Declaration
\*---------------------------------------------------------------------------*/

class cpuLoad
:
    public performanceMetric
{
    // Private Data

        //- Clock object to measure wall-clock time
        //  Local resets are required, hence no use of Time object
        mutable clock clock_;

        //- CPU time object to measure CPU time
        //  Local resets are required, hence no use of Time object
        mutable cpuTime cpuTime_;

        //- Smoothed previous usage value (sentinel -1 = no prior sample)
        mutable scalar prevUsage_;


public:

    //- Runtime type information
    TypeName("cpuLoad");


    // Constructors

        //- Construct from mesh and dictionary
        cpuLoad
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cpuLoad() = default;


    // Member Functions

    // Evaluation

        //- Sample local usage (updates previous samples)
        scalar sampleLocalUsage() const;

        //- Return true after resetting time counters
        //  NOTE: This may need to be refactored if more complex actions are
        //  needed later
        virtual bool resetTimeCounters();


    // I/O

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace performanceMetrics
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
