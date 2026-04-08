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
    Foam::performanceMetrics::cellNumber

Description
    The \c cellNumber load model measures the discrepancy between the actual
    number of cells on a processor and the ideal distribution.

    It does so by computing the absolute value of the relative deviation
    by |1 - (nCells/idealNcells)|, and it compares each processor's cell count
    to the ideal count, and the maximum difference (whether a processor is
    over-loaded or under-loaded) defines the imbalance metric used for
    load balancing decisions.

Usage
    Example fragment in \c constant/dynamicLoadBalanceDict:
    \verbatim
    performanceMetric  cellNumber;
    \endverbatim

    where the entries mean:
    \table
      Property          | Description              | Type | Reqd | Deflt
      performanceMetric | Model name: cellNumber   | word | yes  | -
    \endtable

SourceFiles
    cellNumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_performanceMetrics_cellNumber_H
#define Foam_performanceMetrics_cellNumber_H

#include "performanceMetric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{

/*---------------------------------------------------------------------------*\
                            Class cellNumber Declaration
\*---------------------------------------------------------------------------*/

class cellNumber
:
    public performanceMetric
{
    // Private Data

        //- Smoothed previous usage value (sentinel -1 = no prior sample)
        mutable scalar prevUsage_;


public:

    //- Runtime type information
    TypeName("cellNumber");


    // Constructors

        //- Construct from mesh and dictionary
        cellNumber
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellNumber() = default;


    // Member Functions

    // Evaluation

        //- Sample local usage (updates previous samples if smoothing > 0)
        scalar sampleLocalUsage() const;

        //- Return true after resetting counters (no-op for this metric)
        virtual bool resetTimeCounters();

        //- Return per-processor weights, uniform across processors
        //- (equal target cell count per rank)
        virtual tmp<scalarField> processorWeights() const;


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
