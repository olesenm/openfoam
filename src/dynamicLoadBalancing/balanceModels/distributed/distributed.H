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
    Foam::balanceModels::distributed

Description
    The \c distributed balance model balances the mesh and field distribution
    by computing a new decomposition distribution from the supplied cell weights
    and triggering the processor distributor to send and receive the mesh and
    field parts.

    The new distribution is obtained by injecting per-processor target loads as
    the \c processorWeights coefficient in the decomposition dictionary and
    re-running the selected decomposition method. The decomposition dictionary
    is the one selected on the command line via \c -decomposeParDict, falling
    back to \c system/decomposeParDict when the option is not used:
    \table
      Method               | processorWeights  | Load balancing
      scotch / ptscotch    | yes               | supported
      metis                | yes               | supported
      hierarchical         | yes               | supported
      simple               | yes (approximate) | supported
      kahip                | no                | not supported
      manual               | n/a               | not supported
    \endtable

Usage
    Example fragment in \c constant/dynamicLoadBalanceDict:
    \verbatim
    balanceModel  distributed;
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                  | Type | Reqd | Deflt
      balanceModel | Model name: distributed      | word | yes  | -
    \endtable

SourceFiles
    distributed.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_balanceModels_distributed_H
#define Foam_balanceModels_distributed_H

#include "balanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace balanceModels
{

/*---------------------------------------------------------------------------*\
                            Class distributed Declaration
\*---------------------------------------------------------------------------*/

class distributed
:
    public balanceModel
{
    // Private Data

        //- Local dictionary for run-time distribution parameters
        autoPtr<IOdictionary> distributionDictPtr_;


public:

    //- Runtime type information
    TypeName("distributed");


    // Constructors

        //- Construct from components
        distributed
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~distributed() = default;


    // Evaluation

        //- Balance the load by redistributing the mesh and fields.
        //  Returns true if the redistribution moved any cell, false if the
        //  decomposition left every cell on its current rank.
        virtual bool balance(const scalarField& weights);


    // I-O

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace balanceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
