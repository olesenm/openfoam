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
    Foam::dynamicLoadBalancer

Description
    Runtime-selectable hook that performs dynamic load balancing of the mesh
    and fields according to the specified performance-metric, load and balance
    models.

Usage
    Minimal example by using \c constant/dynamicLoadBalanceDict:
    \verbatim
        // Mandatory entries
        balanceModel        <word>;  // distributed;
        performanceMetric   <word>;  // cellNumber;
        frequency           <Function1<label>>;   // 10;
        imbalanceThreshold  <Function1<scalar>>;  // 0.1;

        // Optional entries
        active              <bool>;
        skipBalance         <Function1<label>>;
        smoothing           <scalar>;
        minUsage            <scalar>;

        chronicImbalanceSettings
        {
            // Mandatory entries
            recomputeImbalanceTolerance  <scalar>;
            backoffPeriod                <scalar>;
        }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      active       | Enable/disable the hook            | bool | no   | true
      balanceModel | Name of specified balance model    | word | yes  | -
      performanceMetric | Name of specified performance metric | word | yes | -
      frequency | Frequency (in time step/iteration index) when      <!--
            --> the dynamic load balancing is performed | F1<label> | yes | -
      imbalanceThreshold | Threshold value that defines when load balancing <!--
                --> should be triggered              | F1<scalar> | yes | -
      skipBalance | Skip load balancing if true      | F1<label>  | no  | -
      smoothing    | Smoothing factor for usage sampling | scalar | no | -1
      minUsage     | Minimum usage (to avoid zero weights) | scalar | no | 1e-3
      chronicImbalanceSettings | Settings for chronic-imbalance backoff <!--
                   -->                                 | dictionary | no | -
      recomputeImbalanceTolerance | Tolerance for recomputing imbalance <!--
                   -->                                 | scalar | cond | -
      backoffPeriod | Backoff period between balance attempts [s] <!--
                   -->                                 | scalar | cond | -
    \endtable

    Options for the \c balanceModel entry:
    \verbatim
      distributed  | Distributes the mesh in a decentralised manner
    \endverbatim

    Options for the \c performanceMetric entry:
    \verbatim
      cellNumber   | Identifies imbalance based on cell-number distribution
      cpuLoad      | Identifies imbalance based on CPU load distribution
    \endverbatim

Note
  - The dynamic load balancer requires running in parallel.
  - Limited to the finite-volume meshes only.
  - The \c frequency variable ignore values less than 1 and treats them as 1.

SourceFiles
    dynamicLoadBalancer.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dynamicLoadBalancer_H
#define Foam_dynamicLoadBalancer_H

#include "fvMesh.H"
#include "dynamicLoadBalancingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class balanceModel;
class performanceMetric;
template<class Type> class Function1;

/*---------------------------------------------------------------------------*\
                    Class dynamicLoadBalancer Declaration
\*---------------------------------------------------------------------------*/

class dynamicLoadBalancer
:
    public dynamicLoadBalancingModel
{
    // Private Data

        //- Top-level dictionary for the user settings
        dictionary dict_;

        //- Load balance model that manages the load balance across processors
        //- by redistributing the mesh and fields on request
        autoPtr<balanceModel> balanceModelPtr_;

        //- Performance metric that defines the criteria for load balancing
        autoPtr<performanceMetric> performanceMetricPtr_;

        //- Frequency (time step/iteration index) at which the load balancing is
        //- performed
        autoPtr<Function1<label>> frequencyPtr_;

        //- Skip load balancing if true (0=false, non-zero=true)
        //  For chronic imbalance situations, repeated load balancing can be
        //  counterproductive; this option allows to skip balancing even though
        //  the imbalance criteria are met
        autoPtr<Function1<label>> skipBalancePtr_;

        //- Tolerance for recomputing imbalance for chronic imbalance
        //  Only consulted when checkChronicImbalance_ is true.
        scalar recomputeImbalanceTolerance_;

        //- Backoff period for load balancing attempts for chronic imbalance [s]
        //  Only consulted when checkChronicImbalance_ is true.
        scalar backoffPeriod_;

        //- Flag to indicate if chronic imbalance checking is enabled
        bool checkChronicImbalance_;

        //- The time of the last successful load balance attempt
        scalar lastBalanceTime_;

        //- The current time-step/iteration index
        label curTimeIndex_;


    // Private Member Functions

        //- Return true if the dynamic load balancing should be triggered
        bool needsUpdate() const;


public:

    //- Runtime type information
    TypeName("dynamicLoadBalancer");


    // Constructors

        //- Construct from the mesh and dictionary
        dynamicLoadBalancer(fvMesh& mesh, const dictionary& dict);


    // Generated Methods

        //- No copy construct
        dynamicLoadBalancer(const dynamicLoadBalancer&) = delete;

        //- No copy assignment
        void operator=(const dynamicLoadBalancer&) = delete;


    //- Destructor
    virtual ~dynamicLoadBalancer();


    // Member Functions

    // Access

        //- Return reference to the top-level user dictionary
        const dictionary& dict() const noexcept { return dict_; }


    // Evaluation

        //- Perform the dynamic load balancing and return the status
        bool balance();

        //- Execute the hook and return mesh-changed status
        virtual bool execute();


    // I/O

        //- Read the top-level dictionary
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
