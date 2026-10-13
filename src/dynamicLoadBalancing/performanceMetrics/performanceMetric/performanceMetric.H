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

Namespace
    Foam::performanceMetrics

Description
    A namespace for various performance model implementations for the dynamic
    load balancing.

Class
    Foam::performanceMetric

Description
    A base class for the performance models.

SourceFiles
    performanceMetric.C
    performanceMetricNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_performanceMetric_H
#define Foam_performanceMetric_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
template<class Type> class Function1;

/*---------------------------------------------------------------------------*\
                       Class performanceMetric Declaration
\*---------------------------------------------------------------------------*/

class performanceMetric
{
    // Private Data

    // Input

        //- Const reference to the finite-volume mesh
        const fvMesh& mesh_;

        //- Threshold value that defines when load balancing should be triggered
        autoPtr<Function1<scalar>> imbalanceThresholdPtr_;

        //- Minimum usage (to avoid zero weights)
        scalar minUsage_;

        //- Smoothing factor for usage
        scalar smoothing_;


    // Output

        //- Local usage estimate value
        mutable scalar localUsage_;

        //- Current imbalance estimate value
        mutable scalar imbalanceValue_;

        //- Previous imbalance estimate value
        mutable scalar prevImbalanceValue_;


public:

    //- Runtime type information
    TypeName("performanceMetric");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            performanceMetric,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected performance model
        static autoPtr<performanceMetric> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from mesh and dictionary
        performanceMetric
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        performanceMetric(const performanceMetric&) = delete;

        //- No copy assignment
        void operator=(const performanceMetric&) = delete;


    //- Destructor
    virtual ~performanceMetric();


    // Member Functions

    // Access - input

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept { return mesh_; }

        //- Return the minimum usage
        scalar minUsage() const noexcept { return minUsage_; }

        //- Return the smoothing factor
        scalar smoothing() const noexcept { return smoothing_; }


    // Access - output

        //- Return the local usage value
        scalar localUsage() const noexcept { return localUsage_; }

        //- Return the model value of the performance imbalance
        //  Scaled [0,1], where 0 means perfectly balanced
        scalar imbalanceValue() const noexcept { return imbalanceValue_; }

        //- Return the previous imbalance value
        scalar prevImbalanceValue() const noexcept
        {
            return prevImbalanceValue_;
        }


    // Evaluation

        //- Measure the balance
        scalar calcImbalance() const;

        //- Return true if the performance across processors is deemed
        //- imbalanced
        bool isImbalanced() const;

        //- Return per-processor weights, proportional to inverse usage
        virtual tmp<scalarField> processorWeights() const;


    // Evaluation - to be implemented in derived classes

        //- Sample local usage (updates previous samples)
        virtual scalar sampleLocalUsage() const = 0;

        //- Return true after resetting time counters
        //  NOTE: This may need to be refactored if more complex actions are
        //  needed later
        virtual bool resetTimeCounters() = 0;


    // Bookkeeping

        //- Snapshot the current imbalance value as the previous one.
        //  Must be called by the orchestrator *after* any chronic-imbalance
        //  comparison that reads prevImbalanceValue() so the next time-step
        //  observes a one-step-delayed shadow.
        void updatePrevImbalanceValue() const noexcept
        {
            prevImbalanceValue_ = imbalanceValue_;
        }


    // I-O

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
