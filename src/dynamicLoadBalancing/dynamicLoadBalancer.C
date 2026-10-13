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

#include "dynamicLoadBalancer.H"
#include "balanceModel.H"
#include "performanceMetric.H"
#include "profilingPstream.H"
#include "Function1.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicLoadBalancer, 0);
    addToRunTimeSelectionTable
    (
        dynamicLoadBalancingModel,
        dynamicLoadBalancer,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::dynamicLoadBalancer::needsUpdate() const
{
    // Requires that  isImbalanced() has been called before

    const scalar t = mesh().time().value();
    const label timeIndex = mesh().time().timeIndex();
    const label freq = Foam::max(label(1), frequencyPtr_().value(t));

    bool skipBalance = false;
    if (skipBalancePtr_)
    {
        skipBalance = bool(skipBalancePtr_().value(t));
    }

    bool needsUpdate =
        (
            UPstream::nProcs() > 1
         && timeIndex > 1
         && curTimeIndex_ != timeIndex
         && (freq <= 1 || (timeIndex % freq == 0))
         && !skipBalance
        );

    if (needsUpdate && checkChronicImbalance_)
    {
        const scalar tol = recomputeImbalanceTolerance_;
        const scalar backoffPeriod = backoffPeriod_;

        // Already reduced across processors
        const scalar currImbl = performanceMetricPtr_->imbalanceValue();
        const scalar prevImbl = performanceMetricPtr_->prevImbalanceValue();

        if (Foam::mag(currImbl - prevImbl) < tol)
        {
            if ((t - lastBalanceTime_) < backoffPeriod)
            {
                Info<< "    Chronic imbalance: skipping balancing" << endl;

                needsUpdate = false;
            }
        }
    }

    return needsUpdate;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicLoadBalancer::dynamicLoadBalancer
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    dynamicLoadBalancingModel(mesh, dict),
    dict_(dict),
    balanceModelPtr_(balanceModel::New(mesh, dict_)),
    performanceMetricPtr_(performanceMetric::New(mesh, dict_)),
    recomputeImbalanceTolerance_(0),
    backoffPeriod_(0),
    checkChronicImbalance_(false),
    lastBalanceTime_(0.0),
    curTimeIndex_(-1)
{
    profilingPstream::enable();
    read(dict_);

    if (!UPstream::parRun())
    {
        WarningInFunction
            << "Dynamic load balancing requires running in parallel."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicLoadBalancer::~dynamicLoadBalancer()
{
    profilingPstream::disable();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicLoadBalancer::read(const dictionary& dict)
{
    dict_ = dict;

    if
    (
        !balanceModelPtr_->read(dict_)
     || !performanceMetricPtr_->read(dict_)
    )
    {
        return false;
    }

    skipBalancePtr_ =
        Function1<label>::NewIfPresent("skipBalance", dict_, &mesh());

    recomputeImbalanceTolerance_ = 0;
    backoffPeriod_ = 0;
    checkChronicImbalance_ = false;

    if
    (
        const auto& chronic = dict_.subOrEmptyDict("chronicImbalanceSettings");
       !chronic.empty()
    )
    {
        const bool hasTol  = chronic.found("recomputeImbalanceTolerance");
        const bool hasBack = chronic.found("backoffPeriod");

        if (hasTol && hasBack)
        {
            recomputeImbalanceTolerance_ =
                chronic.getScalar("recomputeImbalanceTolerance");
            backoffPeriod_ = chronic.getScalar("backoffPeriod");
            checkChronicImbalance_ = true;
        }
        else
        {
            WarningInFunction
                << "chronicImbalanceSettings is partially populated; both "
                << "'recomputeImbalanceTolerance' and 'backoffPeriod' must be "
                << "supplied. Disabling chronic-imbalance checking." << nl
                << "    recomputeImbalanceTolerance present: " << hasTol << nl
                << "    backoffPeriod present:               " << hasBack
                << endl;
        }
    }

    frequencyPtr_ = Function1<label>::New("frequency", dict_, &mesh());

    bool ok = true;
    ok = ok && (frequencyPtr_);

    return ok;
}


bool Foam::dynamicLoadBalancer::balance()
{
    const bool imbalanced = performanceMetricPtr_->isImbalanced();
    const bool due = needsUpdate();
    const bool anyDueImbalanced = returnReduceOr(due && imbalanced);

    if (!anyDueImbalanced)
    {
        Info<< "Dynamic load balancing: not performed." << endl;

        // Snapshot current imbalance for the next-step chronic check.
        performanceMetricPtr_->updatePrevImbalanceValue();
        return false;
    }

    curTimeIndex_ = mesh().time().timeIndex();
    lastBalanceTime_ = mesh().time().value();

    Info<< "Dynamic load balancing: will be performed." << endl;

    auto meshPhi = mesh().setPhi();
    const IOobjectOption::writeOption meshPhiWriteOpt
    (
        meshPhi ? meshPhi.cref().writeOpt() : IOobject::NO_WRITE
    );

    mesh().clearMeshPhi();
    bool balanced =
        balanceModelPtr_->balance(performanceMetricPtr_->processorWeights());
    mesh().createMeshPhi();

    // DLB cleared 'meshPhi' before redistributing a moving mesh and recreated
    // it as a zero field. The next solver operation therefore lacked the
    // swept-volume flux from the mesh motion that had just completed.
    // The distributor preserves old points for moving meshes. Recompute the
    // swept-volume flux on the new partition before the solver resumes.
    if (mesh().moving())
    {
        mesh().updateGeom();
    }

    auto newMeshPhi = mesh().setPhi();
    if (newMeshPhi && meshPhiWriteOpt != IOobject::NO_WRITE)
    {
        newMeshPhi.ref().writeOpt(meshPhiWriteOpt);
    }

    Info<< "Dynamic load balancing: performed." << endl;

    performanceMetricPtr_->resetTimeCounters();

    // Snapshot current imbalance for the next-step chronic check, after
    // the chronic comparison in needsUpdate() has already consumed the
    // previous value.
    performanceMetricPtr_->updatePrevImbalanceValue();

    return balanced;
}


bool Foam::dynamicLoadBalancer::execute()
{
    return balance();
}


// ************************************************************************* //
