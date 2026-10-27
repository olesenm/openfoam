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

#include "distributed.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "dynamicRefineFvMesh.H"
#include "dynamicFvMesh.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "ListOps.H"
#include "cloud.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace balanceModels
{
    defineTypeNameAndDebug(distributed, 0);
    addToRunTimeSelectionTable(balanceModel, distributed, dictionary);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{
//- Largest-remainder transformation of scalar weights into integer
//  weights summing to 'targetSum', preserving the relative ratios of
//  positive entries as closely as possible. Non-positive entries map
//  to zero. Tie-breaking on equal remainders is by ascending index
//  (stable sort by Foam::sortedOrder).
Foam::labelList transformToIntegers
(
    const Foam::scalarField& w,
    const Foam::label targetSum = 1000000   // larger -> smaller rel. error
)
{
    using namespace Foam;

    const label n = w.size();
    labelList out(n, Zero);

    if (n == 0 || targetSum <= 0)
    {
        return out;
    }

    // Sum of positive weights
    scalar sumW = 0;
    for (label i = 0; i < n; ++i)
    {
        if (w[i] > SMALL)
        {
            sumW += w[i];
        }
    }

    if (sumW <= SMALL)
    {
        return out;
    }

    // Base allocation by floor; fractional remainders kept for top-up
    scalarField rem(n, Foam::zero{});
    label assigned = 0;
    for (label i = 0; i < n; ++i)
    {
        const scalar wn = (w[i] > SMALL) ? w[i]/sumW : scalar(0);
        const scalar raw = wn*scalar(targetSum);
        out[i] = label(std::floor(raw));
        rem[i] = raw - scalar(out[i]);
        assigned += out[i];
    }

    // Distribute remaining seats to the largest remainders.
    // Foam::sortedOrder is stable; tie-break is ascending index.
    labelList order;
    Foam::sortedOrder(rem, order, typename UList<scalar>::greater(rem));

    const label left = targetSum - assigned;
    for (label k = 0; k < left; ++k)
    {
        ++out[order[k]];
    }

    return out;
}


// Write processorWeights into both "coeffs" (fallback) and the resolved
// "<method>Coeffs" block. The coeffs-name resolution is generalised for all
// decomposition methods via methodName + "Coeffs".
void setProcessorWeights
(
    Foam::dictionary& decompDict,
    const Foam::labelUList& weights
)
{
    using namespace Foam;

    // Always populate the generic fallback block.
    dictionary& coeffs = decompDict.subDictOrAdd("coeffs");
    coeffs.set("processorWeights", weights);

    if (word methodName; decompDict.readIfPresent("method", methodName))
    {
        word coeffsName(methodName + "Coeffs");

        // ptscotch uses 'scotchCoeffs'
        if (methodName == "ptscotch")
        {
            coeffsName = "scotchCoeffs";
        }

        const bool existed = (decompDict.findDict(coeffsName) != nullptr);

        auto& methodCoeffs = decompDict.subDictOrAdd(coeffsName);

        // Newly created block: seed it from the generic "coeffs" so that
        // required entries kept there are not orphaned (see Orphaning guard).
        if (!existed)
        {
            methodCoeffs.merge(coeffs);
        }

        methodCoeffs.set("processorWeights", weights);
    }
}
}  // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::balanceModels::distributed::distributed
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    balanceModel(mesh, dict)
{
    // Resolve the decomposeParDict. Honour a -decomposeParDict override that
    // the host application registered on the run-time registry (the balance
    // model has no argList access of its own). Fall back to the default
    // system/decomposeParDict location when no override is registered.
    const auto* registeredDictPtr =
        mesh.time().cfindObject<IOdictionary>("decomposeParDict");

    IOobject decomposeParIO
    (
        "decomposeParDict",
        mesh.time().system(),
        mesh.time(),
        (registeredDictPtr ? IOobject::NO_READ : IOobject::READ_IF_PRESENT),
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    const bool hasDecomposeParDict =
        registeredDictPtr
     || decomposeParIO.typeHeaderOk<IOdictionary>(false, false);

    if (!returnReduceAnd(hasDecomposeParDict))
    {
        FatalIOErrorInFunction(decomposeParIO.objectPath())
            << "Cannot find a decomposeParDict for dynamic load balancing."
            << exit(FatalIOError);
    }

    // Use old value as fallback
    distributionDictPtr_.reset
    (
        new IOdictionary(decomposeParIO, registeredDictPtr)
    );

    word methodName;
    if (!distributionDictPtr_->readIfPresent("method", methodName))
    {
        FatalIOErrorInFunction(*distributionDictPtr_)
            << "Required entry 'method' is missing from decomposeParDict."
            << exit(FatalIOError);
    }

    const labelList defaultWeights(UPstream::nProcs(), 1);

    // Warn early if the selected decomposition method does not use
    // processorWeights: DLB cannot rebalance such 'kahip' or 'manual'. 'kahip'
    // uses only cell weights; 'manual' is an explicit prescribed list.
    if (methodName == "kahip" || methodName == "manual")
    {
        WarningInFunction
            << "Decomposition method '" << methodName
            << "' does not use processorWeights; the '" << typeName
            << "' balance model cannot rebalance this case." << nl
            << "    Use scotch, ptscotch, metis, hierarchical or simple"
            << " for dynamic load balancing." << endl;
    }


    // Non-destructive update: write processorWeights into both "coeffs" and
    // the method-specific "<method>Coeffs" block so the decomposer sees the
    // weights regardless of which subdict it reads.
    setProcessorWeights(*distributionDictPtr_, defaultWeights);

    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::balanceModels::distributed::read(const dictionary& dict)
{
    if (!balanceModel::read(dict))
    {
        return false;
    }

    return true;
}


bool Foam::balanceModels::distributed::balance(const scalarField& procWeights)
{
    if (debug)
    {
        Pout<< "[BEFORE BALANCE] Number of cells: " << mesh().nCells() << endl;
    }

    // Check if the size of procWeights matches nProcs
    if (procWeights.size() != UPstream::nProcs())
    {
        FatalErrorInFunction
            << "processorWeights size: " << procWeights.size()
            << ", does not match nProcs: " << UPstream::nProcs() << "."
            << exit(FatalError);
    }


    // Check if weights are valid.
    // procWeights is already a globally-replicated per-rank vector
    // (length nProcs, identical on every rank); use Field-local sum/min
    // here, NOT gSum/gMin which would erroneously multiply by nProcs.
    scalar sumW = sum(procWeights);
    scalar minW = Foam::min(procWeights);
    if (!std::isfinite(sumW) || sumW <= SMALL || minW <= SMALL)
    {
        FatalErrorInFunction
            << "Invalid processorWeights "
            << "(sum:" << sumW << ", min:" << minW << ")."
            << exit(FatalError);
    }


    const word regionName(mesh().name());
    const label nDomains =
        decompositionMethod::nDomains(*distributionDictPtr_, regionName);

    if (nDomains != UPstream::nProcs())
    {
        FatalIOErrorInFunction(*distributionDictPtr_)
            << "The decomposition domain count for region "
            << regionName << " is " << nDomains
            << ", but the active MPI communicator has "
            << UPstream::nProcs() << " ranks." << nl
            << exit(FatalIOError);
    }


    // Update the internal decomposition dictionary with processor weights.
    // Write to both "coeffs" and "<method>Coeffs" so the decomposer sees the
    // weights regardless of which subdict it reads.
    setProcessorWeights
    (
        *distributionDictPtr_,
        transformToIntegers(procWeights)
    );


    // Create decomposition model.
    //
    // Note: a fresh decompositionMethod is constructed on every balance()
    // call by design. Some derived methods (notably ptscotchDecomp) keep
    // a value-copy of the coefficients dictionary internally, so a cached
    // instance would freeze the processorWeights at first construction and
    // miss all subsequent updates from the imbalance metric. If a future
    // method is known to hold coeffsDict_ by const-reference only, this
    // could be cached as a member of distributed.
    auto decomposerPtr_ = decompositionMethod::New
    (
        *distributionDictPtr_,
        regionName
    );


    // Create new decomposition distribution, new processor assignment per cell
    const labelList distribution
    (
        decomposerPtr_->decompose
        (
            mesh(),
            tmp<scalarField>::New()  // cell weights - redundant
        )
    );


    // Detect whether the new decomposition actually moves any cell off its
    // current rank. If every cell stays put on every rank the redistribution
    // is a no-op and the mesh is reported as unchanged.
    bool localMoved = false;
    for (const label procNo : distribution)
    {
        if (procNo != UPstream::myProcNo())
        {
            localMoved = true;
            break;
        }
    }
    if (!returnReduceOr(localMoved))
    {
        return false;
    }


    // Snapshot the global positions of any Lagrangian clouds before the mesh
    // is redistributed. Cloud::autoMap (invoked during the redistribution
    // mapping) relies on these pre-change positions to relocate particles
    // correctly; storing them afterwards would sample an already-changed mesh.
    for (cloud* cloudPtr : mesh().lookupClass<cloud>())
    {
        cloudPtr->storeGlobalPositions();
    }


    // Create mesh distribution engine
    fvMeshDistribute distributor(mesh());


    // Do actual sending/receiving of mesh
    autoPtr<mapDistributePolyMesh> map = distributor.distribute(distribution);

    if (auto* refineMeshPtr = isA_constCast<dynamicRefineFvMesh>(mesh()))
    {
        hexRef8& meshCutter = refineMeshPtr->meshCutter();

        meshCutter.distribute(map());
    }

    // Update mesh-local state after redistribution.
    if (auto* dynamicMeshPtr = isA_constCast<dynamicFvMesh>(mesh()))
    {
        dynamicMeshPtr->distribute(map());
    }


    // Synchronize processor patch neighbour values for all registered fields.
    // Only coupled patches (processor/processorCyclic) need MPI exchange after
    // redistribution. Non-coupled patches must NOT be re-evaluated here as
    // some, e.g. fixedFluxPressure, require solver state not yet available.
    {
        fvMesh& m = mesh();

        // Batch all fields into a single communication round: post every
        // coupled-patch exchange, wait once, then copy all received values
        // back. This replaces the previous per-field wait, cutting the number
        // of comm rounds to one.
        auto initCoupled = [](auto fieldPtrs)
        {
            for (auto* fPtr : fieldPtrs)
            {
                auto& bf = fPtr->boundaryFieldRef();
                for (auto& bfp : bf)
                {
                    if (bfp.coupled())
                    {
                        bfp.initEvaluate(Pstream::defaultCommsType);
                    }
                }
            }
        };

        auto evaluateCoupled = [](auto fieldPtrs)
        {
            for (auto* fPtr : fieldPtrs)
            {
                auto& bf = fPtr->boundaryFieldRef();
                for (auto& bfp : bf)
                {
                    if (bfp.coupled())
                    {
                        bfp.evaluate(Pstream::defaultCommsType);
                    }
                }
            }
        };

        const auto startOfRequests = UPstream::nRequests();

        do
        {
            #undef  doLocalCode
            #define doLocalCode(Method)                                        \
            {                                                                  \
                Method(m.lookupClass<volScalarField>());                       \
                Method(m.lookupClass<volVectorField>());                       \
                Method(m.lookupClass<volSphericalTensorField>());              \
                Method(m.lookupClass<volSymmTensorField>());                   \
                Method(m.lookupClass<volTensorField>());                       \
                Method(m.lookupClass<surfaceScalarField>());                   \
                Method(m.lookupClass<surfaceVectorField>());                   \
                Method(m.lookupClass<surfaceSphericalTensorField>());          \
                Method(m.lookupClass<surfaceSymmTensorField>());               \
                Method(m.lookupClass<surfaceTensorField>());                   \
            }

            // Phase 1: post the exchange for every registered field. Surface
            // fields (notably phi) carry processor-patch face values that must
            // also be exchanged after redistribution; otherwise the first
            // pressure-correction step sees stale fluxes on processor patches.
            doLocalCode(initCoupled);

            // Phase 2: complete all posted requests in a single wait.
            UPstream::waitRequests(startOfRequests);

            // Phase 3: copy the received neighbour values into boundary fields.
            doLocalCode(evaluateCoupled);

            #undef doLocalCode
        }
        while (false);
    }

    if (debug)
    {
        Pout<< "[AFTER BALANCE] Number of cells: " << mesh().nCells() << endl;
    }

    return true;
}


// ************************************************************************* //
