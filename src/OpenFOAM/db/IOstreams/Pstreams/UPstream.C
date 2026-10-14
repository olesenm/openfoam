/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

Note
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "debug.H"
#include "registerSwitch.H"
#include "dictionary.H"
#include "SHA1.H"
#include "OSspecific.H"  // for hostName()
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);
}

const Foam::Enum
<
    Foam::UPstream::commsTypes
>
Foam::UPstream::commsTypeNames
({
    { commsTypes::buffered, "buffered" },   // "buffered"
    { commsTypes::scheduled, "scheduled" },
    { commsTypes::nonBlocking, "nonBlocking" },  // "immediate"
    // compatibility names
    { commsTypes::buffered, "blocking" },
});


// * * * * * * * * * * * * * Controls Information  * * * * * * * * * * * * * //

void Foam::UPstream::printNodeCommsControl(Ostream& os)
{
    if (UPstream::nodeCommsControl_ > 0)
    {
        if (UPstream::usingNodeComms(UPstream::worldComm))
        {
            os << "on [";
        }
        else
        {
            os << "off [";
        }
        if (UPstream::nodeCommsMin_ > 2)
        {
            os  << "min=" << UPstream::nodeCommsMin_ << ",";
        }
        os  << "type=";

        //  1: split by hostname [default]
        //  2: split by shared
        //  >=4: (debug/manual) split with given number per node
        if (UPstream::nodeCommsControl_ >= 4)
        {
            os  << UPstream::nodeCommsControl_;
        }
        else if (UPstream::nodeCommsControl_ == 2)
        {
            os  << "shared";
        }
        else
        {
            os  << "host";
        }
        os  << "]";
    }
    else
    {
        os  << "disabled";
    }
    os  << " (" << UPstream::nProcs(UPstream::worldComm) << " ranks, "
        << UPstream::numNodes() << " nodes)";
}


void Foam::UPstream::printTopoControl(Ostream& os)
{
    unsigned count = 0;

    if (UPstream::topologyControl_ > 0)
    {
        #undef  PrintControl
        #define PrintControl(Ctrl, Name)                      \
        if (UPstream::usingTopoControl(topoControls::Ctrl))   \
        {                                                     \
            os << (count++ ? ' ' : '(') << Name;              \
        }

        PrintControl(broadcast, "broadcast");
        PrintControl(reduce, "reduce");
        PrintControl(gather, "gather");
        PrintControl(combine, "combine");
        PrintControl(mapGather, "mapGather");
        PrintControl(gatherList, "gatherList");

        #undef PrintControl
    }

    if (count)
    {
        os << ')';  // End the list
    }
    else
    {
        os << "none";
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UPstream::setParRun(const int nProcs, const bool haveThreads)
{
    parRun_ = (nProcs > 0);
    haveThreads_ = haveThreads;

    label comm = -1;
    labelRange singleProc(1);

    // Redo communicators that were created during static initialisation.
    // When parRun == true, redo with MPI components
    // When parRun == false, just redo in case of future changes

    if (!parRun_)
    {
        // These are already correct from the static initialisation,
        // but just in case of future changes

        // Using (world, self) ordering
        freeCommunicator(UPstream::commSelf());
        freeCommunicator(UPstream::commGlobal());

        // 0: COMM_WORLD : commWorld() / commGlobal()
        comm = newCommunicator(-1, singleProc, false);
        if (comm != UPstream::commGlobal())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-global:" << UPstream::commGlobal()
                << Foam::exit(FatalError);
        }

        // 1: COMM_SELF
        comm = newCommunicator(-2, singleProc, false);
        if (comm != UPstream::commSelf())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-self:" << UPstream::commSelf()
                << Foam::exit(FatalError);
        }

        Pout.prefix().clear();
        Perr.prefix().clear();
    }
    else
    {
        // Redo communicators that were created during static initialisation
        // but this time with MPI components

        // Using (world, self) ordering
        freeCommunicator(UPstream::commSelf());
        freeCommunicator(UPstream::commGlobal());

        // 0: COMM_WORLD : commWorld() / commGlobal()
        comm = newCommunicator(-1, labelRange(nProcs), true);
        if (comm != UPstream::commGlobal())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-global:" << UPstream::commGlobal()
                << Foam::exit(FatalError);
        }

        const int globalRanki = UPstream::myProcNo(UPstream::commGlobal());

        // 1: COMM_SELF
        // - Processor number wrt world communicator
        singleProc.start() = globalRanki;
        comm = newCommunicator(-2, singleProc, true);
        if (comm != UPstream::commSelf())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-self:" << UPstream::commSelf()
                << Foam::exit(FatalError);
        }

        Pout.prefix() = '[' + std::to_string(globalRanki) + "] ";
        Perr.prefix() = Pout.prefix();
    }

    if (debug)
    {
        Perr<< "UPstream::setParRun :"
            << " nProcs:" << nProcs
            << " haveThreads:" << haveThreads
            << endl;
    }
}


Foam::label Foam::UPstream::getAvailableCommIndex(const label parentIndex)
{
    label index;
    if (!freeComms_.empty())
    {
        // LIFO pop
        index = freeComms_.back();
        freeComms_.pop_back();

        // Reset existing
        myProcNo_[index] = -1;
        parentComm_[index] = parentIndex;

        procIDs_[index].clear();
        // Sizing and filling are demand-driven
        linearCommunication_[index].clear();
        treeCommunication_[index].clear();
    }
    else
    {
        // Extend storage
        index = parentComm_.size();

        myProcNo_.push_back(-1);
        parentComm_.push_back(parentIndex);

        procIDs_.emplace_back();
        // Sizing and filling are demand-driven
        linearCommunication_.emplace_back(index);
        treeCommunication_.emplace_back(index);
    }

    // Set the communication pattern
    linearCommunication_[index].linear(true);
    treeCommunication_[index].linear(false);

    return index;
}


Foam::label Foam::UPstream::newCommunicator
(
    const label parentIndex,
    const labelRange& subRanks,
    const bool withComponents
)
{
    const label index = getAvailableCommIndex(parentIndex);

    if (debug)
    {
        Perr<< "Allocate communicator ["
            << index << "] from [" << parentIndex
            << "] ranks : " << subRanks << nl
            << endl;
    }

    // Initially treat as master,
    // overwritten by allocateCommunicatorComponents
    myProcNo_[index] = UPstream::masterNo();
    auto& procIds = procIDs_[index];

    // The selected sub-ranks.
    // - transcribe from label to int
    // - already in monotonic order
    if
    (
        (withComponents && UPstream::parRun())
      ? (parentIndex < 0 || subRanks.contains(myProcNo_[parentIndex]))
      : !subRanks.empty()
    )
    {
        procIds.resize_nocopy(subRanks.size());
        std::iota(procIds.begin(), procIds.end(), subRanks.start());
    }
    else
    {
        // Not involved
        procIds.clear();
    }

    if (withComponents && UPstream::parRun())
    {
        allocateCommunicatorComponents(parentIndex, index);
    }

    return index;
}


Foam::label Foam::UPstream::newCommunicator
(
    const label parentIndex,
    const labelUList& subRanks,
    const bool withComponents
)
{
    const label index = getAvailableCommIndex(parentIndex);

    if (debug)
    {
        Perr<< "Allocate communicator ["
            << index << "] from [" << parentIndex
            << "] ranks : " << flatOutput(subRanks) << nl
            << endl;
    }

    // Initially treat as master,
    // overwritten by allocateCommunicatorComponents
    myProcNo_[index] = UPstream::masterNo();
    auto& procIds = procIDs_[index];

    // The selected sub-ranks.
    // - transcribe from label to int
    // - sort into monotonic order (if needed)
    if
    (
        (withComponents && UPstream::parRun())
      ? (parentIndex < 0 || subRanks.contains(myProcNo_[parentIndex]))
      : !subRanks.empty()
    )
    {
        procIds.resize_nocopy(subRanks.size());

        label count = 0;
        bool monotonicOrder = true;
        for (const auto ranki : subRanks)
        {
            if (ranki < 0)
            {
                continue;
            }
            // Could also flag/ignore out-of-range ranks
            // (ranki >= numProcs)

            if (monotonicOrder && count)
            {
                monotonicOrder = (procIds[count-1] < ranki);
            }

            procIds[count] = ranki;
            ++count;
        }

        if (!monotonicOrder)
        {
            auto last = procIds.begin() + count;
            std::sort(procIds.begin(), last);
            last = std::unique(procIds.begin(), last);
            count = label(last - procIds.begin());
        }

        procIds.resize(count);
    }
    else
    {
        // Not involved
        procIds.clear();
    }

    if (withComponents && UPstream::parRun())
    {
        allocateCommunicatorComponents(parentIndex, index);
    }

    return index;
}


Foam::label Foam::UPstream::dupCommunicator
(
    const label parentIndex
)
{
    #ifdef FULLDEBUG
    if (FOAM_UNLIKELY(parentIndex < 0))
    {
        // Failed sanity check
        FatalErrorInFunction
            << "Attempted to duplicate an invalid communicator: "
            << parentIndex
            << Foam::exit(FatalError);
    }
    #endif

    const label index = getAvailableCommIndex(parentIndex);

    if (debug)
    {
        Perr<< "Duplicate communicator ["
            << index << "] from [" << parentIndex << "]" << endl;
    }

    // Initially treat as unknown,
    // overwritten by dupCommunicatorComponents
    myProcNo_[index] = -1;
    procIDs_[index].clear();

    if (UPstream::parRun())
    {
        dupCommunicatorComponents(parentIndex, index);
    }

    return index;
}


Foam::label Foam::UPstream::splitCommunicator
(
    const label parentIndex,
    const int colour,
    const bool two_step
)
{
    #ifdef FULLDEBUG
    if (FOAM_UNLIKELY(parentIndex < 0))
    {
        // Failed sanity check
        FatalErrorInFunction
            << "Attempted to split an invalid communicator: "
            << parentIndex
            << Foam::exit(FatalError);
    }
    #endif

    const label index = getAvailableCommIndex(parentIndex);

    if (debug)
    {
        Perr<< "Split communicator ["
            << index << "] from [" << parentIndex
            << "] using colour=" << colour
            << " (two_step=" << two_step << ")" << endl;
    }

    // Initially treat as unknown,
    // overwritten by splitCommunicatorComponents
    myProcNo_[index] = -1;
    procIDs_[index].clear();

    if (UPstream::parRun())
    {
        splitCommunicatorComponents(parentIndex, index, colour, two_step);
    }

    return index;
}


bool Foam::UPstream::setHostCommunicators(const int numPerNode)
{
    // Uses the world communicator (not global communicator)

    // Skip if non-parallel
    if (!UPstream::parRun())
    {
        numNodes_ = 1;
        return false;
    }

    if (FOAM_UNLIKELY(commInterNode_ >= 0 || commLocalNode_ >= 0))
    {
        // Failed sanity check
        FatalErrorInFunction
            << "Node communicator(s) already created!" << endl
            << Foam::abort(FatalError);
        return false;
    }

    commInterNode_ = getAvailableCommIndex(constWorldComm_);
    commLocalNode_ = getAvailableCommIndex(constWorldComm_);

    // Overwritten later
    myProcNo_[commInterNode_] = UPstream::masterNo();
    myProcNo_[commLocalNode_] = UPstream::masterNo();

    // Sorted order, purely cosmetic
    if (commLocalNode_ < commInterNode_)
    {
        std::swap(commLocalNode_, commInterNode_);
    }

    if (debug)
    {
        Perr<< "Allocating node communicators "
            << commInterNode_ << ", " << commLocalNode_
            << " on parent : " << constWorldComm_ << nl
            << endl;
    }

    const int worldRank = UPstream::myProcNo(constWorldComm_);
    const int worldSize = UPstream::nProcs(constWorldComm_);

    if (numPerNode > 1)
    {
        // Manual splitting based on given number of ranks per node
        const int myNodeId = (worldRank/numPerNode);

        // Establish the topology
        {
            DynamicList<int> nodeGroup(numPerNode);
            DynamicList<int> nodeLeader(1+worldSize/numPerNode);

            for (int proci = 0; proci < worldSize; ++proci)
            {
                if (myNodeId == (proci/numPerNode))
                {
                    nodeGroup.push_back(proci);
                }

                if ((proci % numPerNode) == 0)
                {
                    // Local rank 0 is a node leader
                    nodeLeader.push_back(proci);
                }
            }

            procIDs_[commInterNode_] = std::move(nodeLeader);
            procIDs_[commLocalNode_] = std::move(nodeGroup);
        }
    }
    else
    {
        // Determine inter-host/inter-host grouping based on the SHA1 of the
        // hostnames. This allows a single initial Allgather to establish
        // the overall topology. The alternative is to use MPI_Split_comm_type()
        // on SHARED and then MPI_Comm_split() on the leader ranks.

        // Could also add lowercase etc, but since hostName()
        // will be consistent within the same node, there is no need.
        const SHA1Digest myDigest(SHA1(hostName()).digest());

        List<SHA1Digest> digests(worldSize);
        digests[worldRank] = myDigest;

        // The fixed-length digest allows use of MPI_Allgather.
        UPstream::mpiAllGather
        (
            digests.data_bytes(),       // Send/Rev
            SHA1Digest::size_bytes(),   // Num send/recv data per rank
            UPstream::constWorldComm_
        );

        // Establish the topology
        {
            DynamicList<int> nodeGroup(64);
            DynamicList<int> nodeLeader(64);
            DynamicList<SHA1Digest> uniqDigests(64);

            for (int proci = 0; proci < worldSize; ++proci)
            {
                const auto& dig = digests[proci];

                if (myDigest == dig)
                {
                    nodeGroup.push_back(proci);
                }

                if (!uniqDigests.contains(dig))
                {
                    // First appearance of host
                    uniqDigests.push_back(dig);
                    nodeLeader.push_back(proci);
                }
            }

            procIDs_[commInterNode_] = std::move(nodeLeader);
            procIDs_[commLocalNode_] = std::move(nodeGroup);
        }
    }


    // Capture the size (number of nodes) before doing anything further
    numNodes_ = procIDs_[commInterNode_].size();

    // ~~~~~~~~~
    // IMPORTANT
    // ~~~~~~~~~
    // Always retain knowledge of the inter-node leaders,
    // even if this process is not on that communicator.
    // This will help when constructing topology-aware communication.

    // Allocate backend MPI components
    allocateCommunicatorComponents(constWorldComm_, commInterNode_);
    allocateCommunicatorComponents(constWorldComm_, commLocalNode_);

    return true;
}


void Foam::UPstream::freeCommunicator
(
    const label communicator,
    const bool withComponents
)
{
    // Filter out any placeholders
    if (communicator < 0)
    {
        return;
    }

    if (debug)
    {
        Perr<< "Communicators : Freeing communicator " << communicator
            << " parent: " << parentComm_[communicator]
            << " myProcNo: " << myProcNo_[communicator]
            << endl;
    }

    if (withComponents && parRun())
    {
        freeCommunicatorComponents(communicator);
    }

    myProcNo_[communicator] = -1;
    parentComm_[communicator] = -1;
    //procIDs_[communicator].clear();
    linearCommunication_[communicator].clear();
    treeCommunication_[communicator].clear();

    // LIFO push
    freeComms_.push_back(communicator);
}


int Foam::UPstream::baseProcNo(label comm, int procID)
{
    while (UPstream::parent(comm) >= 0 && procID >= 0)
    {
        const auto& parentRanks = UPstream::procID(comm);
        procID = parentRanks[procID];
        comm = parent(comm);
    }

    return procID;
}


Foam::label Foam::UPstream::procNo(const label comm, const int baseProcID)
{
    const auto& parentRanks = UPstream::procID(comm);
    label parentComm = UPstream::parent(comm);

    int procID = baseProcID;

    if (parentComm >= 0)
    {
        procID = UPstream::procNo(parentComm, baseProcID);
    }

    return parentRanks.find(procID);
}


Foam::label Foam::UPstream::procNo
(
    const label comm,
    const label currentComm,
    const int currentProcID
)
{
    label physProcID = UPstream::baseProcNo(currentComm, currentProcID);
    return UPstream::procNo(comm, physProcID);
}


const Foam::UPstream::commsStructList&
Foam::UPstream::linearCommunication(int communicator)
{
    // linear = true

    if (linearCommunication_[communicator].empty())
    {
        linearCommunication_[communicator].init(communicator, true);
    }
    // Probably not needed
    // else
    // {
    //     linearCommunication_[communicator].linear(true);
    // }

    return linearCommunication_[communicator];
}


const Foam::UPstream::commsStructList&
Foam::UPstream::treeCommunication(int communicator)
{
    // linear = false

    if (treeCommunication_[communicator].empty())
    {
        treeCommunication_[communicator].init(communicator, false);
    }
    // Probably not needed
    // else
    // {
    //     treeCommunication_[communicator].linear(false);
    // }

    return treeCommunication_[communicator];
}


void Foam::UPstream::printCommTree
(
    int communicator,
    bool linear
)
{
    const auto& comms = UPstream::whichCommunication(communicator, linear);

    if (UPstream::master(communicator))
    {
        comms.printGraph(Info());
    }
}


bool Foam::UPstream::usingNodeComms(const int communicator)
{
    // Starting point must be "real" world-communicator
    // ("real" means without any local trickery with worldComm)
    // Avoid corner cases:
    // - everthing is on one node
    // - everthing is on different nodes

    return
    (
        parRun_ && (constWorldComm_ == communicator)
     && (nodeCommsControl_ > 0)

        // More than one node and above defined threshold
     && (numNodes_ > 1) && (numNodes_ >= nodeCommsMin_)
        // Some processes do share nodes
     && (numNodes_ < procIDs_[constWorldComm_].size())

        // Extra paranoid (guard against calling during startup)
     && (commInterNode_ > constWorldComm_)
     && (commLocalNode_ > constWorldComm_)
    );
}


const Foam::List<int>& Foam::UPstream::interNode_offsets()
{
    static std::unique_ptr<List<int>> singleton;

    if (!singleton)
    {
        // Extra paranoid (guard against calling during startup)
        if
        (
            (commInterNode_ <= constWorldComm_)
         || (commInterNode_ >= procIDs_.size())
        )
        {
            return List<int>::null();
        }

        singleton = std::make_unique<List<int>>();
        auto& offsets = *singleton;

        const auto& procs = procIDs_[commInterNode_];

        // The procIDs_ are already the offsets, but missing the end offset
        if (!procs.empty())
        {
            const auto count = procs.size();

            offsets.resize(count+1);
            std::copy_n
            (
                procs.begin(),
                count,
                offsets.begin()
            );
            offsets[count] = UPstream::nProcs(constWorldComm_);
        }
    }

    return *singleton;
}


const Foam::UPstream::rangeType& Foam::UPstream::localNode_parentProcs()
{
    static UPstream::rangeType singleton;

    if (singleton.empty())
    {
        // The inter-node offsets [in const world comm] also include a
        // startup guard
        const auto& offsets = UPstream::interNode_offsets();

        const auto nodei =
            Foam::findLower
            (
                offsets,
                // My place within const world comm
                UPstream::myProcNo(constWorldComm_)+1
            );

        if (nodei >= 0)
        {
            singleton.reset
            (
                offsets[nodei],
                offsets[nodei+1] - offsets[nodei]
            );
        }
    }

    return singleton;
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::UPstream::parRun_(false);

bool Foam::UPstream::haveThreads_(false);

bool Foam::UPstream::noInitialCommDup_(false);

int Foam::UPstream::msgType_(1);


Foam::wordList Foam::UPstream::allWorlds_(Foam::one{}, "");
Foam::labelList Foam::UPstream::worldIDs_(Foam::one{}, 0);


Foam::DynamicList<int> Foam::UPstream::myProcNo_(16);
Foam::DynamicList<Foam::List<int>> Foam::UPstream::procIDs_(16);

Foam::DynamicList<Foam::label> Foam::UPstream::parentComm_(16);
Foam::DynamicList<Foam::label> Foam::UPstream::freeComms_;

Foam::DynamicList<Foam::UPstream::commsStructList>
Foam::UPstream::linearCommunication_(16);

Foam::DynamicList<Foam::UPstream::commsStructList>
Foam::UPstream::treeCommunication_(16);


int Foam::UPstream::constWorldComm_(0);
int Foam::UPstream::commInterNode_(-1);
int Foam::UPstream::commLocalNode_(-1);
int Foam::UPstream::numNodes_(1);

Foam::label Foam::UPstream::worldComm(0);  // Initially same as constWorldComm_
Foam::label Foam::UPstream::warnComm(-1);


// Predefine world and self communicator slots.
// These are overwritten in parallel mode (by UPstream::setParRun())
const int nPredefinedComm = []()
{
    // 0: COMM_WORLD : commGlobal(), constWorldComm_, worldComm
    (void) Foam::UPstream::newCommunicator(-1, Foam::labelRange(1), false);

    // 1: COMM_SELF
    (void) Foam::UPstream::newCommunicator(-2, Foam::labelRange(1), false);

    return Foam::UPstream::nComms();
}();


int Foam::UPstream::nodeCommsControl_
(
    Foam::debug::optimisationSwitch("nodeComms", 1)
);
registerOptSwitch
(
    "nodeComms",
    int,
    Foam::UPstream::nodeCommsControl_
);

int Foam::UPstream::nodeCommsMin_
(
    Foam::debug::optimisationSwitch("nodeComms.min", 0)
);
registerOptSwitch
(
    "nodeComms.min",
    int,
    Foam::UPstream::nodeCommsMin_
);

int Foam::UPstream::topologyControl_
(
    Foam::debug::optimisationSwitch("topoControl", 0)
);
registerOptSwitch
(
    "topoControl",
    int,
    Foam::UPstream::topologyControl_
);

bool Foam::UPstream::floatTransfer
(
    Foam::debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitch
(
    "floatTransfer",
    bool,
    Foam::UPstream::floatTransfer
);

int Foam::UPstream::nProcsSimpleSum
(
    Foam::debug::optimisationSwitch("nProcsSimpleSum", 0)
);
registerOptSwitch
(
    "nProcsSimpleSum",
    int,
    Foam::UPstream::nProcsSimpleSum
);

int Foam::UPstream::nProcsNonblockingExchange
(
    Foam::debug::optimisationSwitch("nbx.min", 0)
);
registerOptSwitch
(
    "nbx.min",
    int,
    Foam::UPstream::nProcsNonblockingExchange
);

int Foam::UPstream::tuning_NBX_
(
    Foam::debug::optimisationSwitch("nbx.tuning", 0)
);
registerOptSwitch
(
    "nbx.tuning",
    int,
    Foam::UPstream::tuning_NBX_
);


int Foam::UPstream::nPollProcInterfaces
(
    Foam::debug::optimisationSwitch("nPollProcInterfaces", 0)
);
registerOptSwitch
(
    "nPollProcInterfaces",
    int,
    Foam::UPstream::nPollProcInterfaces
);


Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.get
    (
        "commsType",
        Foam::debug::optimisationSwitches()
    )
);


//! \cond file_scope
namespace Foam
{
    //- Registered reader for UPstream::defaultCommsType
    class addcommsTypeToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addcommsTypeToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addcommsTypeToOpt() = default;

        virtual void readData(Foam::Istream& is)
        {
            UPstream::defaultCommsType =
                UPstream::commsTypeNames.read(is);
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os << UPstream::commsTypeNames[UPstream::defaultCommsType];
        }
    };

    addcommsTypeToOpt addcommsTypeToOpt_("commsType");
}
//! \endcond

int Foam::UPstream::maxCommsSize
(
    Foam::debug::optimisationSwitch("maxCommsSize", 0)
);
registerOptSwitch
(
    "maxCommsSize",
    int,
    Foam::UPstream::maxCommsSize
);


const int Foam::UPstream::mpiBufferSize
(
    Foam::debug::optimisationSwitch("mpiBufferSize", 0)
);


int Foam::UPstream::listChunkLength_
(
    Foam::debug::optimisationSwitch("mpi.listChunkLength", 0)
);
registerOptSwitch
(
    "mpi.listChunkLength",
    int,
    Foam::UPstream::listChunkLength_
);


// ************************************************************************* //
