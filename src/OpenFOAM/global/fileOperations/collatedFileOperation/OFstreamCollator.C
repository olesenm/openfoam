/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "OFstreamCollator.H"
#include "OFstream.H"
#include "decomposedBlockData.H"
#include "dictionary.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstreamCollator, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::OFstreamCollator::writeFile
(
    const label comm,
    const word& objectType,
    const fileName& fName,
    const UList<char>& localData,
    const labelUList& recvSizes,
    const UList<std::string_view>& procData,  // optional proc data
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const dictionary& headerEntries
)
{
    if (debug)
    {
        Pout<< "OFstreamCollator : Writing local " << localData.size()
            << " bytes to " << fName << " using comm " << comm
            << " and " << procData.size() << " sub-ranks" << endl;

        forAll(procData, proci)
        {
            Pout<< "    " << proci << " size:"
                << label(procData[proci].size()) << nl;
        }
    }

    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm))
    {
        Foam::mkDir(fName.path());
        osPtr.reset(new OFstream(atomic, fName, streamOpt, append));
        auto& os = *osPtr;

        if (append == IOstreamOption::NO_APPEND)
        {
            // No IOobject so cannot use IOobject::writeHeader

            // FoamFile
            decomposedBlockData::writeHeader
            (
                os,
                streamOpt,      // streamOpt for container
                objectType,
                "",             // note
                "",             // location (leave empty, otherwise inaccurate)
                fName.name(),   // object name
                headerEntries
            );
        }
    }

    // Assuming threaded writing hides any slowness so we
    // can use scheduled communication to send the data to
    // the master processor in order. However can be unstable
    // for some mpi so default is non-blocking.
    const UPstream::commsTypes myCommsType
    (
        // Blocking when buffer size is 0
        Foam::mag
        (
            fileOperations::masterUncollatedFileOperation::
            maxMasterFileBufferSize
        ) < 1
      ? UPstream::commsTypes::scheduled
      : UPstream::commsTypes::nonBlocking
    );


    List<std::streamoff> blockOffsets;  // Optional
    decomposedBlockData::writeBlocks
    (
        comm,
        osPtr,
        blockOffsets,  // or List<std::streamoff>::null()
        localData,
        recvSizes,
        procData,
        myCommsType,
        false       // do not sync return state
    );

    if (osPtr && !osPtr->good())
    {
        FatalIOErrorInFunction(*osPtr)
            << "Failed writing to " << fName << exit(FatalIOError);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Finished writing "
            << localData.size() << " bytes";

        if (UPstream::master(comm))
        {
            off_t total = 0;
            for (const label recv : recvSizes)
            {
                total += recv;
            }
            // Use std::to_string to display long int
            Pout<< " (overall " << std::to_string(total) << ')';
        }
        Pout<< " to " << fName
            << " using comm " << comm << endl;
    }

    return true;
}


void* Foam::OFstreamCollator::writeAll(void *threadarg)
{
    OFstreamCollator& handler = *static_cast<OFstreamCollator*>(threadarg);

    // Consume stack
    while (true)
    {
        std::unique_ptr<writeData> ptr;

        {
            std::lock_guard<std::mutex> guard(handler.mutex_);

            if (handler.objects_.size())
            {
                // FIFO
                ptr.reset(handler.objects_.front());
                handler.objects_.pop_front();
            }
        }

        if (!ptr)
        {
            break;
        }

        writeData& obj = *ptr;

        // Obtain views from storage
        List<std::string_view> procData(obj.procData_.size());
        forAll(procData, proci)
        {
            procData[proci] = obj.procData_[proci].view();
        }

        bool ok = writeFile
        (
            obj.comm_,
            obj.objectType_,
            obj.pathName_,
            obj.localData_,
            obj.sizes_,
            procData,
            obj.streamOpt_,
            obj.atomic_,
            obj.append_,
            obj.headerEntries_
        );

        if (!ok)
        {
            FatalIOErrorInFunction(obj.pathName_)
                << "Failed writing " << obj.pathName_
                << exit(FatalIOError);
        }
        //sleep(1);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Exiting write thread " << endl;
    }

    {
        std::lock_guard<std::mutex> guard(handler.mutex_);
        handler.threadRunning_ = false;
    }

    return nullptr;
}


void Foam::OFstreamCollator::waitForBufferSpace(const off_t wantedSize) const
{
    while (true)
    {
        // The pending output size(s)
        off_t totalSize = 0;

        {
            std::lock_guard<std::mutex> guard(mutex_);
            for (const writeData* ptr : objects_)
            {
                if (ptr) totalSize += ptr->size();
            }
        }

        if
        (
            totalSize == 0
         || (wantedSize >= 0 && (totalSize+wantedSize) <= maxBufferSize_)
        )
        {
            break;
        }

        if (debug)
        {
            std::lock_guard<std::mutex> guard(mutex_);
            Pout<< "OFstreamCollator : Waiting for buffer space."
                << " Currently in use:" << totalSize
                << " limit:" << maxBufferSize_
                << " files:" << objects_.size()
                << endl;
        }

        sleep(5);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstreamCollator::OFstreamCollator(const off_t maxBufferSize)
:
    OFstreamCollator(maxBufferSize, UPstream::worldComm)
{}


Foam::OFstreamCollator::OFstreamCollator
(
    const off_t maxBufferSize,
    const label comm
)
:
    maxBufferSize_(maxBufferSize),
    threadRunning_(false),
    localComm_(comm),
    threadComm_(UPstream::dupCommunicator(localComm_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstreamCollator::~OFstreamCollator()
{
    if (thread_)
    {
        if (debug)
        {
            Pout<< "~OFstreamCollator : Waiting for write thread" << endl;
        }
        thread_->join();
        thread_.reset(nullptr);
    }

    UPstream::freeCommunicator(threadComm_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OFstreamCollator::write
(
    const word& objectType,
    const fileName& fName,
    DynamicList<char>&& localData,
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const bool useThread,
    const dictionary& headerEntries
)
{
    // Determine (on master) sizes to receive. Note: do NOT use thread
    // communicator
    const labelList recvSizes
    (
        UPstream::listGatherValues<label>(localData.size(), localComm_)
    );

    off_t totalSize = 0;
    label maxLocalSize = 0;

    if (UPstream::master(localComm_))
    {
        for (const label recvSize : recvSizes)
        {
            totalSize += recvSize;
            maxLocalSize = Foam::max(maxLocalSize, recvSize);
        }
    }

    // Broadcast the information to everyone
    {
        int64_t sizes[2] =
        {
            static_cast<int64_t>(totalSize),
            static_cast<int64_t>(maxLocalSize)
        };

        UPstream::broadcast(sizes, 2, localComm_);

        totalSize = static_cast<off_t>(sizes[0]);
        maxLocalSize = static_cast<label>(sizes[1]);
    }


    // Determine how things will be gathered and written...

    enum class dispatchModes { DIRECT_WRITE, THREADED_WRITE, FULL_THREADED };

    dispatchModes dispatch(dispatchModes::DIRECT_WRITE);

    if (!useThread || maxBufferSize_ == 0 || maxLocalSize > maxBufferSize_)
    {
        dispatch = dispatchModes::DIRECT_WRITE;
    }
    else if (totalSize <= maxBufferSize_)
    {
        // Total size can be stored locally
        // - gather all data now and only do the writing in the thread

        dispatch = dispatchModes::THREADED_WRITE;
    }
    else
    {
        // Gather data and write in the thread

        dispatch = dispatchModes::FULL_THREADED;

        if (!UPstream::haveThreads())
        {
            WarningInFunction
                << "MPI not initialized with thread support." << nl
                << "    maxThreadFileBufferSize = 0 to disable threading" << nl
                << " or maxThreadFileBufferSize > " << totalSize
                << " to collate before threaded writing." << nl << nl;

            dispatch = dispatchModes::DIRECT_WRITE;
        }
    }


    // -----------
    // Dispatching
    // -----------

    if (dispatch == dispatchModes::DIRECT_WRITE)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather "
                << "(local comm: " << localComm_
                << "); non-thread write of "
                << fName << endl;
        }

        // Direct collating and writing (so master blocks until all written!)
        return writeFile
        (
            localComm_,
            objectType,
            fName,
            localData,
            recvSizes,
            UList<std::string_view>::null(),  // dummy proc data
            streamOpt,
            atomic,
            append,
            headerEntries
        );
    }
    else if (dispatch == dispatchModes::THREADED_WRITE)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather "
                << "(local comm: " << localComm_
                << "); thread write of "
                << fName << endl;
        }

        if (UPstream::master(localComm_))
        {
            waitForBufferSpace(totalSize);
        }

        std::unique_ptr<writeData> fileAndDataPtr
        (
            new writeData
            (
                threadComm_,
                objectType,
                fName,
                recvSizes,
                streamOpt,
                atomic,
                append,
                headerEntries
            )
        );
        auto& fileAndData = *fileAndDataPtr;

        List<List<char>>& procData = fileAndData.procData_;

        // Receive from these procs
        DynamicList<int> recvProcs;

        if (UPstream::master(localComm_))
        {
            // Move in local data (master only!)
            fileAndData.transfer(localData);

            // Storage for receive data
            procData.resize_nocopy(UPstream::nProcs(localComm_));

            // Sorted by message size
            labelList order(Foam::sortedOrder(recvSizes));
            recvProcs.reserve_exact(order.size());

            // Want to receive large messages first. Ignore empty slots
            forAllReverse(order, i)
            {
                const label proci = order[i];

                // Ignore empty slots and don't try to receive from self
                if (recvSizes[proci] > 0 && proci != UPstream::masterNo())
                {
                    recvProcs.push_back(proci);
                }
            }
        }
        else if (UPstream::is_subrank(localComm_))
        {
            // Requires a size for decomposedBlockData::writeBlocks() logic
            procData.resize_nocopy(UPstream::nProcs(localComm_));
        }


        // Gather all data onto master. Is done in local communicator since
        // not in write thread.
        const auto startOfRequests = UPstream::nRequests();

        const int messageTag = (UPstream::msgType() + 256);

        if (UPstream::master(localComm_))
        {
            // Receive from these procs (non-empty slots)
            for (const int proci : recvProcs)
            {
                auto& slot = procData[proci];
                slot.resize_nocopy(recvSizes[proci]);

                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    slot.data_bytes(),
                    slot.size_bytes(),
                    messageTag,
                    localComm_
                );
            }
        }
        else if (UPstream::is_subrank(localComm_) && !localData.empty())
        {
            // Send to content to master
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    UPstream::masterNo(),
                    localData.cdata_bytes(),
                    localData.size_bytes(),
                    messageTag,
                    localComm_
                )
            )
            {
                FatalErrorInFunction
                    << "Failure to send message (size: "
                    << localData.size() << ") to master" << nl
                    << Foam::abort(FatalError);
            }
        }

        UPstream::waitRequests(startOfRequests);

        // The localData has been moved (master) or communicated
        localData.clearStorage();


        // Queue up for threading
        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Add to thread buffer (as FIFO), take ownership
            objects_.push_back(fileAndDataPtr.release());

            // Start thread if not running
            if (!threadRunning_)
            {
                if (thread_)
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_->join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread"
                        << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }
    else if (dispatch == dispatchModes::FULL_THREADED)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : thread gather; thread write "
                << "(thread comm: " << threadComm_
                << ") of " << fName << endl;
        }

        if (UPstream::master(localComm_))
        {
            waitForBufferSpace(localData.size());
        }

        std::unique_ptr<writeData> fileAndDataPtr
        (
            new writeData
            (
                threadComm_,
                objectType,
                fName,
                recvSizes,
                streamOpt,
                atomic,
                append,
                headerEntries
            )
        );

        // Move in local data (all procs)
        fileAndDataPtr->transfer(localData);


        // Queue up for threading
        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Add to thread buffer (as FIFO), take ownership
            objects_.push_back(fileAndDataPtr.release());

            // Note: no proc data provided
            // so it will trigger communication inside the thread!!!

            if (!threadRunning_)
            {
                if (thread_)
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_->join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread" << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }

    FatalErrorInFunction
        << "Unknown dispatch mode: " << int(dispatch)
        << " - programming error?" << abort(FatalError);

    return false;
}


void Foam::OFstreamCollator::waitAll()
{
    // Wait for all buffer space to be available
    // - ie, wait for all jobs to finish

    if (UPstream::master(localComm_))
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : waiting for thread to have consumed all"
                << endl;
        }
        waitForBufferSpace(-1);
    }
}


// ************************************************************************* //
