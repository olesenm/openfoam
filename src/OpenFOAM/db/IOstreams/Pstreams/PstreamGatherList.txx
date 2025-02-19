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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    values[UPstream::myProcNo(comm)].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use broadcast to distribute the data.

\*---------------------------------------------------------------------------*/

#include "contiguous.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::Pstream::gatherList_tree_algorithm
(
    UList<T>& values,
    const int tag,
    const label communicator
)
{
    if (FOAM_UNLIKELY(!UPstream::is_parallel(communicator)))
    {
        // Nothing to do
        return;
    }
    else
    {
        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        const label myProci = UPstream::myProcNo(communicator);

        // Communication order
        const auto& comms = UPstream::whichCommunication(communicator);
        // if (comms.empty()) return;  // extra safety?
        const auto& myComm = comms[myProci];


        // Local buffer for send/recv of contiguous
        [[maybe_unused]] DynamicList<T> buffer;

        // Presize buffer
        if constexpr (is_contiguous_v<T>)
        {
            label maxCount = 0;

            for (const auto belowID : myComm.below())
            {
                auto count = comms[belowID].allBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            if (myComm.above() >= 0)
            {
                auto count = myComm.allBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            buffer.reserve(maxCount + 1);
        }


        // Receive from my downstairs neighbours
        for (const auto belowID : myComm.below())
        {
            const auto& leaves = comms[belowID].allBelow();

            if constexpr (is_contiguous_v<T>)
            {
                if (leaves.empty())
                {
                    // Receive directly into destination
                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        belowID,
                        values[belowID],
                        tag,
                        communicator
                    );
                }
                else
                {
                    // Receive via intermediate buffer
                    buffer.resize_nocopy(leaves.size() + 1);

                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        belowID,
                        buffer,
                        tag,
                        communicator
                    );

                    label recvIdx(0);
                    values[belowID] = buffer[recvIdx++];

                    for (const auto leafID : leaves)
                    {
                        values[leafID] = buffer[recvIdx++];
                    }
                }
            }
            else
            {
                IPstream fromBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,  // bufsize
                    tag,
                    communicator
                );
                fromBelow >> values[belowID];

                if (debug & 2)
                {
                    Perr<< " received through "
                        << belowID << " data from:" << belowID
                        << " data:" << values[belowID] << endl;
                }

                // Receive from all other processors below belowID
                for (const auto leafID : leaves)
                {
                    fromBelow >> values[leafID];

                    if (debug & 2)
                    {
                        Perr<< " received through "
                            << belowID << " data from:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }

        // Send up from values:
        // - my own value first
        // - all belowLeaves next
        if (myComm.above() >= 0)
        {
            const auto& leaves = myComm.allBelow();

            if (debug & 2)
            {
                Perr<< " sending to " << myComm.above()
                    << " data from me:" << myProci
                    << " data:" << values[myProci] << endl;
            }

            if constexpr (is_contiguous_v<T>)
            {
                if (leaves.empty())
                {
                    // Send directly
                    UOPstream::write
                    (
                        UPstream::commsTypes::scheduled,
                        myComm.above(),
                        values[myProci],
                        tag,
                        communicator
                    );
                }
                else
                {
                    // Send via intermediate buffer
                    buffer.resize_nocopy(leaves.size() + 1);

                    label sendIdx(0);
                    buffer[sendIdx++] = values[myProci];

                    for (const auto leafID : leaves)
                    {
                        buffer[sendIdx++] = values[leafID];
                    }

                    UOPstream::write
                    (
                        UPstream::commsTypes::scheduled,
                        myComm.above(),
                        buffer,
                        tag,
                        communicator
                    );
                }
            }
            else
            {
                OPstream toAbove
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    0,  // bufsize
                    tag,
                    communicator
                );
                toAbove << values[myProci];

                for (const auto leafID : leaves)
                {
                    if (debug & 2)
                    {
                        Perr<< " sending to "
                            << myComm.above() << " data from:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                    toAbove << values[leafID];
                }
            }
        }
    }
}


template<class T>
void Foam::Pstream::scatterList_tree_algorithm
(
    UList<T>& values,
    const int tag,
    const label communicator
)
{
    if (FOAM_UNLIKELY(!UPstream::is_parallel(communicator)))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Apart from the additional size check, the only difference
        // between scatterList() and using broadcast(List<T>&) or a regular
        // scatter(List<T>&) is that processor-local data is skipped.

        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        const label myProci = UPstream::myProcNo(communicator);

        // Communication order
        const auto& comms = UPstream::whichCommunication(communicator);
        // if (comms.empty()) return;  // extra safety?
        const auto& myComm = comms[myProci];


        // Local buffer for send/recv of contiguous
        [[maybe_unused]] DynamicList<T> buffer;

        // Presize buffer
        if constexpr (is_contiguous_v<T>)
        {
            label maxCount = 0;

            if (myComm.above() >= 0)
            {
                auto count = myComm.allNotBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            for (const auto belowID : myComm.below())
            {
                auto count = comms[belowID].allNotBelow().size();
                maxCount = Foam::max(maxCount, count);
            }

            buffer.reserve(maxCount);
        }


        // Receive from up
        if (myComm.above() >= 0)
        {
            const auto& leaves = myComm.allNotBelow();

            if constexpr (is_contiguous_v<T>)
            {
                buffer.resize_nocopy(leaves.size());

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    buffer,
                    tag,
                    communicator
                );

                label recvIdx(0);
                for (const auto leafID : leaves)
                {
                    values[leafID] = buffer[recvIdx++];
                }
            }
            else
            {
                IPstream fromAbove
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    0,  // bufsize
                    tag,
                    communicator
                );

                for (const auto leafID : leaves)
                {
                    fromAbove >> values[leafID];

                    if (debug & 2)
                    {
                        Perr<< " received through "
                            << myComm.above() << " data for:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }

        // Send to my downstairs neighbours
        forAllReverse(myComm.below(), belowI)
        {
            const auto belowID = myComm.below()[belowI];
            const auto& leaves = comms[belowID].allNotBelow();

            if constexpr (is_contiguous_v<T>)
            {
                buffer.resize_nocopy(leaves.size());

                label sendIdx(0);
                for (const auto leafID : leaves)
                {
                    buffer[sendIdx++] = values[leafID];
                }

                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    buffer,
                    tag,
                    communicator
                );
            }
            else
            {
                OPstream toBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,  // bufsize
                    tag,
                    communicator
                );

                // Send data destined for all other processors below belowID
                for (const auto leafID : leaves)
                {
                    toBelow << values[leafID];

                    if (debug & 2)
                    {
                        Perr<< " sent through "
                            << belowID << " data for:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }
    }
}


template<class T>
void Foam::Pstream::gatherList
(
    UList<T>& values,
    [[maybe_unused]] const int tag,
    const label communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (is_contiguous_v<T>)
    {
        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(communicator)))
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << UPstream::nProcs(communicator) << nl
                << Foam::abort(FatalError);
        }

        // In-place gather for contiguous types
        UPstream::mpiGather
        (
            nullptr,
            values.data_bytes(),
            sizeof(T),
            communicator
        );
    }
    else
    {
        Pstream::gatherList_tree_algorithm(values, tag, communicator);
    }
}


template<class T>
void Foam::Pstream::scatterList
(
    UList<T>& values,
    [[maybe_unused]] const int tag,
    const label communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (is_contiguous_v<T>)
    {
        // In-place scatter for contiguous types
        UPstream::mpiScatter
        (
            nullptr,
            values.data_bytes(),
            sizeof(T),
            communicator
        );
    }
    else
    {
        Pstream::scatterList_tree_algorithm(values, tag, communicator);
    }
}


template<class T>
void Foam::Pstream::allGatherList
(
    UList<T>& values,
    [[maybe_unused]] const int tag,
    const label comm
)
{
    if (!UPstream::is_parallel(comm))
    {
        // Nothing to do
        return;
    }
    else if constexpr (is_contiguous_v<T>)
    {
        if (FOAM_UNLIKELY(values.size() < UPstream::nProcs(comm)))
        {
            FatalErrorInFunction
                << "List of values is too small:" << values.size()
                << " vs numProcs:" << UPstream::nProcs(comm) << nl
                << Foam::abort(FatalError);
        }

        UPstream::mpiAllGather(values.data_bytes(), sizeof(T), comm);
    }
    else
    {
        Pstream::gatherList(values, tag, comm);
        Pstream::scatterList(values, tag, comm);
    }
}


// ************************************************************************* //
