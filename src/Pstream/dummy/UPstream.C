/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::UPstream::vendorCodes Foam::UPstream::vendorCode() noexcept
{
    return UPstream::vendorCodes::invalid;
}


void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{}


bool Foam::UPstream::initNull()
{
    WarningInFunction
        << "The dummy Pstream library cannot be used in parallel mode"
        << endl;

    return false;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    FatalErrorInFunction
        << "The dummy Pstream library cannot be used in parallel mode"
        << endl
        << Foam::exit(FatalError);

    return false;
}


bool Foam::UPstream::setSharedMemoryCommunicators()
{
    return false;
}


void Foam::UPstream::shutdown(int errNo)
{}


void Foam::UPstream::exit(int errNo)
{
    // No MPI - just exit
    std::exit(errNo);
}


void Foam::UPstream::abort(int errNo)
{
    // No MPI - just abort
    std::abort();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::allocateCommunicatorComponents
(
    const label parentIndex,
    const label index
)
{}


void Foam::UPstream::dupCommunicatorComponents
(
    const label parentIndex,
    const label index
)
{}

void Foam::UPstream::splitCommunicatorComponents
(
    const label parentIndex,
    const label index,
    int colour,
    const bool two_step
)
{}


void Foam::UPstream::freeCommunicatorComponents(const label index)
{}


void Foam::UPstream::barrier(const int communicator, UPstream::Request* req)
{}


void Foam::UPstream::send_done
(
    const int toProc,
    const int communicator,
    const int tag
)
{}


int Foam::UPstream::wait_done
(
    const int fromProc,
    const int communicator,
    const int tag
)
{
    return -1;
}


std::pair<int,int64_t>
Foam::UPstream::probeMessage
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    const int tag,
    const int communicator
)
{
    return {-1, 0};
}


void Foam::UPstream::probeMessages
(
    DynamicList<int64_t>& messageSizes,
    const UList<int>& fromProcs,
    const int tag,
    const int communicator
)
{
    messageSizes.clear();
}


// ************************************************************************* //
