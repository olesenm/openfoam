/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "contiguous.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Single values to/from a list

template<class T>
Foam::List<T> Foam::Pstream::listGatherValues
(
    const T& localValue,
    const int communicator,
    [[maybe_unused]] const int tag
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(communicator) as well?
        List<T> allValues(1);
        allValues[0] = localValue;
        return allValues;
    }
    else if constexpr (is_contiguous_v<T>)
    {
        // UPstream version is contiguous only
        return UPstream::listGatherValues(localValue, communicator);
    }
    else
    {
        // Non-contiguous: using all-to-one/one-to-all communication

        // Non-trivial to manage non-blocking gather without a
        // PEX/NBX approach (eg, PstreamBuffers).
        // Leave with simple exchange for now

        List<T> allValues;
        if (UPstream::master(communicator))
        {
            allValues.resize(UPstream::nProcs(communicator));

            for (auto proci : UPstream::subProcs(communicator))
            {
                IPstream::recv(allValues[proci], proci, tag, communicator);
            }

            allValues[0] = localValue;
        }
        else if (UPstream::is_rank(communicator))
        {
            OPstream::send(localValue, UPstream::masterNo(), tag, communicator);
        }

        return allValues;
    }
}


template<class T>
T Foam::Pstream::listScatterValues
(
    const UList<T>& allValues,
    const int communicator,
    [[maybe_unused]] const int tag
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return first value
        // TBD: only when UPstream::is_rank(communicator) as well?

        if (!allValues.empty())
        {
            return allValues[0];
        }

        return T{};  // Fallback value
    }
    else if constexpr (is_contiguous_v<T>)
    {
        // UPstream version is contiguous only
        return UPstream::listScatterValues(allValues, communicator);
    }
    else
    {
        // Non-contiguous: using all-to-one/one-to-all communication

        T localValue{};

        if (UPstream::master(communicator))
        {
            const auto numProc = UPstream::nProcs(communicator);

            if (FOAM_UNLIKELY(allValues.size() < numProc))
            {
                FatalErrorInFunction
                    << "Attempting to send " << allValues.size()
                    << " values to " << numProc << " processors" << endl
                    << Foam::abort(FatalError);
            }

            const auto startOfRequests = UPstream::nRequests();

            List<DynamicList<char>> sendBuffers(numProc);

            for (auto proci : UPstream::subProcs(communicator))
            {
                UOPstream toProc
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendBuffers[proci],
                    tag,
                    communicator
                );
                toProc << allValues[proci];
            }

            // Wait for outstanding requests
            UPstream::waitRequests(startOfRequests);

            return allValues[0];
        }
        else if (UPstream::is_rank(communicator))
        {
            IPstream::recv(localValue, UPstream::masterNo(), tag, communicator);
        }

        return localValue;
    }
}


// ************************************************************************* //
