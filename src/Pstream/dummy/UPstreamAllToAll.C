/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
    Copyright (C) 2026 Mark Olesen
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
#include "Map.H"

#include <cinttypes>
#include <cstring>  // memmove

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Type)                                          \
void Foam::UPstream::allToAll_n                                               \
(                                                                             \
    const UList<Type>& sendData,                                              \
    UList<Type>& recvData,                                                    \
    const int count,                                                          \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    if (count < 1) recvData = Foam::zero{};                                   \
    else recvData.deepCopy(sendData);                                         \
}

Pstream_CommonRoutines(int32_t);
Pstream_CommonRoutines(int64_t);

#undef Pstream_CommonRoutines


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Type)                                          \
void Foam::UPstream::allToAllConsensus_n                                      \
(                                                                             \
    const UList<Type>& sendData,                                              \
    UList<Type>& recvData,                                                    \
    const int count,                                                          \
    const int tag,                                                            \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    if (count < 1) recvData = Foam::zero{};                                   \
    else recvData.deepCopy(sendData);                                         \
}                                                                             \
                                                                              \
void Foam::UPstream::allToAllConsensus                                        \
(                                                                             \
    const Map<Type>& sendData,                                                \
    Map<Type>& recvData,                                                      \
    const int tag,                                                            \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    recvData = sendData;                                                      \
}


Pstream_CommonRoutines(int32_t);
Pstream_CommonRoutines(int64_t);

#undef Pstream_CommonRoutines


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Type)                                          \
void Foam::UPstream::allToAllv                                                \
(                                                                             \
    const Type* sendData,                                                     \
    const UList<int>& sendCounts,                                             \
    const UList<int>& sendOffsets,                                            \
    Type* recvData,                                                           \
    const UList<int>& recvCounts,                                             \
    const UList<int>& recvOffsets,                                            \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    if (recvCounts[0] != sendCounts[0])                                       \
    {                                                                         \
        FatalErrorInFunction                                                  \
            << "Number to send " << sendCounts[0]                             \
            << " does not equal number to receive " << recvCounts[0]          \
            << Foam::abort(FatalError);                                       \
    }                                                                         \
    std::memmove(recvData, sendData, recvCounts[0]*sizeof(Type));             \
}


// Unused: Pstream_CommonRoutines(char);

#undef Pstream_CommonRoutines

// ************************************************************************* //
