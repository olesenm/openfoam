/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2025 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Request::Request() noexcept
:
    UPstream::Request(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Request::good() const noexcept
{
    return false;
}


void Foam::UPstream::Request::reset() noexcept
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// Foam::UPstream::Request
// Foam::UPstream::Request::lookup(int req)
// {
//     return UPstream::Request(nullptr);
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int Foam::UPstream::nRequests() noexcept { return 0; }

void Foam::UPstream::resetRequests(int n) {}

void Foam::UPstream::addRequest(UPstream::Request&) {}

void Foam::UPstream::cancelRequest(int i) {}
void Foam::UPstream::cancelRequest(UPstream::Request&) {}
void Foam::UPstream::cancelRequests(UList<UPstream::Request>&) {}

void Foam::UPstream::removeRequests(int pos, int len) {}

void Foam::UPstream::freeRequest(UPstream::Request&) {}
void Foam::UPstream::freeRequests(UList<UPstream::Request>&) {}

void Foam::UPstream::waitRequests(int pos, int len) {}
void Foam::UPstream::waitRequests(UList<UPstream::Request>&) {}

bool Foam::UPstream::waitAnyRequest(int pos, int len) { return false; }

bool Foam::UPstream::waitSomeRequests
(
    int pos,
    int len,
    DynamicList<int>* indices
)
{
    if (indices) indices->clear();
    return false;
}

bool Foam::UPstream::waitSomeRequests
(
    UList<UPstream::Request>& requests,
    DynamicList<int>* indices
)
{
    if (indices) indices->clear();
    return false;
}

int Foam::UPstream::waitAnyRequest(UList<UPstream::Request>&)
{
    return -1;
}

void Foam::UPstream::waitRequest(int i) {}
void Foam::UPstream::waitRequest(UPstream::Request&) {}

bool Foam::UPstream::activeRequest(int i) { return false; }
bool Foam::UPstream::activeRequest(const UPstream::Request&) { return false; }

bool Foam::UPstream::finishedRequest(int i) { return true; }
bool Foam::UPstream::finishedRequest(UPstream::Request&) { return true; }

bool Foam::UPstream::finishedRequests(int pos, int len) { return true; }


bool Foam::UPstream::finishedRequests(UList<UPstream::Request>&)
{
    return true;
}


bool Foam::UPstream::finishedRequestPair(int& req0, int& req1)
{
    req0 = -1;
    req1 = -1;
    return true;
}


void Foam::UPstream::waitRequestPair(int& req0, int& req1)
{
    req0 = -1;
    req1 = -1;
}


// ************************************************************************* //
