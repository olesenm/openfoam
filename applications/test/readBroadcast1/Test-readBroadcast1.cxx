/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

Description
    Test file reading with broadcast

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Fstream.H"
#include "Pstream.H"
#include "SpanStream.H"

using namespace Foam;

// Common implementation: readFn is called with IFstream& or fileName
template<class BufferType, class ReadFunc>
int64_t do_reading
(
    BufferType& buffer,
    const fileName& srcName,
    const bool syncPar,
    const bool optFail,
    ReadFunc readFn
)
{
    const int communicator = UPstream::worldComm;

    if (UPstream::master(communicator) || !syncPar)
    {
        if (optFail)
        {
            IFstream ifs(srcName, IOstreamOption::BINARY);

            if (!ifs.good())
            {
                FatalIOErrorInFunction(srcName)
                    << "Cannot open file " << srcName
                    << exit(FatalIOError);
            }

            buffer = readFn(ifs);
        }
        else
        {
            buffer = readFn(srcName);
        }
    }

    int64_t len(buffer.size());

    if (syncPar && UPstream::parRun())
    {
        Pstream::broadcast(len);

        if (UPstream::is_subrank(communicator))
        {
            buffer.resize(len);
        }

        if (len)
        {
            // Only broadcast non-empty content
            UPstream::broadcast(buffer.data(), buffer.size(), communicator);
        }
    }

    return len;
}


int64_t do_reading
(
    std::string& buffer,
    const fileName& srcName,
    const bool syncPar,
    const bool optFail
)
{
    return do_reading
    (
        buffer, srcName, syncPar, optFail,
        [](auto& arg) { return IFstream::readContents_string(arg); }
    );
}


int64_t do_reading
(
    std::vector<char>& buffer,
    const fileName& srcName,
    const bool syncPar,
    const bool optFail
)
{
    return do_reading
    (
        buffer, srcName, syncPar, optFail,
        [](auto& arg) { return IFstream::readContents_vector(arg); }
    );
}


int64_t do_reading
(
    DynamicList<char>& buffer,
    const fileName& srcName,
    const bool syncPar,
    const bool optFail
)
{
    return do_reading
    (
        buffer, srcName, syncPar, optFail,
        [](auto& arg) { return IFstream::readContents_list(arg); }
    );
}


void processLines(std::string_view sv)
{
    ISpanStream is(sv.data(), sv.size());

    Pout<< "input:" << is.capacity() << endl;

    for (string line; is.getLine(line); /*nil*/)
    {
        Pout<< "L:" << is.lineNumber() << ": " << line.c_str() << nl;
    }
}

void processLines(const char* data, size_t len)
{
    std::string_view sv(data, len);
    processLines(sv);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noFunctionObjects();
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("additional information");
    argList::addBoolOption("fail", "fail if file cannot be opened");
    argList::addBoolOption("no-broadcast", "suppress broadcast contents");

    argList::addBoolOption("string", "use std::string");
    argList::addBoolOption("vector", "use std::vector<char>");
    argList::addBoolOption("dynlist", "use DynamicList<char>");

    argList::addNote("Test master-only reading (with broadcast)");

    argList::addArgument("srcFile");

    #include "setRootCase.H"

    const bool syncPar = (UPstream::parRun() && !args.found("no-broadcast"));
    const bool optFail = args.found("fail");

    auto srcName = args.get<fileName>(1);

    if (srcName.has_ext("gz"))
    {
        srcName.remove_ext();
        Info<< "stripping extraneous .gz ending" << endl;
    }

    bool needsTest = true;

    if (args.found("string"))
    {
        std::string buffer;
        needsTest = false;

        Info<< "Using std::string reading" << endl;
        do_reading(buffer, srcName, syncPar, optFail);

        Pout<< "content-size: " << int64_t(buffer.size()) << nl;
        processLines(buffer.data(), buffer.size());
    }

    if (args.found("vector"))
    {
        std::vector<char> buffer;

        needsTest = false;
        Info<< "Using std::vector<char> reading" << endl;
        do_reading(buffer, srcName, syncPar, optFail);

        Pout<< "content-size: " << int64_t(buffer.size()) << nl;
        processLines(buffer.data(), buffer.size());
    }

    if (args.found("dynlist") || needsTest)
    {
        DynamicList<char> buffer;
        needsTest = false;

        Info<< "Using DynamicList<char> reading" << endl;
        do_reading(buffer, srcName, syncPar, optFail);

        Pout<< "content-size: " << int64_t(buffer.size()) << nl;
        processLines(buffer.data(), buffer.size());
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
