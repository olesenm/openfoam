/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

Application
    Test-surfaceReading

Description
    Test basic surface format reading capabilities (and speeds)

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "clockTime.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "MeshedSurfaces.H"
#include "UnsortedMeshedSurfaces.H"
#include "STLReader.H"

using namespace Foam;


template<class T>
Ostream& print_average(Ostream& os, const UList<T>& values)
{
    if (values.empty()) os << "n/a"; else os << average(values);
    return os;
}


template<class T>
Ostream& print_minmax(Ostream& os, const UList<T>& values)
{
    if (values.empty()) os << "n/a"; else os << minMax(values);
    return os;
}

void test_measure(const fileName& importName)
{
    // Dummy Time, without needing a case directory
    autoPtr<Time> timePtr(Time::NewGlobalTime());

    fileName resolvedName(importName);
    resolvedName.toAbsolute();

    IOobject io
    (
        resolvedName,   // absolute path
       *timePtr,
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    triSurfaceMesh surf(io);

    List<vector> centres;
    List<scalar> radiusSqr;
    List<label> regions;
    List<vector> normals;

    // const label n = surf.measure(nullptr, nullptr, nullptr, nullptr);
    // const label n = surf.measure(&centres, &radiusSqr, &regions, &normals);
    // const label n = surf.measure(nullptr, &radiusSqr, nullptr, nullptr);
    const label n = surf.measure(&centres, nullptr, &regions, &normals);
    // const label n = surf.measure(nullptr, nullptr, &regions, nullptr);

    Info<< "Read surface:" << endl;
    surf.writeStats(Info);
    Info<< "measure() returned " << n << " elements" << nl;

    Info<< "  centres  : " << centres.size()
        << "  average = "; print_average(Info, centres) << nl;

    Info<< "  radiusSqr: " << radiusSqr.size()
        << "  min/max = "; print_minmax(Info, radiusSqr) << nl;

    Info<< "  regions  : " << regions.size()
        << "  min/max = "; print_minmax(Info, regions) << nl;

    Info<< "  normals  : " << normals.size()
        << "  average = "; print_average(Info, normals) << nl;

    Info<< endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test basic surface format reading capabilities (and speeds)"
    );

    argList::noParallel();
    argList::addArgument("inputFile");

    argList::addBoolOption
    (
        "triSurface",
        "Use triSurface for read"
    );
    argList::addBoolOption
    (
        "triFace",
        "Use triFace instead of face"
    );
    argList::addBoolOption
    (
        "unsorted",
        "Use UnsortedMeshedSurface instead of MeshedSurface, "
        "or unsorted output (with -triSurface option)"
    );

    argList::addBoolOption
    (
        "measure",
        "Load as triSurfaceMesh (with dummy Time) and report measure()"
    );

    argList::addOption
    (
        "ext",
        "name",
        "Force alternative extension"
    );

    argList::addOption
    (
        "stl-parser",
        "N",
        "ASCII parser type: 0=Flex, 1=Ragel, 2=Manual"
    );

    #include "setRootCase.H"

    const auto importName = args.get<fileName>(1);

    word ext =
    (
        importName.has_ext("gz")
      ? importName.stem().ext()
      : importName.ext()
    );

    // Allow override of extension
    args.readIfPresent("ext", ext);

    args.readIfPresent("stl-parser", fileFormats::STLReader::parserType);

    clockTime timing;

    if (args.found("measure"))
    {
        test_measure(importName);
    }
    else if (args.found("triSurface"))
    {
        triSurface surf(importName, ext);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area        : " << sum(surf.magSf()) << nl << endl;
    }
    else if (args.found("triFace"))
    {
        MeshedSurface<triFace> surf(importName, ext);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area        : " << sum(surf.magSf()) << nl << endl;
    }
    else if (args.found("unsorted"))
    {
        UnsortedMeshedSurface<face> surf(importName, ext);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area        : " << sum(surf.magSf()) << nl << endl;
    }
    else
    {
        MeshedSurface<face> surf(importName, ext);

        Info<< "Read surface:" << endl;
        surf.writeStats(Info);
        Info<< "Area        : " << sum(surf.magSf()) << nl << endl;
    }

    Info<< nl << "Reading took " << timing.elapsedTime() << "s" << nl
        << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
