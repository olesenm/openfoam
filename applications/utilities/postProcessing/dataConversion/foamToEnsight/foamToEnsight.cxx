/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2026 OpenCFD Ltd.
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
    foamToEnsight

Group
    grpPostProcessingUtilitie

Description
    Translate OpenFOAM data to EnSight format.
    An Ensight part is created for cellZones (unzoned cells are "internalMesh")
    and patches.

    - Handles volume fields, dimensioned fields, point fields
    - Handles mesh topology changes.

Usage
    \b foamToEnsight [OPTION]

    Options:
      - \par -ascii
        Write Ensight data in ASCII format instead of "C Binary"

      - \par -fields \<fields\>
        Specify single or multiple fields to write (all by default)
        For example,
        \verbatim
          -fields T
          -fields '(p T U \"alpha.*\")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -nearCellValue
        Use zero-gradient cell values on patches

      - \par -nodeValues
        Force interpolation of values to nodes

      - \par -no-boundary
        Suppress output for all boundary patches

      - \par -no-internal
        Suppress output for internal (volume) mesh

      - \par -no-cellZones
        Suppress cellZone handling

      - \par -no-lagrangian
        Suppress writing lagrangian positions and fields.

      - \par -no-mesh
        Suppress writing the geometry. Can be useful for converting partial
        results for a static geometry.

      - \par -no-point-data
        Suppress conversion of pointFields. No interpolated PointData.

      - \par -noZero
        Exclude the often incomplete initial conditions.

      - \par -index \<start\>
        Use consecutive indexing for \c data/######## files with the
        specified start index.
        Ignore the time index contained in the uniform/time file.

      - \par -name \<subdir\>
        Define sub-directory name to use for Ensight data (default: "EnSight")

      - \par -width \<n\>
        Width of Ensight data subdir (default: 8)

      - \par -cellZones NAME | LIST
        Specify single zone or multiple cell zones (name or regex) to write

      - \par -faceZones NAME | LIST
        Specify single zone or multiple face zones (name or regex) to write

      - \par -patches NAME | LIST
        Specify single patch or multiple patches (name or regex) to write
        For example,
        \verbatim
          -patches top
          -patches '( front \".*back\" )'
        \endverbatim

      - \par -exclude-patches NAME | LIST
        Exclude single or multiple patches (name or regex) from writing.
        For example,
        \verbatim
          -exclude-patches '( inlet_1 inlet_2 "proc.*" )'
        \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "Pstream.H"
#include "HashOps.H"
#include "PtrDynList.H"
#include "regionProperties.H"

#include "fvc.H"
#include "fvMesh.H"
#include "fieldTypes.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "vectorIOField.H"

// file-format/conversion
#include "ensightCase.H"
#include "ensightGeoFile.H"
#include "ensightFaMesh.H"
#include "ensightMesh.H"
#include "ensightOutputCloud.H"
#include "ensightOutputAreaField.H"
#include "ensightOutputVolField.H"

// local files
#include "readFields.H"
#include "writeVolFields.H"
#include "writeDimFields.H"
#include "writePointFields.H"
#include "writeAreaFields.H"

#include "memInfo.H"

#undef foamToEnsight_useTimeIndex

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Translate OpenFOAM data to Ensight format with individual parts"
        " for cellZones, unzoned cells and patches"
    );
    timeSelector::addOptions();

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");

    argList::addVerboseOption();

    #include "addAllRegionOptions.H"
    #include "addAllFaRegionOptions.H"

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of 'C Binary'"
    );
    argList::addOption
    (
        "index",
        "start",
        "Starting index for consecutive number of Ensight data/ files."
        " Ignore the time index contained in the uniform/time file."
        , true  // mark as an advanced option
    );
    argList::addOption
    (
        "name",
        "dir",
        "Directory name for output (default: 'EnSight'),"
        " relative to case dir, or an absolute path."
    );
    argList::addOption
    (
        "base-name",
        "name",
        "The base/stem for output files (default: <case>)"
    );
    argList::addBoolOption
    (
        "no-overwrite",
        "Suppress removal of existing EnSight output directory"
    );
    argList::addOption
    (
        "width",
        "n",
        "Width of Ensight data subdir"
    );
    argList::addBoolOption
    (
        "nearCellValue",
        "Use zero-gradient cell values on patches"
        , true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "nodeValues",
        "Force interpolation of values to nodes"
        , true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "no-boundary",  // noPatches
        "Suppress writing any patches"
    );
    argList::addOptionCompat("no-boundary", {"noPatches", 1806});

    argList::addBoolOption
    (
        "no-internal",
        "Suppress writing the internal mesh"
    );
    argList::addBoolOption
    (
        "no-cellZones",
        "Suppress writing any cellZones"
    );
    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Suppress writing lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 1806});

    argList::addBoolOption
    (
        "no-point-data",
        "Suppress conversion of pointFields, disable -nodeValues"
    );
    argList::addBoolOption
    (
        "no-mesh", // noMesh
        "Suppress writing the geometry."
        " Can be useful for converting partial results for a static geometry"
        , true  // mark as an advanced option
    );

    // Future?
    // argList::addBoolOption
    // (
    //     "one-boundary",  // allPatches
    //     "Combine all patches into a single part"
    // );
    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress output of finite-area mesh/fields",
        true  // mark as an advanced option
    );
    argList::ignoreOptionCompat
    (
        {"finite-area", 2112},  // use -no-finite-area to disable
        false           // bool option, no argument
    );

    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to write\n"
        "Eg, 'inlet' or '(outlet \"inlet.*\")'"
    );
    argList::addOption
    (
        "exclude-patches",
        "wordRes",
        "Exclude single or multiple patches from writing\n"
        "Eg, 'outlet' or '( inlet \".*Wall\" )'"
        , true  // mark as an advanced option
    );
    argList::addOptionCompat("exclude-patches", {"excludePatches", 2112});

    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple face zones to write\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'."
    );
    argList::addOption
    (
        "exclude-faceZones",
        "wordRes",
        "Exclude single or multiple face zones from writing"
        , true  // mark as an advanced option
    );

    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to write (all by default)\n"
        "Eg, 'T' or '( \"U.*\" )'"
    );
    argList::addOption
    (
        "exclude-fields",
        "wordRes",
        "Exclude single or multiple fields",
        true  // mark as an advanced option
    );
    argList::addBoolOption
    (
        "no-fields",
        "Suppress conversion of fields"
    );

    argList::addOption
    (
        "cellZones",
        "wordRes",
        "Specify single or multiple cell zones to write\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'."
    );
    argList::addOptionCompat("cellZones", {"cellZone", 1912});
    argList::addOption
    (
        "exclude-cellZones",
        "wordRes",
        "Exclude single or multiple cell zones from writing"
        , true  // mark as an advanced option
    );

    // Prevent volume BCs from triggering finite-area
    regionModels::allowFaModels(false);

    #include "setRootCase.H"

    // ------------------------------------------------------------------------
    // Configuration

    // Default to binary output, unless otherwise specified
    const IOstreamOption::streamFormat format =
    (
        args.found("ascii")
      ? IOstreamOption::ASCII
      : IOstreamOption::BINARY
    );

    const int optVerbose = args.verbose();
    const bool doBoundary    = !args.found("no-boundary");
    const bool doInternal    = !args.found("no-internal");
    const bool doCellZones   = !args.found("no-cellZones");
    const bool doLagrangian  = !args.found("no-lagrangian");
    const bool doFiniteArea  = !args.found("no-finite-area");
    const bool doPointValues = !args.found("no-point-data");
    const bool nearCellValue = args.found("nearCellValue") && doBoundary;

    // Control for numbering iterations
    label indexingNumber(0);
    const bool doConsecutive = args.readIfPresent("index", indexingNumber);

    // Write the geometry, unless otherwise specified
    bool doGeometry = !args.found("no-mesh");

    if (nearCellValue)
    {
        Info<< "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }
    if (!doPointValues)
    {
        Info<< "Point fields and interpolated point data"
            << " disabled with the '-no-point-data' option"
            << nl;
    }

    //
    // General (case) output options
    //
    ensightCase::options caseOpts(format);

    // Forced point interpolation?
    caseOpts.nodeValues(doPointValues && args.found("nodeValues"));
    caseOpts.width(args.getOrDefault<label>("width", 8));
    caseOpts.overwrite(!args.found("no-overwrite")); // Remove existing?

    // Can also have separate directory for lagrangian
    // caseOpts.separateCloud(true);

    ensightMesh::options writeOpts;
    writeOpts.useBoundaryMesh(doBoundary);
    writeOpts.useInternalMesh(doInternal);
    writeOpts.useCellZones(doCellZones);

    // Include/exclude patches
    if (wordRes list; args.readListIfPresent<wordRe>("patches", list))
    {
        writeOpts.patchInclude(list);
    }
    if (wordRes list; args.readListIfPresent<wordRe>("exclude-patches", list))
    {
        writeOpts.patchExclude(list);
    }

    // Include/exclude cell zones
    if (wordRes list; args.readListIfPresent<wordRe>("cellZones", list))
    {
        writeOpts.cellZoneInclude(list);
    }
    if (wordRes list; args.readListIfPresent<wordRe>("exclude-cellZones", list))
    {
        writeOpts.cellZoneExclude(list);
    }

    // Include/exclude face zones
    if (wordRes list; args.readListIfPresent<wordRe>("faceZones", list))
    {
        writeOpts.faceZoneInclude(list);
    }
    if (wordRes list; args.readListIfPresent<wordRe>("exclude-faceZones", list))
    {
        writeOpts.faceZoneExclude(list);
    }

    // Report the setup
    writeOpts.print(Info);

    // Field selection/deselection
    wordRes includedFields, excludedFields;
    const bool doConvertFields = !args.found("no-fields");
    if (doConvertFields)
    {
        if (args.readListIfPresent<wordRe>("fields", includedFields))
        {
            Info<< "Including fields "
                << flatOutput(includedFields) << nl << endl;
        }
        if (args.readListIfPresent<wordRe>("exclude-fields", excludedFields))
        {
            Info<< "Excluding fields "
                << flatOutput(excludedFields) << nl << endl;
        }
    }
    else
    {
        Info<< "Field conversion disabled with the '-no-fields' option" << nl;
    }

    const wordRes::filter fieldSelector(includedFields, excludedFields);

    // ------------------------------------------------------------------------

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    // Handle volume region selections
    #include "getAllRegionOptions.H"

    // Handle area region selections
    #include "getAllFaRegionOptions.H"

    if (!doFiniteArea)
    {
        areaRegionNames.clear();  // For consistency
    }

    // ------------------------------------------------------------------------
    // Directory management

    // Define directory name to use for output data.
    // The output path is at case level (or global path) only.
    // - For parallel cases, data only written from master

    fileName outputDir(args.globalPath()/"EnSight");
    if (fileName dir; args.readIfPresent("name", dir) && !dir.empty())
    {
        if (dir.isAbsolute())
        {
            outputDir = std::move(dir);
        }
        else
        {
            outputDir = args.globalPath()/dir;
        }
        outputDir.clean();  // Remove unneeded ".."
    }

    // The base name for naming output files - without directory!
    word outputBaseName;
    if (fileName fn; args.readIfPresent("base-name", fn) && !fn.empty())
    {
        fn.expand();
        outputBaseName = fn.name();
    }
    else
    {
        outputBaseName = args.globalCaseName().name();
    }

    // ------------------------------------------------------------------------
    cpuTime timer;
    Info<< "Initial memory " << Foam::memInfo{}.size() << " kB" << endl;

    #include "createNamedMeshes.H"
    #include "createMeshAccounting.H"

    if (UPstream::master())
    {
        Info<< "Converting " << timeDirs.size() << " time steps" << nl;
        // ensCase.printInfo(Info) << endl;
    }

    // Check mesh motion
    #include "checkMeshMoving.H"
    if (hasMovingMesh && !doGeometry)
    {
        Info<< "has moving mesh: ignoring '-no-mesh' option" << endl;
        doGeometry = true;
    }

    // Check lagrangian
    #include "findCloudFields.H"

    // Check field availability
    #include "checkFieldAvailability.H"

    // test the pre-check variable if there is a moving mesh
    // time-set for geometries
    // TODO: split off into separate time-set,
    // but need to verify ensight spec

    Info<< "Startup in "
        << timer.cpuTimeIncrement() << " s, "
        << Foam::memInfo{}.size() << " kB" << nl << endl;


    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        // Index for the Ensight case(s). Continues if not possible
        #include "getTimeIndex.H"

        Info<< "Time [" << timeIndex << "] = " << runTime.timeName() << nl;

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word& regionDir = polyMesh::regionName(regionName);

            if (regionNames.size() > 1)
            {
                Info<< "region=" << regionName << nl;
            }

            auto& mesh = meshes[regioni];

            polyMesh::readUpdateState meshState = mesh.readUpdate();
            const bool moving = (meshState != polyMesh::UNCHANGED);

            // Ensight (fvMesh)
            auto& ensCase = ensightCases[regioni];
            auto& ensMesh = ensightMeshes[regioni];

            ensCase.setTime(timeDirs[timei], timeIndex);

            // Finite-area (optional)
            // Accounting exists for each volume region but may be empty
            auto& ensFaCases = ensightCasesFa[regioni];
            auto& ensFaMeshes = ensightMeshesFa[regioni];

            for (auto& ensFaCase : ensFaCases)
            {
                ensFaCase.setTime(timeDirs[timei], timeIndex);
            }

            // Movement
            if (moving)
            {
                ensMesh.expire();
                ensMesh.correct();

                for (auto& ensFaMesh : ensFaMeshes)
                {
                    ensFaMesh.expire();
                    ensFaMesh.correct();
                }
            }

            if ((timei == 0 || moving) && doGeometry)
            {
                // finite-volume
                {
                    autoPtr<ensightGeoFile> os =
                        ensCase.newGeometry(hasMovingMesh);

                    ensMesh.write(os.ref());
                }

                // finite-area
                forAll(ensFaMeshes, areai)
                {
                    const auto& ensFaCase = ensFaCases[areai];
                    const auto& ensFaMesh = ensFaMeshes[areai];

                    autoPtr<ensightGeoFile> os =
                        ensFaCase.newGeometry(hasMovingMesh);

                    ensFaMesh.write(os.ref());
                }
            }


            // Objects at this time
            IOobjectList objects(mesh, runTime.timeName());

            objects.filterObjects
            (
                availableRegionObjectNames[regioni]
            );

            // Volume, internal, point fields
            #include "convertVolumeFields.H"

            // finite-area
            forAll(ensFaMeshes, areai)
            {
                auto* ensFaCasePtr = ensFaCases.get(areai);
                auto* ensFaMeshPtr = ensFaMeshes.get(areai);

                // The finite-area region objects at this time
                IOobjectList faObjects;

                if (ensFaMeshPtr)  // Cannot really fail
                {
                    const auto& areaMesh = ensFaMeshPtr->mesh();
                    const word& areaName = areaMesh.name();

                    faObjects = IOobjectList
                    (
                        areaMesh,
                        runTime.timeName(),
                        IOobjectOption::NO_REGISTER
                    );

                    faObjects.filterObjects
                    (
                        availableFaRegionObjectNames[regioni](areaName)
                    );
                }

                // The finiteArea fields
                #include "convertAreaFields.H"
            }

            // Lagrangian fields
            #include "convertLagrangian.H"
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << Foam::memInfo{}.size() << " kB" << nl << nl;
    }

    // Write cases
    forAll(ensightCases, regioni)
    {
        // finite-volume
        ensightCases[regioni].write();

        // Finite-area (if any)
        for (const auto& ensFaCase : ensightCasesFa[regioni])
        {
            ensFaCase.write();
        }
    }

    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << Foam::memInfo{}.peak() << " kB (peak)" << nl << endl;

    return 0;
}


// ************************************************************************* //
