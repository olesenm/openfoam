/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    decomposePar

Group
    grpParallelUtilities

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of OpenFOAM.

Usage
    \b decomposePar [OPTIONS]

    Options:
      - \par -allRegions
        Decompose all regions in regionProperties. Does not check for
        existence of processor*.

      - \par -case \<dir\>
        Specify case directory to use (instead of the cwd).

      - \par -cellDist
        Write the cell distribution as a labelList, for use with 'manual'
        decomposition method and as a VTK or volScalarField for visualization.

      - \par -constant
        Include the 'constant/' dir in the times list.

      - \par -copyUniform
        Copy any \a uniform directories too.

      - \par -copyZero
        Copy \a 0 directory to processor* rather than decompose the fields.

      - \par -debug-switch \<name=val\>
        Specify the value of a registered debug switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -decomposeParDict \<file\>
        Use specified file for decomposePar dictionary.

      - \par -dry-run
        Test without writing the decomposition. Changes -cellDist to
        only write VTK output.

      - \par -fields
        Use existing geometry decomposition and convert fields only.

      - \par fileHandler \<handler\>
        Override the file handler type.

      - \par -force
        Remove any existing \a processor subdirectories before decomposing the
        geometry.

      - \par -ifRequired
        Only decompose the geometry if the number of domains has changed from a
        previous decomposition. No \a processor subdirectories will be removed
        unless the \a -force option is also specified. This option can be used
        to avoid redundant geometry decomposition (eg, in scripts), but should
        be used with caution when the underlying (serial) geometry or the
        decomposition method etc. have been changed between decompositions.

      - \par -info-switch \<name=val\>
        Specify the value of a registered info switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -latestTime
        Select the latest time.

      - \par -lib \<name\>
        Additional library or library list to load (can be used multiple times).

      - \par -noFunctionObjects
        Do not execute function objects.

      - \par -noSets
        Skip decomposing cellSets, faceSets, pointSets.

      - \par -noZero
        Exclude the \a 0 dir from the times list.

      - \par -opt-switch \<name=val\>
        Specify the value of a registered optimisation switch (int/bool).
        Default is 1 if the value is omitted. (Can be used multiple times)

      - \par -region \<regionName\>
        Decompose named region. Does not check for existence of processor*.

      - \par -time \<ranges\>
        Override controlDict settings and decompose selected times. Does not
        re-decompose the mesh i.e. does not handle moving mesh or changing
        mesh cases. Eg, ':10,20 40:70 1000:', 'none', etc.

      - \par -verbose
        Additional verbosity.

      - \par -doc
        Display documentation in browser.

      - \par -doc-source
        Display source code in browser.

      - \par -help
        Display short help and exit.

      - \par -help-man
        Display full help (manpage format) and exit.

      - \par -help-notes
        Display help notes (description) and exit.

      - \par -help-full
        Display full help and exit.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "OSspecific.H"
#include "IOobjectList.H"

#include "decompositionModel.H"
#include "domainDecomposition.H"
#include "domainDecompositionDryRun.H"

#include "regionProperties.H"

#include "fieldsDistributor.H"

#include "fvFieldDecomposer.H"
#include "pointFields.H"
#include "pointFieldDecomposer.H"

#include "lagrangianFieldDecomposer.H"

#include "emptyFaPatch.H"
#include "faFieldDecomposer.H"
#include "faMeshDecomposition.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Return cached finite-volume proc addressing,
// or read from facesInstance
FOAM_NO_DANGLING_REFERENCE
const labelIOList& procAddressing
(
    const UPtrList<fvMesh>& procMeshList,
    const label proci,
    const word& name,
    PtrList<labelIOList>& procAddressingList
)
{
    const auto& procMesh = procMeshList[proci];

    // Allow lazy initial sizing
    if (procAddressingList.size() < procMeshList.size())
    {
        procAddressingList.resize(procMeshList.size());
    }

    return procAddressingList.try_emplace
    (
        proci,
        IOobject
        (
            name,
            procMesh.facesInstance(),
            polyMesh::meshSubDir,  // local
            procMesh,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        )
    );
}


// Return cached finite-area proc addressing,
// or read from facesInstance (which is normally just "constant")
FOAM_NO_DANGLING_REFERENCE
const labelIOList& procAddressing
(
    const UPtrList<faMesh>& procMeshList,
    const label proci,
    const word& name,
    PtrList<labelIOList>& procAddressingList
)
{
    const auto& procMesh = procMeshList[proci];

    // Allow lazy initial sizing
    if (procAddressingList.size() < procMeshList.size())
    {
        procAddressingList.resize(procMeshList.size());
    }

    return procAddressingList.try_emplace
    (
        proci,
        IOobject
        (
            name,
            procMesh.facesInstance(),
            faMesh::meshSubDir,  // local
            procMesh,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        )
    );
}


// Return cached processor Time or create
Foam::Time& emplaceTime
(
    PtrList<Time>& procTimes,
    const label proci,
    const argList& args
)
{
    if (!procTimes.test(proci))
    {
        procTimes.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/("processor" + Foam::name(proci)),
                args.allowFunctionObjects(),
                args.allowLibs()
            )
        );
    }
    return procTimes[proci];
}


void decomposeUniform
(
    const bool copyUniform,
    const domainDecomposition& mesh,
    const Time& processorDb,
    const word& regionDir = word::null
)
{
    const Time& runTime = mesh.time();

    // Any uniform data to copy/link?
    const fileName uniformDir(regionDir/"uniform");

    if (fileHandler().isDir(runTime.timePath()/uniformDir))
    {
        Info<< "Detected additional non-decomposed files in "
            << runTime.relativePath(uniformDir) << endl;

        // Bit of trickery to synthesise the correct directory base,
        // e.g. processors4/0.01
        const fileName timePath = fileHandler().objectPath
        (
            IOobject
            (
                "dummy",
                runTime.timeName(),
                processorDb
            ),
            word::null
        ).path();

        // If no fields have been decomposed the destination
        // directory will not have been created so make sure.
        Foam::mkDir(timePath);

        if (copyUniform || mesh.distributed())
        {
            if (!fileHandler().exists(timePath/uniformDir))
            {
                fileHandler().cp
                (
                    runTime.timePath()/uniformDir,
                    timePath/uniformDir
                );
            }
        }
        else
        {
            // Link with relative paths
            string parentPath = string("..")/"..";

            if (!regionDir.empty())
            {
                parentPath = parentPath/"..";
            }

            fileName currentDir(cwd());
            Foam::chDir(timePath);

            if (!fileHandler().exists(uniformDir))
            {
                fileHandler().ln
                (
                    parentPath/runTime.timeName()/uniformDir,
                    uniformDir
                );
            }
            Foam::chDir(currentDir);
        }
    }
}

} // End namespace Foam


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Decompose a mesh and fields of a case for parallel execution"
    );

    argList::noParallel();
    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Alternative decomposePar dictionary file"
    );

    #include "addAllRegionOptions.H"
    #include "addAllFaRegionOptions.H"

    argList::addVerboseOption
    (
        "Additional verbosity. (Can be used multiple times for debugging)"
    );
    argList::addDryRunOption
    (
        "Test without writing the decomposition. "
        "Changes -cellDist to only write VTK output."
    );
    argList::addOption
    (
        "domains",
        "N",
        "Override numberOfSubdomains (-dry-run only)",
        true  // Advanced option
    );
    argList::addOption
    (
        "method",
        "name",
        "Override decomposition method (-dry-run only)",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress finiteArea mesh/field decomposition",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "no-lagrangian",
        "Suppress lagrangian (cloud) decomposition",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "disable-edge-encoding",
        "Emit edgeProcAddressing without encoding edge flips, "
        "as per 2512 and earlier [special use]",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as a labelList - for use with 'manual' "
        "decomposition method and as a volScalarField for visualization."
    );
    argList::addBoolOption
    (
        "copyZero",
        "Copy 0/ directory to processor*/ rather than decompose the fields"
    );
    argList::addBoolOption
    (
        "copyUniform",
        "Copy any uniform/ directories too"
    );
    argList::addBoolOption
    (
        "fields",
        "Use existing geometry decomposition and convert fields only"
    );
    argList::addBoolOption
    (
        "no-fields",
        "Suppress conversion of fields (volume, finite-area, lagrangian)"
    );

    argList::addBoolOption
    (
        "no-sets",
        "Skip decomposing cellSets, faceSets, pointSets"
    );
    argList::addOptionCompat("no-sets", {"noSets", 2106});

    argList::addBoolOption
    (
        "force",
        "Remove existing processor*/ subdirs before decomposing the geometry"
    );
    argList::addBoolOption
    (
        "ifRequired",
        "Only decompose geometry if the number of domains has changed"
    );

    // Allow explicit -constant, have zero from time range
    timeSelector::addOptions(true, false);  // constant(true), zero(false)

    // Prevent volume BCs from triggering finite-area
    regionModels::allowFaModels(false);

    #include "setRootCase.H"

    // ------------------------------------------------------------------------
    // Configuration

    const bool writeCellDist    = args.found("cellDist");

    // Most of these are ignored for dry-run (not triggered anywhere)
    const bool copyZero         = args.found("copyZero");
    const bool copyUniform      = args.found("copyUniform");
    const bool decomposeSets    = !args.found("no-sets");

    const bool decomposeIfRequired = args.found("ifRequired");

    const bool doDecompFields = !args.found("no-fields");
    const bool doFiniteArea = !args.found("no-finite-area");
    const bool doLagrangian = !args.found("no-lagrangian");

    bool decomposeFieldsOnly = args.found("fields");
    bool forceOverwrite      = args.found("force");

    // Special use - emit old (2512 and earlier) edgeProcAddressing format
    // without encoded edge flips.
    if (args.found("disable-edge-encoding"))
    {
        faMeshDecomposition::allowEdgeEncoding(false);
    }

    // Set time from database
    #include "createTime.H"

    // Allow override of time (unless dry-run)
    instantList times;
    if (args.dryRun())
    {
        Info<< "\ndry-run: ignoring -copy*, -fields, -force, time selection"
            << nl;
    }
    else
    {
        if (decomposeFieldsOnly && !doDecompFields)
        {
            FatalErrorIn(args.executable())
                << "Options -fields and -no-fields are mutually exclusive"
                << " ... giving up" << nl
                << exit(FatalError);
        }

        if (!doDecompFields)
        {
            Info<< "Skip decompose of all fields" << nl;
        }
        if (!doFiniteArea)
        {
            Info<< "Skip decompose of finiteArea mesh/fields" << nl;
        }
        if (!doLagrangian)
        {
            Info<< "Skip decompose of lagrangian positions/fields" << nl;
        }

        times = timeSelector::selectIfPresent(runTime, args);
    }


    // Allow override of decomposeParDict location
    fileName decompDictFile;
    if
    (
        args.readIfPresent("decomposeParDict", decompDictFile)
     && !decompDictFile.empty() && !decompDictFile.isAbsolute()
    )
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }

    // Handle volume region selections
    #include "getAllRegionOptions.H"

    // Handle area region selections
    #include "getAllFaRegionOptions.H"

    if (!doFiniteArea)
    {
        areaRegionNames.clear();  // For consistency
    }

    const bool optRegions =
        (regionNames.size() != 1 || regionNames[0] != polyMesh::defaultRegion);

    if (regionNames.size() == 1 && regionNames[0] != polyMesh::defaultRegion)
    {
        Info<< "Using region: " << regionNames[0] << nl << endl;
    }

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = polyMesh::regionName(regionName);

        if (args.dryRun())
        {
            Info<< "dry-run: decomposing mesh " << regionName << nl << nl
                << "Create mesh..." << flush;

            domainDecompositionDryRun decompTest
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                decompDictFile,
                args.getOrDefault<label>("domains", 0),
                args.getOrDefault<word>("method", word::null)
            );

            decompTest.execute(writeCellDist, args.verbose());
            continue;
        }

        Info<< "\n\nDecomposing mesh";
        if (!regionDir.empty())
        {
            Info<< ' ' << regionName;
        }
        Info<< nl << endl;

        // Determine the existing processor count directly
        const label nProcsOld =
            fileHandler().nProcs(runTime.path(), regionDir);

        // Get requested numberOfSubdomains directly from the dictionary.
        // Note: have no mesh yet so cannot use decompositionModel::New
        const label nDomains = decompositionMethod::nDomains
        (
            IOdictionary
            (
                IOobject::selectIO
                (
                    IOobject
                    (
                        decompositionModel::canonicalName,
                        runTime.time().system(),
                        regionDir,  // region (if non-default)
                        runTime,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    decompDictFile
                )
            )
        );

        // Give file handler a chance to determine the output directory
        fileHandler().constCast().nProcs(nDomains);

        if (decomposeFieldsOnly)
        {
            // Sanity check on previously decomposed case
            if (nProcsOld != nDomains)
            {
                FatalErrorIn(args.executable())
                    << "Specified -fields, but the case was decomposed with "
                    << nProcsOld << " domains"
                    << nl
                    << "instead of " << nDomains
                    << " domains as specified in decomposeParDict" << nl
                    << exit(FatalError);
            }
        }
        else if (nProcsOld)
        {
            bool procDirsProblem = true;

            if (decomposeIfRequired && nProcsOld == nDomains)
            {
                // We can reuse the decomposition
                decomposeFieldsOnly = true;
                procDirsProblem = false;
                forceOverwrite = false;

                Info<< "Using existing processor directories" << nl;
            }

            if (optRegions)
            {
                procDirsProblem = false;
                forceOverwrite = false;
            }

            if (forceOverwrite)
            {
                Info<< "Removing " << nProcsOld
                    << " existing processor directories" << endl;

                // Remove existing processors directory
                fileNameList dirs
                (
                    fileHandler().readDir
                    (
                        runTime.path(),
                        fileName::Type::DIRECTORY
                    )
                );
                forAllReverse(dirs, diri)
                {
                    const fileName& d = dirs[diri];

                    label proci = -1;

                    if
                    (
                        d.starts_with("processor")
                     &&
                        (
                            // Collated is "processors"
                            d[9] == 's'

                            // Uncollated has integer(s) after 'processor'
                         || Foam::read(d.substr(9), proci)
                        )
                    )
                    {
                        if (fileHandler().exists(d))
                        {
                            fileHandler().rmDir(d);
                        }
                    }
                }

                procDirsProblem = false;
            }

            if (procDirsProblem)
            {
                FatalErrorIn(args.executable())
                    << "Case is already decomposed with " << nProcsOld
                    << " domains, use the -force option or manually" << nl
                    << "remove processor directories before decomposing. e.g.,"
                    << nl
                    << "    rm -rf " << runTime.path().c_str() << "/processor*"
                    << nl
                    << exit(FatalError);
            }
        }

        Info<< "Create mesh" << endl;
        domainDecomposition mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            decompDictFile
        );
        // Make sure pointMesh gets read as well
        (void)pointMesh::New(mesh, IOobject::READ_IF_PRESENT);


        // Decompose the mesh
        if (!decomposeFieldsOnly)
        {
            mesh.decomposeMesh();
            mesh.writeDecomposition(decomposeSets);

            if (writeCellDist)
            {
                const labelUList& procIds = mesh.cellToProc();

                // Write decomposition for visualization
                mesh.writeVolField("cellDist");
                //TBD: mesh.writeVTK("cellDist");

                // Write decomposition as labelList for use with 'manual'
                // decomposition method.

                IOobject io
                (
                    "cellDecomposition",
                    mesh.facesInstance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                labelIOList::writeContents(io, procIds);

                Info<< nl << "Wrote decomposition to "
                    << io.objectRelPath()
                    << " for use in manual decomposition." << endl;
            }

            fileHandler().flush();
        }

        // Finite area handling
        // - all area regions use the same volume decomposition
        HashPtrTable<faMeshDecomposition> faMeshes;
        HashTable<bool> faMeshEdgeEncoding;

        if (doFiniteArea && !areaRegionNames.empty())
        {
            const word boundaryInst =
                mesh.time().findInstance(mesh.meshDir(), "boundary");

            for (const word& areaName : areaRegionNames)
            {
                autoPtr<faMeshDecomposition> faDecompPtr;

                IOobject io
                (
                    "faBoundary",
                    boundaryInst,
                    faMesh::meshDir(mesh, areaName),
                    mesh.time(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                if (io.typeHeaderOk<faBoundaryMesh>(true))
                {
                    // Always based on the volume decomposition!
                    faDecompPtr = autoPtr<faMeshDecomposition>::New
                    (
                        areaName,
                        mesh,
                        mesh.nProcs(),
                        mesh.model()
                    );
                }

                if (faDecompPtr)
                {
                    if (!decomposeFieldsOnly)
                    {
                        // Decompose the finite-area mesh
                        auto& aMesh = faDecompPtr();
                        Info<< "\nFinite area mesh decomposition: "
                            << areaName << endl;

                        aMesh.decomposeMesh();
                        aMesh.writeDecomposition();

                        // Remember edge encoding used
                        faMeshEdgeEncoding.set
                        (
                            areaName,
                            faMeshDecomposition::allowEdgeEncoding()
                        );
                    }

                    // Cache for subsequent field decomposition
                    faMeshes.set(areaName, std::move(faDecompPtr));
                }
            }
        }


        if (copyZero)
        {
            // Copy the 0/ directory into each of the processor directories
            // with fallback of 0.orig/ directory if necessary.

            fileName inputDir(runTime.path()/"0");

            bool canCopy = fileHandler().isDir(inputDir);
            if (!canCopy)
            {
                // Try with "0.orig" instead
                inputDir.ext("orig");
                canCopy = fileHandler().isDir(inputDir);
            }

            if (canCopy)
            {
                // Avoid copying into the same directory multiple times
                // (collated format). Don't need a hash here.
                fileName prevOutputDir;
                for (label proci = 0; proci < mesh.nProcs(); ++proci)
                {
                    Time processorDb
                    (
                        Time::controlDictName,
                        args.rootPath(),
                        args.caseName()/("processor" + Foam::name(proci)),
                        false,  // No function objects
                        false   // No extra controlDict libs
                    );
                    // processorDb.setTime(runTime);

                    // Get corresponding directory name
                    // (to handle processors/)
                    const fileName outputDir
                    (
                        fileHandler().objectPath
                        (
                            IOobject
                            (
                                word::null, // name
                                "0",        // instance (time == 0)
                                processorDb
                            ),
                            word::null
                        )
                    );

                    if (outputDir != prevOutputDir)
                    {
                        Info<< "Processor " << proci
                            << ": copying \""
                            << inputDir.name() << "/\" to "
                            << runTime.relativePath(outputDir)
                            << endl;

                        fileHandler().cp(inputDir, outputDir);
                        prevOutputDir = outputDir;
                    }
                }
            }
            else
            {
                Info<< "No 0/ or 0.orig/ directory to copy" << nl;
            }
        }
        else
        {
            // Decompose field files, lagrangian, finite-area
            const auto numProcs = mesh.nProcs();

            // Cached processor meshes and maps.
            // These are only preserved if running with multiple times.
            PtrList<Time> processorDbList(numProcs);
            PtrList<fvMesh> procMeshList(numProcs);
            PtrList<labelIOList> faceProcAddressingList(numProcs);
            PtrList<labelIOList> cellProcAddressingList(numProcs);
            PtrList<labelIOList> boundaryProcAddressingList(numProcs);
            PtrList<labelIOList> pointProcAddressingList(numProcs);
            PtrList<labelIOList> pointBoundaryProcAddressingList(numProcs);

            PtrList<fvFieldDecomposer> fieldDecomposerList(numProcs);
            PtrList<pointFieldDecomposer> pointFieldDecomposerList(numProcs);


            // Cached processor meshes and maps.
            // These are only preserved if running with multiple times.
            HashPtrTable<PtrList<faMesh>> procFaMeshes;
            HashPtrTable<PtrList<labelIOList>> faFaceProcAddressing;
            HashPtrTable<PtrList<labelIOList>> faEdgeProcAddressing;
            HashPtrTable<PtrList<labelIOList>> faBoundProcAddressing;
            HashPtrTable<PtrList<faFieldDecomposer>> faFieldDecomposers;

            // Slightly wasteful, but with an *existing* finite-area
            // decomposition must scan edgeProcAddressing (from disk)
            // to know if it uses flip encoding or not.

            if
            (
                doDecompFields
             && !faMeshes.empty() && faMeshEdgeEncoding.empty()
            )
            {
                for (label proci = numProcs-1; proci >= 0; --proci)
                {
                    auto& procTime = emplaceTime(processorDbList, proci, args);

                    forAllConstIters(faMeshes, iter)
                    {
                        const word& areaName = iter.key();

                        if (faMeshEdgeEncoding.contains(areaName))
                        {
                            // Already found the encoding type
                            continue;
                        }

                        IOobject ioAddr
                        (
                            "edgeProcAddressing",
                            procTime.constant(),
                            faMesh::meshDir(regionName, areaName),
                            procTime,
                            IOobject::READ_IF_PRESENT,
                            IOobject::NO_WRITE,
                            IOobject::NO_REGISTER
                        );

                        labelList edgeProcAddr
                        (
                            labelIOList::readContents(ioAddr)
                        );

                        // Look for 0 or -ve values
                        auto i = ListOps::find_if
                        (
                            edgeProcAddr,
                            labelRange::le0()
                        );

                        if (i >= 0)
                        {
                            // A -ve value : definitely uses edge encoding.
                            // A '0' value : only occurs without encoding.
                            faMeshEdgeEncoding.set
                            (
                                areaName,
                                (edgeProcAddr[i] < 0)
                            );
                        }
                    }
                }
            }

            // Report edge-encoding (if disabled)
            if (!faMeshEdgeEncoding.empty())
            {
                bool header = false;
                forAllConstIters(faMeshEdgeEncoding, iter)
                {
                    const auto& areaName = iter.key();
                    const bool encoding = iter.val();

                    if (!encoding)
                    {
                        if (!header)
                        {
                            header = true;
                            Info<< "Area region without edge encoding:" << nl;
                        }
                        Info<< "    " << areaName;
                    }
                }
                if (header)
                {
                    Info<< endl;
                }
            }


            // Loop over all times
            forAll(times, timei)
            {
                runTime.setTime(times[timei], timei);

                Info<< nl << "Time = " << runTime.timeName() << endl;

                // Field objects at this time
                IOobjectList objects;

                if (doDecompFields)
                {
                    // List of volume mesh objects for this instance
                    objects = IOobjectList(mesh, runTime.timeName());

                    // Ignore generated fields: (cellDist)
                    objects.remove("cellDist");

                }

                // The finite-area fields (single or multiple per volume)
                HashTable<IOobjectList> faObjects;

                if (doDecompFields && doFiniteArea && faMeshes.size())
                {
                    // Lists of finite-area fields
                    faObjects.reserve(faMeshes.size());

                    forAllConstIters(faMeshes, iter)
                    {
                        const word& areaName = iter.key();

                        // The finite-area objects for this area region
                        IOobjectList objs
                        (
                            faMesh::Registry(mesh),
                            runTime.timeName(),
                            polyMesh::regionName(areaName),
                            IOobjectOption::NO_REGISTER
                        );

                        if (!objs.empty())
                        {
                            faObjects.emplace_set(areaName, std::move(objs));
                        }
                    }
                }

                // Volume/surface/internal fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                fvFieldDecomposer::fieldsCache volumeFieldCache;
                if (doDecompFields)
                {
                    volumeFieldCache.readAllFields(mesh, objects);
                }


                // Point fields
                // ~~~~~~~~~~~~

                // Read decomposed pointMesh
                const pointMesh& pMesh =
                    pointMesh::New(mesh, IOobject::READ_IF_PRESENT);

                pointFieldDecomposer::fieldsCache pointFieldCache;
                if (doDecompFields)
                {
                    pointFieldCache.readAllFields(pMesh, objects);
                }


                // Lagrangian fields
                // ~~~~~~~~~~~~~~~~~

                fileNameList cloudDirs;
                if (doDecompFields && doLagrangian)
                {
                    cloudDirs = fileHandler().readDir
                    (
                        runTime.timePath()/cloud::prefix,
                        fileName::DIRECTORY
                    );
                }

                // Particles
                PtrList<Cloud<indexedParticle>> lagrangianPositions
                (
                    cloudDirs.size()
                );
                // Particles per cell
                PtrList<List<SLList<indexedParticle*>*>> cellParticles
                (
                    cloudDirs.size()
                );

                lagrangianFieldDecomposer::fieldsCache lagrangianFieldCache
                (
                    cloudDirs.size()
                );

                label cloudI = 0;

                for (const fileName& cloudDir : cloudDirs)
                {
                    IOobjectList cloudObjects
                    (
                        mesh,
                        runTime.timeName(),
                        cloud::prefix/cloudDir,
                        IOobject::NO_REGISTER
                    );

                    // Note: look up "positions" for backwards compatibility
                    if
                    (
                        cloudObjects.found("coordinates")
                     || cloudObjects.found("positions")
                    )
                    {
                        // Read lagrangian particles
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~

                        Info<< "Identified lagrangian data set: "
                            << cloudDir << endl;

                        lagrangianPositions.set
                        (
                            cloudI,
                            new Cloud<indexedParticle>
                            (
                                mesh,
                                cloudDir,
                                false
                            )
                        );


                        // Sort particles per cell
                        // ~~~~~~~~~~~~~~~~~~~~~~~

                        cellParticles.set
                        (
                            cloudI,
                            new List<SLList<indexedParticle*>*>
                            (
                                mesh.nCells(),
                                static_cast<SLList<indexedParticle*>*>(nullptr)
                            )
                        );

                        label i = 0;

                        for (indexedParticle& p : lagrangianPositions[cloudI])
                        {
                            p.index() = i++;

                            label celli = p.cell();

                            // Check
                            if (celli < 0 || celli >= mesh.nCells())
                            {
                                FatalErrorIn(args.executable())
                                    << "Illegal cell number " << celli
                                    << " for particle with index "
                                    << p.index()
                                    << " at position "
                                    << p.position() << nl
                                    << "Cell number should be between 0 and "
                                    << mesh.nCells()-1 << nl
                                    << "On this mesh the particle should"
                                    << " be in cell "
                                    << mesh.findCell(p.position())
                                    << exit(FatalError);
                            }

                            if (!cellParticles[cloudI][celli])
                            {
                                cellParticles[cloudI][celli] =
                                    new SLList<indexedParticle*>();
                            }

                            cellParticles[cloudI][celli]->append(&p);
                        }

                        // Read fields
                        // ~~~~~~~~~~~

                        IOobjectList lagrangianObjects
                        (
                            mesh,
                            runTime.timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            IOobject::NO_REGISTER
                        );

                        lagrangianFieldCache.readAllFields
                        (
                            cloudI,
                            lagrangianObjects
                        );

                        ++cloudI;
                    }
                }

                lagrangianPositions.resize(cloudI);
                cellParticles.resize(cloudI);
                lagrangianFieldCache.resize(cloudI);


                // Finite-area (area/edge) fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                HashPtrTable<faFieldDecomposer::fieldsCache> areaFieldCaches;

                if (doDecompFields && doFiniteArea)
                {
                    forAllConstIters(faObjects, iter)
                    {
                        const word& areaName = iter.key();
                        const auto& objs = iter.val();

                        if
                        (
                            const auto meshIter = faMeshes.cfind(areaName);
                            (meshIter.good() && !objs.empty())
                        )
                        {
                            const faMesh& aMesh = *(meshIter.val());
                            auto& cache = areaFieldCaches.emplace_set(areaName);

                            cache.readAllFields(aMesh, objs);
                        }
                    }
                }

                Info<< endl;

                // Split the fields over processors
                for
                (
                    label proci = 0;
                    doDecompFields && proci < numProcs;
                    ++proci
                )
                {
                    Info<< "Processor " << proci << ": field transfer" << endl;

                    auto& processorDb =
                        emplaceTime(processorDbList, proci, args);
                    processorDb.setTime(runTime);

                    // Read the mesh
                    const auto& procMesh = procMeshList.try_emplace
                    (
                        proci,
                        IOobject
                        (
                            regionName,
                            processorDb.timeName(),
                            processorDb
                        )
                    );

                    const auto& faceProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "faceProcAddressing",
                        faceProcAddressingList
                    );

                    const auto& cellProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "cellProcAddressing",
                        cellProcAddressingList
                    );

                    const auto& boundaryProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "boundaryProcAddressing",
                        boundaryProcAddressingList
                    );


                    // FV fields: volume, surface, internal
                    {
                        if (!fieldDecomposerList.test(proci))
                        {
                            fieldDecomposerList.set
                            (
                                proci,
                                new fvFieldDecomposer
                                (
                                    mesh,
                                    procMesh,
                                    faceProcAddressing,
                                    cellProcAddressing,
                                    boundaryProcAddressing
                                )
                            );
                        }

                        volumeFieldCache.decomposeAllFields
                        (
                            fieldDecomposerList[proci]
                        );

                        if (times.size() == 1)
                        {
                            // Clear cached decomposer
                            fieldDecomposerList.set(proci, nullptr);
                        }
                    }


                    // Point fields
                    if (!pointFieldCache.empty())
                    {
                        const auto& pointProcAddressing = procAddressing
                        (
                            procMeshList,
                            proci,
                            "pointProcAddressing",
                            pointProcAddressingList
                        );

                        const pointMesh& procPMesh =
                            pointMesh::New(procMesh, IOobject::READ_IF_PRESENT);

                        if (!pointBoundaryProcAddressingList.test(proci))
                        {
                            pointBoundaryProcAddressingList.set
                            (
                                proci,
                                autoPtr<labelIOList>::New
                                (
                                    IOobject
                                    (
                                        "boundaryProcAddressing",
                                        procMesh.facesInstance(),
                                        polyMesh::meshSubDir
                                       /pointMesh::meshSubDir,
                                        procPMesh.thisDb(),
                                        IOobject::READ_IF_PRESENT,
                                        IOobject::NO_WRITE,
                                        IOobject::NO_REGISTER
                                    ),
                                    boundaryProcAddressing
                                )
                            );
                        }
                        const auto& pointBoundaryProcAddressing =
                            pointBoundaryProcAddressingList[proci];


                        if (!pointFieldDecomposerList.test(proci))
                        {
                            pointFieldDecomposerList.set
                            (
                                proci,
                                new pointFieldDecomposer
                                (
                                    pMesh,
                                    procPMesh,
                                    pointProcAddressing,
                                    pointBoundaryProcAddressing
                                )
                            );
                        }

                        pointFieldCache.decomposeAllFields
                        (
                            pointFieldDecomposerList[proci]
                        );

                        if (times.size() == 1)
                        {
                            // Early deletion
                            pointBoundaryProcAddressingList.set
                            (
                                proci,
                                nullptr
                            );
                            pointProcAddressingList.set(proci, nullptr);
                            pointFieldDecomposerList.set(proci, nullptr);
                        }
                    }


                    // If there is lagrangian data write it out
                    forAll(lagrangianPositions, cloudi)
                    {
                        if (lagrangianPositions[cloudi].size())
                        {
                            lagrangianFieldDecomposer fieldDecomposer
                            (
                                mesh,
                                procMesh,
                                faceProcAddressing,
                                cellProcAddressing,
                                cloudDirs[cloudi],
                                lagrangianPositions[cloudi],
                                cellParticles[cloudi]
                            );

                            // Lagrangian fields
                            lagrangianFieldCache.decomposeAllFields
                            (
                                cloudi,
                                cloudDirs[cloudi],
                                fieldDecomposer
                            );
                        }
                    }

                    if (doDecompFields)
                    {
                        // Decompose "uniform" directory in the time region
                        // directory
                        decomposeUniform
                        (
                            copyUniform, mesh, processorDb, regionDir
                        );

                        // For a multi-region case, also decompose "uniform"
                        // directory in the time directory
                        if (regionNames.size() > 1 && regioni == 0)
                        {
                            decomposeUniform(copyUniform, mesh, processorDb);
                        }
                    }

                    if (times.size() == 1)
                    {
                        // Early deletion
                        boundaryProcAddressingList.set(proci, nullptr);
                        cellProcAddressingList.set(proci, nullptr);
                        faceProcAddressingList.set(proci, nullptr);
                    }

                    // Finite-area fields
                    for (const auto& iter : areaFieldCaches.csorted())
                    {
                        const word& areaName = iter.key();
                        const auto& areaCache = *(iter.val());

                        // Serial mesh:
                        const faMesh& aMesh = *(faMeshes[areaName]);

                        // List of proc meshes:
                        auto& faProcMeshList =
                            procFaMeshes.try_emplace(areaName, numProcs);

                        auto& procFaMesh =
                            faProcMeshList
                            .try_emplace(proci, areaName, procMesh);

                        const auto& faFaceProcAddr =
                            procAddressing
                            (
                                faProcMeshList,
                                proci,
                                "faceProcAddressing",
                                faFaceProcAddressing.try_emplace(areaName)
                            );

                        const auto& faBoundProcAddr =
                            procAddressing
                            (
                                faProcMeshList,
                                proci,
                                "boundaryProcAddressing",
                                faBoundProcAddressing.try_emplace(areaName)
                            );

                        const auto& faEdgeProcAddr =
                            procAddressing
                            (
                                faProcMeshList,
                                proci,
                                "edgeProcAddressing",
                                faEdgeProcAddressing.try_emplace(areaName)
                            );


                        auto& faFieldDecomposerList =
                            faFieldDecomposers.try_emplace(areaName, numProcs);

                        if (!faFieldDecomposerList.test(proci))
                        {
                            faFieldDecomposerList.emplace
                            (
                                proci,
                                //
                                aMesh,
                                procFaMesh,
                                faEdgeProcAddr,
                                faFaceProcAddr,
                                faBoundProcAddr,
                                // noEdgeEncoding
                                (!faMeshEdgeEncoding.lookup(areaName, true))
                            );
                        }

                        auto& fieldDecomposer = faFieldDecomposerList[proci];

                        areaCache.decomposeAllFields
                        (
                            fieldDecomposer,
                            args.verbose()  // report
                        );
                    }

                    // We have cached all the constant mesh data for the current
                    // processor. This is only important if running with
                    // multiple times, otherwise it is just extra storage.
                    if (times.size() == 1)
                    {
                        forAllIters(faFieldDecomposers, iter)
                        {
                            iter.val()->set(proci, nullptr);
                        }
                        forAllIters(faEdgeProcAddressing, iter)
                        {
                            iter.val()->set(proci, nullptr);
                        }
                        forAllIters(faFaceProcAddressing, iter)
                        {
                            iter.val()->set(proci, nullptr);
                        }
                        forAllIters(faBoundProcAddressing, iter)
                        {
                            iter.val()->set(proci, nullptr);
                        }
                        forAllIters(procFaMeshes, iter)
                        {
                            iter.val()->set(proci, nullptr);
                        }

                        procMeshList.set(proci, nullptr);
                        processorDbList.set(proci, nullptr);
                    }
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
