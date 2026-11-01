/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "triSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "edgeHashes.H"
#include "triSurfaceFields.H"
#include "Time.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceMesh, 0);
    addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);
}

Foam::word Foam::triSurfaceMesh::meshSubDir = "triSurface";


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Helper function for isSurfaceClosed
bool addFaceToEdge
(
    const Foam::edge& e,
    Foam::EdgeMap<Foam::label>& facesPerEdge
)
{
    using namespace Foam;

    //auto& count = facesPerEdge(e, 0);  // lookup or new entry
    //if (count == 2)
    //{
    //    return false;
    //}
    //++count;
    //return true;

    if (auto fnd = facesPerEdge.find(e); fnd.good())
    {
        auto& count = fnd.val();
        const int dir = edge::compare(e, fnd.key());
        if (dir == 1)
        {
            // Incorrect order. Mark with special value
            count = -1;
        }
        else if (dir == 0)
        {
            FatalErrorInFunction
                << "Incorrect matched edge " << fnd.key()
                << " to edge " << e << exit(FatalError);
            return false;
        }
        else if (count != -1)
        {
            if (count == 2)
            {
                return false;
            }
            ++count;
        }
    }
    else
    {
        facesPerEdge.insert(e, 1);
    }

    return true;
}


// Calculate the bounding sphere for a triFace
inline Foam::scalar faceBoundingSphere
(
    const Foam::point& faceCentre,
    const Foam::triFace& tri,
    const Foam::UList<Foam::point>& points
)
{
    using namespace Foam;

    const scalar a = faceCentre.distSqr(points[tri.a()]);
    const scalar b = faceCentre.distSqr(points[tri.b()]);
    const scalar c = faceCentre.distSqr(points[tri.c()]);

    // Add a bit to make sure all points are tested inside
    return (std::max(std::max(a, b), c) + Foam::sqr(SMALL));
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triSurfaceMesh::isSurfaceClosed() const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::isSurfaceClosed:"
            << " surface:" << searchableSurface::name()
            << " determining closedness for surface with "
            << triSurface::size() << " triangles" << endl;
    }

    const triSurface& s = *this;
    const auto& pts = s.points();

/*
    // Construct pointFaces. Let's hope surface has compact point
    // numbering ...
    labelListList pointFaces;
    invertManyToMany(pts.size(), s, pointFaces);

    // Loop over all faces surrounding point. Count edges emanating from point.
    // Every edge should be used by two faces exactly.
    // To prevent doing work twice per edge only look at edges to higher
    // point
    EdgeMap<label> facesPerEdge;
    forAll(pointFaces, pointi)
    {
        const auto& pFaces = pointFaces[pointi];

        facesPerEdge.clear();
        for (const label facei : pFaces)
        {
            const auto& f = s[facei];
            const label fp = f.find(pointi);

            // Something weird: if I expand the code of addFaceToEdge in both
            // below instances it gives a segmentation violation on some
            // surfaces. Compiler (4.3.2) problem?


            // Forward edge
            const label nextPointi = f[f.fcIndex(fp)];

            if (nextPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, nextPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    if (debug)
                    {
                        Pout<< "triSurfaceMesh::isSurfaceClosed :"
                            << " surface is open" << endl;
                    }
                    return false;
                }
            }

            // Reverse edge
            const label prevPointi = f[f.rcIndex(fp)];

            if (prevPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, prevPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    if (debug)
                    {
                        Pout<< "triSurfaceMesh::isSurfaceClosed :"
                            << " surface is open" << endl;
                    }
                    return false;
                }
            }
        }

        // Check for any edges used only once.
        forAllConstIters(facesPerEdge, iter)
        {
            if (iter.val() == -1)
            {
                // Special value for incorrect orientation
                if (debug)
                {
                    Pout<< "triSurfaceMesh::isSurfaceClosed :"
                        << " surface is closed but has inconsistent"
                        << " face orientation" << endl;
                }
                WarningInFunction << "Surface " << searchableSurface::name()
                    << " is closed but has inconsistent face orientation"
                    << " at edge " << pts[iter.key().first()]
                    << pts[iter.key().second()] << endl;
                return false;
            }
            else if (iter.val() != 2)
            {
                if (debug)
                {
                    Pout<< "triSurfaceMesh::isSurfaceClosed :"
                        << " surface is open" << endl;
                }
                return false;
            }
        }
    }
*/

    EdgeMap<label> facesPerEdge;
    facesPerEdge.reserve(s.size());

    for (const auto& f : s)
    {
        forAll(f, fp)
        {
            // Count number of occurences of edge between fp and fp+1
            const bool okFace = addFaceToEdge(f.edge(fp), facesPerEdge);

            if (!okFace)
            {
                if (debug)
                {
                    Pout<< "triSurfaceMesh::isSurfaceClosed :"
                        << " surface:" << searchableSurface::name()
                        << " surface is non-manifold" << endl;
                }
                return false;
            }
        }
    }


    // Check for any edges with inconsistent normal orientation,
    // open surfaces, edges used only once etc.

    edge warnEdge(-1, -1);
    label nWarnings = 0;

    forAllConstIters(facesPerEdge, iter)
    {
        if (iter.val() < 0)  // Special value (incorrect orientation)
        {
            if (!nWarnings)
            {
                warnEdge = iter.key();
            }
            ++nWarnings;
        }
        else if (iter.val() != 2)
        {
            if (debug)
            {
                Pout<< "triSurfaceMesh::isSurfaceClosed :"
                    << " surface:" << searchableSurface::name()
                    << " surface is open" << endl;
            }
            // not closed
            return false;
        }
    }

    if (nWarnings)
    {
        WarningInFunction
            << "Surface " << searchableSurface::name()
            << " is closed but has inconsistent face orientation at "
            << nWarnings << " edges." << nl
            << "    first edge: "
            << pts[warnEdge.first()] << pts[warnEdge.second()] << nl
            << "    This means it probably cannot be used"
            << " for inside/outside queries." << endl;
    }

    if (debug)
    {
        Pout<< "triSurfaceMesh::isSurfaceClosed :"
            << " surface:" << searchableSurface::name()
            << " surface is closed" << endl;
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
    searchableSurface(io),
    objectRegistry
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(s),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    // Local bound box (no reduction)
    bounds() = treeBoundBox(triSurface::points());
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(static_cast<const searchableSurface&>(*this), dictionary::null),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    // Local bound box (no reduction)
    bounds() = treeBoundBox(triSurface::points());
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(static_cast<const searchableSurface&>(*this), dict),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    // Read with supplied file name instead of objectPath/filePath
    if (dict.readIfPresent("file", fName_, keyType::LITERAL))
    {
        fName_ = triSurface::relativeFilePath
        (
            static_cast<const searchableSurface&>(*this),
            fName_,
            true
        );
    }

    // Report optional scale factor (eg, mm to m),
    // which was already applied during triSurface construction
    if (scalar factor; dict.readIfPresent("scale", factor) && factor > 0)
    {
        Info<< searchableSurface::name()
            << " : using scale " << factor << endl;
    }

    // Local bound box (no reduction)
    bounds() = treeBoundBox(triSurface::points());

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const readAction r)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    const bool searchGlobal
    (
        r == readAction::localOrGlobal
     || r == readAction::masterOnly
    );

    // Check IO flags
    if (io.isAnyRead())
    {
        const fileName actualFile
        (
            searchGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        // const fileName actualFile
        // (
        //     io.filePath(searchGlobal, typeName)
        // );

        if (debug)
        {
            Pout<< "triSurfaceMesh(const IOobject& io) :"
                << " loading surface " << io.objectPath()
                << " local filePath:" << io.localFilePath(typeName)
                << " from:" << actualFile << endl;
        }

        if (searchGlobal && Pstream::parRun())
        {
            // Check where surface was found
            const fileName localFile(io.localFilePath(typeName));

            if
            (
                r == masterOnly
             && ((actualFile.empty() || actualFile != localFile))
            )
            {
                // Found undecomposed surface. Load on master only
                if (Pstream::master())
                {
                    triSurface s2(actualFile);
                    triSurface::transfer(s2);
                }
                Pstream::broadcast(triSurface::patches());
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
            else
            {
                // Read on all processors
                triSurface s2(actualFile);
                triSurface::transfer(s2);
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
        }
        else
        {
            // Read on all processors
            triSurface s2(actualFile);
            triSurface::transfer(s2);
            if (debug)
            {
                Pout<< "triSurfaceMesh(const IOobject& io) :"
                    << " loaded triangles:" << triSurface::size() << endl;
            }
        }
    }

    // Local bound box (no reduction)
    bounds() = treeBoundBox(triSurface::points());
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict,
    const readAction r
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    const bool searchGlobal
    (
        r == readAction::localOrGlobal
     || r == readAction::masterOnly
    );

    if (io.isAnyRead())
    {
        // Surface type (optional)
        const word surfType(dict.getOrDefault<word>("fileType", word::null));

        // Scale factor (optional)
        const scalar scaleFactor(dict.getOrDefault<scalar>("scale", 0));

        fileName actualFile
        (
            searchGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        // const fileName actualFile
        // (
        //     io.filePath(searchGlobal, typeName)
        // );

        // Reading from supplied file name instead of objectPath/filePath
        if (dict.readIfPresent("file", fName_, keyType::LITERAL))
        {
            fName_ = relativeFilePath
            (
                static_cast<const searchableSurface&>(*this),
                fName_,
                searchGlobal
            );
            actualFile = fName_;
        }

        if (debug)
        {
            Pout<< "triSurfaceMesh(const IOobject& io, const dictionary&) :"
                << " loading surface " << io.objectPath()
                << " local filePath:" << io.localFilePath(typeName)
                << " from:" << actualFile << endl;
        }


        if (searchGlobal && Pstream::parRun())
        {
            // Check where surface was found. Bit tricky:
            // - master will have actualFile (in parent directory)
            //   different from localFilePath (in processor0/)
            // - slave might have actualFile empty and localFile empty

            const fileName localFile(io.localFilePath(typeName));

            if
            (
                r == masterOnly
             && ((actualFile.empty() || actualFile != localFile))
            )
            {
                // Surface not loaded from processor directories -> undecomposed
                // surface. Load on master only
                if (Pstream::master())
                {
                    triSurface s2(actualFile, surfType, scaleFactor);
                    triSurface::transfer(s2);
                }
                Pstream::broadcast(triSurface::patches());
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
            else
            {
                // Read on all processors
                triSurface s2(actualFile, surfType, scaleFactor);
                triSurface::transfer(s2);
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
        }
        else
        {
            // Read on all processors
            triSurface s2(actualFile, surfType, scaleFactor);
            triSurface::transfer(s2);
            if (debug)
            {
                Pout<< "triSurfaceMesh(const IOobject& io) :"
                    << " loaded triangles:" << triSurface::size() << endl;
            }
        }


        // Report optional scale factor (eg, mm to m),
        // which was already applied during triSurface reading
        if (scaleFactor > 0)
        {
            Info<< searchableSurface::name()
                << " : using scale " << scaleFactor << endl;
        }
    }

    // Local bound box (no reduction)
    bounds() = treeBoundBox(triSurface::points());

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{
    clearOut();
}


void Foam::triSurfaceMesh::clearOut()
{
    // Do not clear closedness status
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::coordinates() const
{
    auto tpts = tmp<pointField>::New();
    triSurface::centres(tpts.ref());
    return tpts;
}


void Foam::triSurfaceMesh::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    const triSurface& s = *this;
    const auto& pts = s.points();

    centres.resize_nocopy(s.size());
    radiusSqr.resize_nocopy(s.size());

    triSurface::centres(centres);

    forAll(s, facei)
    {
        radiusSqr[facei] = faceBoundingSphere(centres[facei], s[facei], pts);
    }

    if (debug)
    {
        Pout<< "triSurfaceMesh::boundingSpheres :"
            << " surface:" << searchableSurface::name()
            << " returning " << centres.size() << " values"
            << endl;
    }
}


Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::points() const
{
    return triSurface::points();
}


bool Foam::triSurfaceMesh::overlaps(const boundBox& bb) const
{
    return tree().overlaps(treeBoundBox(bb));
}


void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::movePoints :"
            << " surface:" << searchableSurface::name()
            << " moving at time " << objectRegistry::time().timeName()
            << endl;
    }

    // Preserve topological point status (surfaceClosed_, outsideVolType_)

    // Update local information (instance, event number)
    searchableSurface::instance() = objectRegistry::time().timeName();
    objectRegistry::instance() = searchableSurface::instance();

    const label event = getEvent();
    searchableSurface::eventNo() = event;
    objectRegistry::eventNo() = searchableSurface::eventNo();

    // Clear additional addressing
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::movePoints(newPoints);

    // Local bound box (no reduction)
    bounds() = treeBoundBox(triSurface::points());

    if (debug)
    {
        Pout<< "triSurfaceMesh::movePoints: finished moving points" << endl;
    }
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::triSurfaceMesh::edgeTree() const
{
    if (!edgeTree_)
    {
        if (debug)
        {
            Pout<< "triSurfaceMesh::edgeTree :"
                << " surface:" << searchableSurface::name()
                << " constructing tree for " << nEdges() - nInternalEdges()
                << " boundary edges" << endl;
        }

        // Boundary edges
        const labelRange bEdges(nInternalEdges(), nBoundaryEdges());

        // Random number generator. Bit dodgy since not exactly random ;-)
        constexpr label rand_seed = 65431;

        treeBoundBox bb(point::zero);

        if (bEdges.size())
        {
            label nPoints;
            PatchTools::calcBounds(*this, bb, nPoints);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.

            bb.inflate_random(rand_seed, 1e-4, ROOTVSMALL);
        }


        if (debug)
        {
            Pout<< "triSurfaceMesh::edgeTree : "
                << "calculating edge tree for bb:" << bb << endl;
        }

        const scalar oldTol = indexedOctreeBase::perturbTol(tolerance());

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                // Boundary edges
                treeDataEdge(edges(), localPoints(), bEdges),

                bb,                 // bb
                maxTreeDepth(),     // maxLevel
                10,                 // leafsize
                3.0                 // duplicity
            )
        );

        indexedOctreeBase::perturbTol(oldTol);

        if (debug)
        {
            Pout<< "triSurfaceMesh::edgeTree :"
                << " finished constructing tree for "
                << nEdges() - nInternalEdges()
                << " boundary edges" << endl;
        }
    }

    return *edgeTree_;
}


const Foam::wordList& Foam::triSurfaceMesh::regions() const
{
    if (regions_.empty())
    {
        regions_.resize(patches().size());
        forAll(regions_, regioni)
        {
            regions_[regioni] = patches()[regioni].name();
        }
    }
    return regions_;
}


bool Foam::triSurfaceMesh::hasVolumeType() const
{
    if (surfaceClosed_ == -1)
    {
        if (isSurfaceClosed())
        {
            surfaceClosed_ = 1;
        }
        else
        {
            surfaceClosed_ = 0;
        }
    }

    return surfaceClosed_ == 1;
}


Foam::volumeType Foam::triSurfaceMesh::outsideVolumeType() const
{
    if (outsideVolType_ == volumeType::UNKNOWN)
    {
        // Get point outside bounds()
        const point outsidePt(bounds().max() + 0.5*bounds().span());

        if (debug)
        {
            Pout<< "triSurfaceMesh::outsideVolumeType :"
                << " surface:" << searchableSurface::name()
                << " triggering outsidePoint:" << outsidePt
                << " orientation" << endl;
        }

        //outsideVolType_ = tree().shapes().getVolumeType(tree(), outsidePt);
        // Note: do not use tree directly so e.g. distributedTriSurfaceMesh
        //       has opportunity to intercept
        List<volumeType> outsideVolTypes;
        getVolumeType(pointField(1, outsidePt), outsideVolTypes);
        outsideVolType_ = outsideVolTypes[0];

        if (debug)
        {
            Pout<< "triSurfaceMesh::outsideVolumeType :"
                << " finished outsidePoint:" << outsidePt
                << " orientation:" << volumeType::names[outsideVolType_]
                << endl;
        }
    }

    return outsideVolType_;
}


void Foam::triSurfaceMesh::flip()
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::flip :"
            << " surface:" << searchableSurface::name()
            << " with current orientation "
            << volumeType::names[outsideVolType_]
            << endl;
    }

    // Don't bother getting nearest etc. Just flip the triangles.

    // triSurface
    {
        triSurface& s = *this;
        s.clearOut();
        for (auto& tri : s)
        {
            tri.flip();
        }
    }

    // triSurfaceRegionSearch (if cached volume type)
    triSurfaceRegionSearch::flip();

    // edge tree not relevant

    if (hasVolumeType())
    {
        // outsideVolType_
        if (outsideVolType_ == volumeType::INSIDE)
        {
            outsideVolType_ = volumeType::OUTSIDE;
        }
        else if (outsideVolType_ == volumeType::OUTSIDE)
        {
            outsideVolType_ = volumeType::INSIDE;
        }
    }
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        scalar maxDistSqr = 0;
        if (!samples.empty())
        {
            if (nearestDistSqr.empty())
            {
                maxDistSqr = searchableSurface::defaultMaxDistSqr();
            }
            else
            {
                maxDistSqr = Foam::max(nearestDistSqr);
            }
        }

        Pout<< "triSurfaceMesh::findNearest :"
            << " surface:" << searchableSurface::name()
            << " trying to find nearest for " << samples.size()
            << " samples with max sphere " << std::sqrt(maxDistSqr) << endl;
    }

    triSurfaceSearch::findNearest(samples, nearestDistSqr, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findNearest :"
            << " finished trying to find nearest for " << samples.size()
            << " samples" << endl;
    }
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const labelList& regionIndices,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        scalar maxDistSqr = 0;
        if (!samples.empty())
        {
            if (nearestDistSqr.empty())
            {
                maxDistSqr = searchableSurface::defaultMaxDistSqr();
            }
            else
            {
                maxDistSqr = Foam::max(nearestDistSqr);
            }
        }

        Pout<< "triSurfaceMesh::findNearest :"
            << " surface:" << searchableSurface::name()
            << " trying to find nearest and region for " << samples.size()
            << " samples with max sphere " << std::sqrt(maxDistSqr) << endl;
    }

    triSurfaceRegionSearch::findNearest
    (
        samples,
        nearestDistSqr,
        regionIndices,
        info
    );

    if (debug)
    {
        Pout<< "triSurfaceMesh::findNearest :"
            << " finished trying to find nearest and region for "
            << samples.size() << " samples" << endl;
    }
}


void Foam::triSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLine :"
            << " surface:" << searchableSurface::name()
            << " intersecting with "
            << start.size() << " rays" << endl;
    }
    triSurfaceSearch::findLine(start, end, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLine :"
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::triSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAny :"
            << " surface:" << searchableSurface::name()
            << " intersecting with "
            << start.size() << " rays" << endl;
    }
    triSurfaceSearch::findLineAny(start, end, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAny :"
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::triSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAll :"
            << " surface:" << searchableSurface::name()
            << " intersecting with "
            << start.size() << " rays" << endl;
    }
    triSurfaceSearch::findLineAll(start, end, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAll :"
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::triSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::getRegion :"
            << " surface:" << searchableSurface::name()
            << " getting region for "
            << info.size() << " triangles" << endl;
    }

    region.resize_nocopy(info.size());

    const triSurface& s = *this;
    const label numFaces = s.size();

    forAll(info, i)
    {
        const label facei = info[i].index();

        region[i] =
        (
            (info[i].hit() && facei >= 0 && facei < numFaces)
          ? s[facei].region() : -1
        );
    }

    if (debug)
    {
        Pout<< "triSurfaceMesh::getRegion :"
            << " finished getting region for "
            << info.size() << " triangles" << endl;
    }
}


namespace
{

// Search neighbour triangles for one with a better quality
template<class FaceFaces>
Foam::label betterQualityNeighbour
(
    const Foam::triSurface& s,
    Foam::label facei,
    const FaceFaces& faceFaces
)
{
    using namespace Foam;

    label bestNbr = -1;

    if (facei >= 0 && facei < s.size())
    {
        const auto& pts = s.points();

        auto qual = s[facei].tri(pts).quality();

        for (label nbri : faceFaces[facei])
        {
            auto nbrQual = s[nbri].tri(pts).quality();

            if (qual < nbrQual)
            {
                qual = nbrQual;
                bestNbr = nbri;
            }
        }
    }

    return bestNbr;
}

} // End anonymous namespace


void Foam::triSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::getNormal :"
            << " surface:" << searchableSurface::name()
            << " getting normal for "
            << info.size() << " triangles" << endl;
    }

    const triSurface& s = *this;
    const auto& pts = s.points();
    const label numFaces = s.size();

    normal.resize_nocopy(info.size());

    // Note: keep this code synchronized with the measure() method

    // Track which faces (if any) need to be adjusted for poor quality
    bitSet fixQuality;

    if (minQuality_ >= 0)
    {
        fixQuality.resize(numFaces);
    }

    forAll(info, i)
    {
        if
        (
            const label facei = info[i].index();
            (info[i].hit() && facei >= 0 && facei < numFaces)
        )
        {
            // Uncached
            normal[i] = s[facei].unitNormal(pts);

            if (minQuality_ >= 0)
            {
                scalar qual = s[facei].tri(pts).quality();

                if (qual < minQuality_)
                {
                    fixQuality.set(facei);
                }
            }
        }
        else
        {
            // No useful information possible
            normal[i] = Foam::zero{};
        }
    }

    // Make sure we don't use triangles with low quality since normal is
    // not reliable. Only apply to ones that we already know to be poor
    // (checked above) and otherwise avoid creating face-face addressing.

    label nAdjusted = 0;

    if (fixQuality.any())
    {
        const auto& faceFaces = s.faceFaces();

        forAll(info, i)
        {
            if
            (
                const label facei = info[i].index();
                (info[i].hit() && fixQuality.test(facei))
            )
            {
                label bestNbr = betterQualityNeighbour(s, facei, faceFaces);

                if (bestNbr >= 0)
                {
                    normal[i] = s[bestNbr].unitNormal(pts);
                    ++nAdjusted;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "triSurfaceMesh::getNormal :"
            << " finished getting normal for "
            << info.size() << " triangles (adjusted "
            << nAdjusted << " for quality)" << endl;
    }
}


Foam::label Foam::triSurfaceMesh::measure
(
    List<vector>* centres,
    List<scalar>* boundingRadiusSqr,
    List<label>* regions,
    List<vector>* normals
) const
{
    const triSurface& s = *this;
    const auto& pts = s.points();

    const label output_size = s.size();

    bool usedParams = false;

    // Extract face centres
    if (centres)
    {
        usedParams = true;
        s.centres(*centres);
    }

    // Extract bounding sphere radius squared
    if (boundingRadiusSqr)
    {
        usedParams = true;
        auto& radiusSqr = *boundingRadiusSqr;
        radiusSqr.resize_nocopy(output_size);

        if (centres)
        {
            // Reuse extracted face centres
            forAll(s, facei)
            {
                radiusSqr[facei] = faceBoundingSphere
                (
                    (*centres)[facei],
                    s[facei],
                    pts
                );
            }
        }
        else
        {
            // Generate face centre on the fly
            forAll(s, facei)
            {
                radiusSqr[facei] = faceBoundingSphere
                (
                    s[facei].centre(pts),
                    s[facei],
                    pts
                );
            }
        }
    }

    // Extract regions
    if (regions)
    {
        usedParams = true;
        s.regionIDs(*regions);
    }

    // Extract normals.
    // Note: keep this code synchronized with the getNormal() method
    if (normals)
    {
        usedParams = true;
        auto& normal = *normals;
        normal.resize_nocopy(output_size);

        // Track which faces (if any) need to be adjusted for poor quality
        bitSet fixQuality;

        if (minQuality_ >= 0)
        {
            fixQuality.resize(output_size);
        }

        forAll(s, facei)
        {
            // Uncached
            normal[facei] = s[facei].unitNormal(pts);

            if (minQuality_ >= 0)
            {
                scalar qual = s[facei].tri(pts).quality();

                if (qual < minQuality_)
                {
                    fixQuality.set(facei);
                }
            }
        }

        // Make sure we don't use triangles with low quality since normal is
        // not reliable. Only apply to ones that we already know to be poor
        // (checked above).

        if (fixQuality.any())
        {
            const auto& faceFaces = s.faceFaces();

            for (label facei : fixQuality)
            {
                label bestNbr = betterQualityNeighbour(s, facei, faceFaces);

                if (bestNbr >= 0)
                {
                    normal[facei] = s[bestNbr].unitNormal(pts);
                }
            }
        }
    }

    return (usedParams ? output_size : 0);
}


void Foam::triSurfaceMesh::setField(const labelList& values)
{
    auto* fldPtr = getObjectPtr<triSurfaceLabelField>("values");

    if (fldPtr)
    {
        (*fldPtr).field() = values;
    }
    else
    {
        fldPtr = new triSurfaceLabelField
        (
            IOobject
            (
                "values",
                objectRegistry::time().timeName(),  // instance
                meshSubDir,                         // local
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            *this,
            dimless,
            labelField(values)
        );

        // Store field on triMesh
        fldPtr->store();
    }
    if (debug)
    {
        Pout<< "triSurfaceMesh::setField :"
            << " surface:" << searchableSurface::name()
            << " finished setting field for "
            << values.size() << " triangles" << endl;
    }
}


void Foam::triSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    const auto* fldPtr = getObjectPtr<triSurfaceLabelField>("values");

    if (fldPtr)
    {
        const auto& fld = *fldPtr;

        values.setSize(info.size());

        forAll(info, i)
        {
            if (info[i].hit())
            {
                values[i] = fld[info[i].index()];
            }
        }
    }
    if (debug)
    {
        Pout<< "triSurfaceMesh::setField :"
            << " surface:" << searchableSurface::name()
            << " finished getting field for "
            << info.size() << " triangles" << endl;
    }
}


void Foam::triSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    const scalar oldTol = indexedOctreeBase::perturbTol(tolerance());

    if (debug)
    {
        Pout<< "triSurfaceMesh::getVolumeType :"
            << " surface:" << searchableSurface::name()
            << " finding orientation for " << points.size()
            << " samples" << endl;
    }

    volType.resize_nocopy(points.size());

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        if (tree().bb().contains(pt))
        {
            // Use cached volume type per each tree node
            volType[pointi] = tree().getVolumeType(pt);
        }
        else if (hasVolumeType())
        {
            // Precalculate and cache value for this outside point
            if (outsideVolType_ == volumeType::UNKNOWN)
            {
                outsideVolType_ = tree().shapes().getVolumeType(tree(), pt);
            }
            volType[pointi] = outsideVolType_;
        }
        else
        {
            // Have to calculate directly as outside the octree
            volType[pointi] = tree().shapes().getVolumeType(tree(), pt);
        }
    }

    indexedOctreeBase::perturbTol(oldTol);

    if (debug)
    {
        Pout<< "triSurfaceMesh::getVolumeType :"
            << " finished finding orientation for " << points.size()
            << " samples" << endl;
    }
}


bool Foam::triSurfaceMesh::writeObject
(
    IOstreamOption,
    const bool writeOnProc
) const
{
    const Time& runTime = searchableSurface::time();
    const fileName& instance = searchableSurface::instance();

    if
    (
        instance != runTime.timeName()
     && instance != runTime.system()
     && instance != runTime.caseSystem()
     && instance != runTime.constant()
     && instance != runTime.caseConstant()
    )
    {
        const_cast<triSurfaceMesh&>(*this).searchableSurface::instance() =
            runTime.timeName();
        const_cast<triSurfaceMesh&>(*this).objectRegistry::instance() =
            runTime.timeName();
    }

    fileName fullPath;
    if (fName_.size())
    {
        // Override file name

        fullPath = fName_;

        fullPath.expand();
        if (!fullPath.isAbsolute())
        {
            // Add directory from regIOobject
            fullPath = searchableSurface::objectPath().path()/fullPath;
        }
    }
    else
    {
        fullPath = searchableSurface::objectPath();
    }

    if (!mkDir(fullPath.path()))
    {
        return false;
    }

    triSurface::write(fullPath);

    if (!isFile(fullPath))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
