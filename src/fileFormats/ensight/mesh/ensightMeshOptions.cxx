/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "ensightMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::options::options() noexcept
:
    lazy_(false),
    internal_(true),
    boundary_(true),
    cellZones_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightMesh::options::useFaceZones() const noexcept
{
    return !faceZoneInclude_.empty();
}


void Foam::ensightMesh::options::reset()
{
    internal_ = true;
    boundary_ = true;
    cellZones_ = true;
    patchInclude_.clear();
    patchExclude_.clear();
    cellZoneInclude_.clear();
    cellZoneExclude_.clear();
    faceZoneInclude_.clear();
    faceZoneExclude_.clear();
}


bool Foam::ensightMesh::options::lazy(bool on) noexcept
{
    bool old(lazy_);
    lazy_ = on;
    return old;
}


bool Foam::ensightMesh::options::useInternalMesh(bool on) noexcept
{
    bool old(internal_);
    internal_ = on;
    return old;
}


bool Foam::ensightMesh::options::useBoundaryMesh(bool on)
{
    bool old(boundary_);
    boundary_ = on;

    if (!boundary_ && anyPatchFilters())
    {
        patchInclude_.clear();
        patchExclude_.clear();

        WarningInFunction
            << "Deactivating boundary, removed old patch selection"
            << endl;
    }

    return old;
}


bool Foam::ensightMesh::options::useCellZones(bool on)
{
    bool old(cellZones_);
    cellZones_ = on;

    if (!cellZones_ && anyCellZoneFilters())
    {
        cellZoneInclude_.clear();
        cellZoneExclude_.clear();

        WarningInFunction
            << "Deactivating cellZones, removed old zone selection"
            << endl;
    }

    return old;
}


void Foam::ensightMesh::options::patchInclude
(
    const UList<wordRe>& patterns
)
{
    if (boundary_)
    {
        patchInclude_ = wordRes(patterns);
    }
    else
    {
        patchInclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, boundary is disabled" << endl;
        }
    }
}


void Foam::ensightMesh::options::patchInclude
(
    List<wordRe>&& patterns
)
{
    if (boundary_)
    {
        patchInclude_ = wordRes(std::move(patterns));
    }
    else
    {
        patchInclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, boundary is disabled" << endl;
        }
    }
}


void Foam::ensightMesh::options::patchExclude
(
    const UList<wordRe>& patterns
)
{
    if (boundary_)
    {
        patchExclude_ = wordRes(patterns);
    }
    else
    {
        patchExclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, boundary is disabled" << endl;
        }
    }
}


void Foam::ensightMesh::options::patchExclude
(
    List<wordRe>&& patterns
)
{
    if (boundary_)
    {
        patchExclude_ = wordRes(std::move(patterns));
    }
    else
    {
        patchExclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, boundary is disabled" << endl;
        }
    }
}


void Foam::ensightMesh::options::cellZoneInclude
(
    const UList<wordRe>& patterns
)
{
    if (cellZones_)
    {
        cellZoneInclude_ = wordRes(patterns);
    }
    else
    {
        cellZoneInclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, cell zones are disabled" << endl;
        }
    }
}


void Foam::ensightMesh::options::cellZoneInclude
(
    List<wordRe>&& patterns
)
{
    if (cellZones_)
    {
        cellZoneInclude_ = wordRes(std::move(patterns));
    }
    else
    {
        cellZoneInclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, cell zones are disabled" << endl;
        }
    }
}


void Foam::ensightMesh::options::cellZoneExclude
(
    const UList<wordRe>& patterns
)
{
    if (cellZones_)
    {
        cellZoneExclude_ = wordRes(patterns);
    }
    else
    {
        cellZoneExclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, cell zones are disabled" << endl;
        }
    }
}


void Foam::ensightMesh::options::cellZoneExclude
(
    List<wordRe>&& patterns
)
{
    if (cellZones_)
    {
        cellZoneExclude_ = wordRes(std::move(patterns));
    }
    else
    {
        cellZoneExclude_.clear();

        if (!patterns.empty())
        {
            WarningInFunction
                << "Ignoring selection, cell zones are disabled"
                << endl;
        }
    }
}


void Foam::ensightMesh::options::faceZoneInclude
(
    const UList<wordRe>& patterns
)
{
    faceZoneInclude_ = wordRes(patterns);
}


void Foam::ensightMesh::options::faceZoneInclude
(
    List<wordRe>&& patterns
)
{
    faceZoneInclude_ = wordRes(std::move(patterns));
}


void Foam::ensightMesh::options::faceZoneExclude
(
    const UList<wordRe>& patterns
)
{
    faceZoneExclude_ = wordRes(patterns);
}


void Foam::ensightMesh::options::faceZoneExclude
(
    List<wordRe>&& patterns
)
{
    faceZoneExclude_ = wordRes(std::move(patterns));
}


void Foam::ensightMesh::options::print(Ostream& os) const
{
    os << "internal: " << Switch::name(internal_) << nl;

    os << "cellZones: " << Switch::name(cellZones_) << nl;
    if (cellZones_ && anyCellZoneFilters())
    {
        os.incrIndent();
        if (!cellZoneInclude_.empty())
        {
            os.writeKeyword("include") << flatOutput(cellZoneInclude_) << nl;
        }
        if (!cellZoneExclude_.empty())
        {
            os.writeKeyword("exclude") << flatOutput(cellZoneExclude_) << nl;
        }
        os.decrIndent();
    }

    os << "boundary: " << Switch::name(boundary_) << nl;
    if (boundary_ && anyPatchFilters())
    {
        os.incrIndent();
        if (!patchInclude_.empty())
        {
            os.writeKeyword("include") << flatOutput(patchInclude_) << nl;
        }
        if (!patchExclude_.empty())
        {
            os.writeKeyword("exclude") << flatOutput(patchExclude_) << nl;
        }
        os.decrIndent();
    }

    os << "faceZones: " << Switch::name(useFaceZones()) << nl;
    if (useFaceZones() && anyFaceZoneFilters())
    {
        os.incrIndent();
        if (!faceZoneInclude_.empty())
        {
            os.writeKeyword("include") << flatOutput(faceZoneInclude_) << nl;
        }
        if (!faceZoneExclude_.empty())
        {
            os.writeKeyword("exclude") << flatOutput(faceZoneExclude_) << nl;
        }
        os.decrIndent();
    }
}


// ************************************************************************* //
