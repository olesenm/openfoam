/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "nearestToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearestToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, nearestToCell, word);
    addToRunTimeSelectionTable(topoSetSource, nearestToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, nearestToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, nearestToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        nearestToCell,
        word,
        nearest
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        nearestToCell,
        istream,
        nearest
    );
}


Foam::topoSetSource::addToUsageTable Foam::nearestToCell::usage_
(
    nearestToCell::typeName,
    "\n    Usage: nearestToCell (pt0 .. ptn)\n\n"
    "    Select the nearest cell for each of the points pt0 ..ptn\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearestToCell::combine(topoSet& set, const bool add) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);

    // Nearest cell index (locally)
    labelList objectIndices(points_.size(), Foam::zero{});

    // The nearest selection key. Construct as missed sample
    // == (large distance, proc = max_ranks)
    List<Tuple2<scalar, int>> nearest
    (
        points_.size(),
        { VGREAT, UPstream::max_ranks() }
    );

    forAll(points_, pointi)
    {
        const auto& sample = points_[pointi];

        if
        (
            label celli = mesh_.findNearestCell(sample);
            celli >= 0
        )
        {
            const point& cc = mesh_.cellCentres()[celli];

            objectIndices[pointi] = celli;
            nearest[pointi] = { sample.distSqr(cc), myProci };
        }
    }

    // Find rank with nearest - globally consistent
    if (UPstream::vendorCode() == UPstream::vendorCodes::intelmpi)
    {
        // Avoid slow MINLOC on intelmpi (OCT-2026)
        Pstream::listCombineReduce
        (
            nearest, minFirstEqOp{}, UPstream::msgType(), communicator
        );
    }
    else
    {
        UPstream::minLocation(nearest, communicator);
    }

    forAll(nearest, samplei)
    {
        if (myProci == nearest[samplei].second())
        {
            label celli = objectIndices[samplei];

            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearestToCell::nearestToCell
(
    const polyMesh& mesh,
    const pointField& points
)
:
    topoSetCellSource(mesh),
    points_(points)
{}


Foam::nearestToCell::nearestToCell
(
    const polyMesh& mesh,
    pointField&& points
)
:
    topoSetCellSource(mesh),
    points_(std::move(points))
{}


Foam::nearestToCell::nearestToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    points_(dict.get<pointField>("points"))
{}


Foam::nearestToCell::nearestToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    points_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearestToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells nearest to " << points_ << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells nearest to " << points_ << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
