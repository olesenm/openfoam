/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "findRefCell.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::setRefCell
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    if (fieldRef.needReference() || forceReference)
    {
        const auto& mesh = field.mesh();
        const auto communicator = mesh.comm();

        const word refCellName = field.name() + "RefCell";
        const word refPointName = field.name() + "RefPoint";
        const word refValueName = field.name() + "RefValue";

        const entry* refCellFnd = dict.findEntry(refCellName);
        const entry* refPointFnd = dict.findEntry(refPointName);

        int n = 0;
        if (refCellFnd) n++;
        if (refPointFnd) n++;
        if (n != 1)
        {
            FatalIOErrorInFunction(dict)
                << "Unable to set reference cell for field " << field.name()
                << nl
                << "    Please supply either " << refPointName
                << " or " << refCellName << nl << exit(FatalIOError);
        }

        if (refCellFnd)
        {
            refCellFnd->stream() >> refCelli;

            if (UPstream::master(communicator))
            {
                if (refCelli < 0 || refCelli >= mesh.nCells())
                {
                    FatalIOErrorInFunction(dict)
                        << "Illegal master cellID " << refCelli
                        << ". Should be 0.." << mesh.nCells()
                        << nl
                        << "    Preferentially use "
                        << refPointName << " to supply a location instead" << nl
                        << exit(FatalIOError);
                }
            }
            else
            {
                refCelli = -1;
            }
        }
        else if (refPointFnd)
        {
            const point refPt(refPointFnd->stream());

            // Try fast approximate search avoiding octree construction
            refCelli = mesh.findCell(refPt, polyMesh::FACE_PLANES);

            const auto myProci = UPstream::myProcNo(communicator);
            const auto numProc = UPstream::nProcs(communicator);

            int bestProci = UPstream::find_first(refCelli >= 0, communicator);
            if (bestProci < 0)
            {
                // Nobody found it - try with cell tet-decompositon
                refCelli = mesh.findCell(refPt);
                bestProci = UPstream::find_first(refCelli >= 0, communicator);
            }

            // The refCell is only valid on the best domain
            if (bestProci != myProci)
            {
                refCelli = -1;
            }

            // But also possible that nobody found it
            if (bestProci < 0)
            {
                // Try nearest cell to the reference point
                Tuple2<scalar, int> nearest(VGREAT, UPstream::max_ranks());

                pointIndexHit pHit =
                    mesh.cellTree().findNearest(refPt, Foam::sqr(GREAT));

                if (const auto& p = pHit.point(); pHit.hit())
                {
                    refCelli = pHit.index();
                    nearest = { refPt.distSqr(p), myProci };
                }
                else
                {
                    refCelli = -1;
                }

                bestProci = UPstream::minLocation(nearest, communicator);

                // The refCell is only valid on the best domain
                if (myProci != bestProci)
                {
                    refCelli = -1;
                }

                if (bestProci < 0 || bestProci >= numProc)
                {
                    // Give up
                    FatalIOErrorInFunction(dict)
                        << "Reference point " << refPointName
                        << " : " << refPt
                        << " not found in any cell on this domain" << nl
                        << exit(FatalIOError);
                }
                else
                {
                    IOWarningInFunction(dict)
                        << "Reference point " << refPointName
                        << " not found in any cell on any domain" << nl
                        << "    Reverting to nearest cell, which is on rank="
                        << bestProci << nl;
                }
            }
        }
        else
        {
            // Programming/logic error...
            FatalErrorInFunction << abort(FatalError);
        }

        dict.readEntry(refValueName, refValue);

        return true;
    }
    else
    {
        refCelli = -1;
    }

    return false;
}


bool Foam::setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    return setRefCell(field, field, dict, refCelli, refValue, forceReference);
}


Foam::scalar Foam::getRefCellValue
(
    const volScalarField& field,
    const label refCelli
)
{
    const label nCells = field.mesh().nCells();

    #ifdef FULLDEBUG
    if (refCelli >= nCells)
    {
        FatalErrorInFunction
            << "Illegal reference cellID " << refCelli
            << ". Mesh has " << nCells << ". cells."
            << exit(FatalError);
    }
    #endif

    scalar value = (refCelli >= 0 && refCelli < nCells) ? field[refCelli] : 0;

    // Currently distributing the value to all processors. This is generally
    // not needed since only the processor holding the reference cell needs
    // it. Tdb.
    return returnReduce(value, sumOp<scalar>());
}


// ************************************************************************* //
