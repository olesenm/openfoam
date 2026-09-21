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

        const word refCellName = field.name() + "RefCell";
        const word refPointName = field.name() + "RefPoint";
        const word refValueName = field.name() + "RefValue";

        const entry* refCellFnd = dict.findEntry(refCellName);
        const entry* refPointFnd = dict.findEntry(refPointName);

        label n = 0;
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
            if (Pstream::master(mesh.comm()))
            {
                refCellFnd->stream() >> refCelli;

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
        else
        {
            const point refPt(refPointFnd->stream());

            // Try fast approximate search avoiding octree construction
            refCelli = mesh.findCell(refPt, polyMesh::FACE_PLANES);

            const label communicator = field.mesh().comm();

            label ranking = UPstream::find_first((refCelli >= 0), communicator);
            if (ranking < 0)
            {
                // Nobody found it - try with cell tet-decompositon
                refCelli = field.mesh().findCell(refPt);
                ranking = UPstream::find_first((refCelli >= 0), communicator);
            }

            // Lowest rank takes the lead
            if (ranking != UPstream::myProcNo(communicator))
            {
                refCelli = -1;
            }
            // But also possible that nobody found it
            if (ranking < 0)
            {
                // Nobody found it - try nearest

                const auto& globalCells =
                    mesh.globalData().globalMeshCellAddr();

                // Find the nearest cell to the reference point
                pointIndexHit nearInfo =
                    mesh.cellTree().findNearest(refPt, Foam::sqr(GREAT));

                // Convert the cell to global numbering
                if (nearInfo.hit())
                {
                    nearInfo.setIndex(globalCells.toGlobal(nearInfo.index()));
                }

                // Reduce to find the globally nearest cell
                Foam::reduce
                (
                    nearInfo,
                    [&refPt](pointIndexHit& a, const pointIndexHit& b) -> auto
                    {
                        if (!a.hit())
                        {
                            if (b.hit())
                            {
                                return b;
                            }
                        }
                        else if
                        (
                            b.hit()
                        && magSqr(refPt - a.point()) > magSqr(refPt - b.point())
                        )
                        {
                            return b;
                        }
                        return a;
                    },
                    UPstream::msgType(),
                    mesh.comm()
                );

                if (!nearInfo.hit())
                {
                    // Give up
                    FatalIOErrorInFunction(dict)
                        << "Reference point " << refPointName
                        << " not found in any cell on this domain" << nl
                        << exit(FatalIOError);
                }
                else
                {
                    IOWarningInFunction(dict)
                        << "Reference point " << refPointName
                        << " not found in any cell on any domain" << nl
                        << "    Reverting to nearest cell at "
                        << nearInfo.point() << nl;

                    // Set ref cell on the domain that holds it
                    if (globalCells.isLocal(nearInfo.index()))
                    {
                        refCelli = globalCells.toLocal(nearInfo.index());
                    }
                    else
                    {
                        refCelli = -1;
                    }
                }
            }
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
    #ifdef FULLDEBUG
    if (refCelli >= field.mesh().nCells())
    {
        FatalErrorInFunction
            << "Illegal reference cellID " << refCelli
            << ". Mesh has " << field.mesh().nCells() << ". cells."
            << exit(FatalError);
    }
    #endif
    scalar refCellValue = (refCelli >= 0 ? field[refCelli] : 0.0);

    // Currently distributing the value to all processors. This is generally
    // not needed since only the processor holding the reference cell needs
    // it. Tdb.
    return returnReduce(refCellValue, sumOp<scalar>());
}


// ************************************************************************* //
