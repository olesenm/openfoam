/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2026 OpenCFD Ltd.
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

#include "foamVtuWriter.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::vtuWriter::writeUniform
(
    const word& fieldName,
    const Type& val
)
{
    label nValues(0);

    // These are local counts - the backend does the rest
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
        nValues = cellSlab_.size();
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
        nValues = pointSlab_.size();
    }
    else
    {
        reportBadState
        (
            FatalErrorInFunction,
            outputState::CELL_DATA,
            outputState::POINT_DATA
        )
            << " for uniform field " << fieldName << nl << endl
            << exit(FatalError);

        return;
    }

    vtk::fileWriter::writeUniform<Type>(fieldName, val, nValues);
}


template<class Type>
void Foam::vtk::vtuWriter::writeCellData
(
    const word& fieldName,
    const UList<Type>& field
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::CELL_DATA)
            << " for field " << fieldName << nl << endl
            << exit(FatalError);
    }

    const labelUList& cellMap = vtuCells_.cellMap();

    // Needs adjustment??
    // // Use the cellMap?
    // const bool useMap
    // (
    //     vtuCells_.useCellMap()
    //  || (cellMap.size() != cellSlab_.size())
    //  || (cellMap.size() != field.size())
    // );

    this->beginDataArray<Type>(fieldName, nTotalCells());

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), field, cellMap);
    }
    else
    {
        vtk::writeList(format(), field, cellMap);
    }

    this->endDataArray();
}


template<class Type>
void Foam::vtk::vtuWriter::writePointData
(
    const word& fieldName,
    const UList<Type>& field
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << fieldName << nl << endl
            << exit(FatalError);
    }

    this->beginDataArray<Type>(fieldName, nTotalPoints());

    if (parallel_)
    {
        if (returnReduceOr(vtuCells_.merged()))
        {
            // With pointMap...
            const UIndirectList<Type> values(field, vtuCells_.pointMap());

            vtk::writeListParallel(format_.ref(), values);
        }
        else
        {
            // No pointMap...
            vtk::writeListParallel(format_.ref(), field);
        }
    }
    else
    {
        // No pointMap for non-parallel...
        vtk::writeList(format(), field);
    }

    this->endDataArray();
}


template<class Type>
void Foam::vtk::vtuWriter::writePointData
(
    const word& fieldName,
    const UList<Type>& field,
    const UList<Type>& extra
)
{
    if
    (
        isNull(extra)
     || (parallel_ ? returnReduceAnd(extra.empty()) : extra.empty())
    )
    {
        // Can write without the "extra" data...
        writePointData(fieldName, field);
        return;
    }

    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << fieldName << nl << endl
            << exit(FatalError);
    }

    this->beginDataArray<Type>(fieldName, nTotalPoints());

    if (parallel_)
    {
        if (returnReduceOr(vtuCells_.merged()))
        {
            // With pointMap...
            const UIndirectList<Type> values(field, vtuCells_.pointMap());

            vtk::writeListsParallel(format_.ref(), values, extra);
        }
        else
        {
            // No pointMap...
            vtk::writeListsParallel(format_.ref(), field, extra);
        }
    }
    else
    {
        // No pointMap for non-parallel...
        vtk::writeLists(format(), field, extra);
    }

    this->endDataArray();
}


// ************************************************************************* //
