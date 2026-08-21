/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "cellShape.H"
#include "hexCell.H"
#include "tetCell.H"
#include "degenerateMatcher.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellShape::cellShape(const hexCell& labels, bool doCollapse)
:
    labelList(labels),
    m(cellModel::ptr(cellModel::HEX))
{
    if (doCollapse)
    {
        collapse();
    }
}


Foam::cellShape::cellShape(const tetCell& labels, bool doCollapse)
:
    labelList(labels),
    m(cellModel::ptr(cellModel::TET))
{
    if (doCollapse)
    {
        collapse();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellShape::reset(const hexCell& labels, bool doCollapse)
{
    static_cast<labelList&>(*this) = labels;
    m = cellModel::ptr(cellModel::HEX);

    if (doCollapse)
    {
        collapse();
    }
}


void Foam::cellShape::reset(const tetCell& labels, bool doCollapse)
{
    static_cast<labelList&>(*this) = labels;
    m = cellModel::ptr(cellModel::TET);

    if (doCollapse)
    {
        collapse();
    }
}


void Foam::cellShape::collapse()
{
    operator=(degenerateMatcher::match(*this));
}


// ************************************************************************* //
