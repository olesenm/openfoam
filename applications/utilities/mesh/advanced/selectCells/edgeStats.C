/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "edgeStats.H"
#include "Time.H"
#include "polyMesh.H"
#include "Ostream.H"
#include "twoDPointCorrector.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::edgeStats::edgeTol_ = 1e-3;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::edgeStats::getNormalDir
(
    const twoDPointCorrector* correct2DPtr
) const
{
    if (correct2DPtr)
    {
        const vector& normal = correct2DPtr->planeNormal();

        if (mag(normal.x()) > 1-edgeTol_)
        {
            return vector::X;
        }
        else if (mag(normal.y()) > 1-edgeTol_)
        {
            return vector::Y;
        }
        else if (mag(normal.z()) > 1-edgeTol_)
        {
            return vector::Z;
        }
    }

    return direction(3);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeStats::edgeStats(const polyMesh& mesh)
:
    mesh_(mesh),
    normalDir_(3)
{
    IOobject motionObj
    (
        "motionProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    if (motionObj.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Reading " << mesh.time().constant() / "motionProperties"
            << endl << endl;

        IOdictionary motionProperties(motionObj);

        if (motionProperties.get<bool>("twoDMotion"))
        {
            Info<< "Correcting for 2D motion" << endl << endl;

            twoDPointCorrector correct2D(mesh);

            normalDir_ = getNormalDir(&correct2D);
        }
    }
}


// Construct from components
Foam::edgeStats::edgeStats
(
    const polyMesh& mesh,
    const twoDPointCorrector* correct2DPtr
)
:
    mesh_(mesh),
    normalDir_(getNormalDir(correct2DPtr))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::edgeStats::minLen(Ostream& os) const
{
    label nAny(0);
    label nX(0);
    label nY(0);
    label nZ(0);

    scalarMinMax limitsAny(GREAT, -GREAT);
    scalarMinMax limitsX(limitsAny);
    scalarMinMax limitsY(limitsAny);
    scalarMinMax limitsZ(limitsAny);

    const edgeList& edges = mesh_.edges();

    for (const edge& e : edges)
    {
        vector eVec(e.vec(mesh_.points()));

        scalar eMag = mag(eVec);

        eVec /= eMag;

        if (mag(eVec.x()) > 1-edgeTol_)
        {
            limitsX.add(eMag);
            nX++;
        }
        else if (mag(eVec.y()) > 1-edgeTol_)
        {
            limitsY.add(eMag);
            nY++;
        }
        else if (mag(eVec.z()) > 1-edgeTol_)
        {
            limitsZ.add(eMag);
            nZ++;
        }
        else
        {
            limitsAny.add(eMag);
            nAny++;
        }
    }

    // serial information
    os  << "Mesh bounding box:" << treeBoundBox(mesh_.points()) << nl << nl
        << "Mesh edge statistics:" << nl
        << "    x aligned :  number:" << nX
        << "\tminLen:" << limitsX.min() << "\tmaxLen:" << limitsX.max() << nl
        << "    y aligned :  number:" << nY
        << "\tminLen:" << limitsY.min() << "\tmaxLen:" << limitsY.max() << nl
        << "    z aligned :  number:" << nZ
        << "\tminLen:" << limitsZ.min() << "\tmaxLen:" << limitsZ.max() << nl
        << "    other     :  number:" << nAny
        << "\tminLen:" << limitsAny.min()
        << "\tmaxLen:" << limitsAny.max() << nl << endl;

    if (normalDir_ == vector::X)
    {
        return Foam::min
        (
            limitsAny.min(),
            Foam::min(limitsY.min(), limitsZ.min())
        );
    }
    else if (normalDir_ == vector::Y)
    {
        return Foam::min
        (
            limitsAny.min(),
            Foam::min(limitsX.min(), limitsZ.min())
        );
    }
    else if (normalDir_ == vector::Z)
    {
        return Foam::min
        (
            limitsAny.min(),
            Foam::min(limitsX.min(), limitsY.min())
        );
    }
    else
    {
        return Foam::min
        (
            limitsAny.min(),
            Foam::min
            (
                limitsX.min(),
                Foam::min(limitsY.min(), limitsZ.min())
            )
        );
    }
}


// ************************************************************************* //
