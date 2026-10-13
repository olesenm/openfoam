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

Class
    Foam::multiSolidBodyMotionSolver

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    multiSolidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_multiSolidBodyMotionSolver_H
#define Foam_multiSolidBodyMotionSolver_H

#include "points0MotionSolver.H"
#include "solidBodyMotionFunction.H"
#include "zoneMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiSolidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class multiSolidBodyMotionSolver
:
    public points0MotionSolver
{
    // Private Data

        //- The solid-body motion control functions
        PtrList<solidBodyMotionFunction> motionControls_;

        //- Points to move, per motion function and cellZone selection
        PtrList<zoneMotion> zoneMotions_;


public:

    //- Runtime type information
    TypeName("multiSolidBodyMotionSolver");


    // Generated Methods

        //- No copy construct
        multiSolidBodyMotionSolver(const multiSolidBodyMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const multiSolidBodyMotionSolver&) = delete;


    // Constructors

        //- Construct from mesh and dictionary
        multiSolidBodyMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~multiSolidBodyMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Update local data after mesh redistribution
        virtual void distribute(const mapDistributePolyMesh& map);

        //- Solve for motion
        virtual void solve()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
