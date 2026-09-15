/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::domainDecompositionDryRun

Description
    Testing of domain decomposition for finite-volume meshes

SourceFiles
    domainDecompositionDryRun.cxx
    domainDecompositionDryRunWrite.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_domainDecompositionDryRun_H
#define Foam_domainDecompositionDryRun_H

#include "fvMesh.H"

namespace Foam
{

// Forward Declarations
template<class T> class CompactListList;

/*---------------------------------------------------------------------------*\
                  Class domainDecompositionDryRun Declaration
\*---------------------------------------------------------------------------*/

class domainDecompositionDryRun
{
    // Private Data

        //- The mesh
        fvMesh mesh_;

        //- Optional non-standard file for decomposeParDict
        const fileName decompDictFile_;

        //- Optional override for numberOfSubdomains in decomposeParDict
        const label nDomainsOverride_;

        //- Optional override for method in decomposeParDict
        const word methodNameOverride_;


    // Private Members

        //- Write decomposition as volScalarField for visualization
        void writeVolField
        (
            const word& timeName,
            const labelUList& cellToProc
        ) const;

        //- Write cell distribution as VTU file (binary)
        void writeVTK
        (
            const fileName& file,
            const labelUList& cellToProc
        ) const;

        //- Write originating mesh as scotch ".grf" and ".xyz" files and
        //- the target decomposition as a ".map"
        void writeScotch
        (
            const fileName& file,
            const CompactListList<label>& cellCells,
            const labelUList& cellToProc
        ) const;


public:

    // Constructors

        //- Construct from components.
        domainDecompositionDryRun
        (
            //! The IOobject for mesh
            const IOobject& io,
            //! optional non-standard location for the decomposeParDict file
            const fileName& decompDictFile = "",
            //! optional override value for numberOfSubdomains
            const label nDomains = 0,
            //! optional override value for decomposition method
            const word& methodName = ""
        );


    //- Destructor
    ~domainDecompositionDryRun() = default;


    // Member Functions

        //- The mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Perform dry-run
        void execute
        (
            //! Write VTK (optionally scotch output, when verbose > 1)
            const bool writeCellDist,
            //! With verbose = 1, report extra statistics.
            //! With verbose > 1, write scotch files when writeCellDist is true
            int verbose = 0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
