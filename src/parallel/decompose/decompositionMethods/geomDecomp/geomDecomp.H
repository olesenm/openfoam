/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::geomDecomp

Description
    Base for geometrical domain decomposition methods

    Base coefficients:
    \table
        Property  | Description                             | Required | Default
        n         | (nx ny nz)                              | yes |
        order     | order of operation                      | no  | xyz
        delta     | delta (jitter) for rotation matrix      | no  | 0.001
        rotation  | coordinate rotation                     | no  |
        processorWeights | Relative target load per domain  | no  | uniform
    \endtable

Note
  - For \c processorWeights, only the ratios matter since \c processorWeights
    is internally normalised; it only provides best-effort approximation.

SourceFiles
    geomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_geomDecomp_H
#define Foam_geomDecomp_H

#include "decompositionMethod.H"
#include "cartesianCS.H"
#include "Vector.H"
#include "tmp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class geomDecomp Declaration
\*---------------------------------------------------------------------------*/

class geomDecomp
:
    public decompositionMethod
{
    // Private Data

        //- Small delta (default: 0.001) to avoid staggering when
        //- mesh itself is aligned with x/y/z
        scalar delta_;

        //- Local coordinates, or delta rotation tensor
        coordSystem::cartesian csys_;


    // Private Member Functions

        //- Convert ordering string ("xyz") into list of components.
        //  Checks for bad entries, but no check for duplicate entries.
        void setOrder();

        //- Read input values and initialize rotation matrix
        void readCoeffs();


protected:

    // Protected Data

        //- The divisions
        Vector<label> n_;

        //- Decomposition order in terms of components (optional)
        Vector<direction> order_;

        //- Coefficients for all derived methods
        const dictionary& coeffsDict_;

        //- Relative target load per domain
        scalarList procWeights_;


    // Protected Member Functions

        //- Apply delta (jitter) or rotation to coordinates
        tmp<pointField> adjustPoints(const pointField&) const;

        //- Check that mesh directions are compatible with decomposition
        void checkDecompositionDirections(const Vector<label>&) const;

        //- Return the normalized processor-weight fraction over the current
        //- processor subtree, not necessarily over all processors
        scalarList binWeightFractions
        (
            const label procOffset,
            const label stride,
            const label nBins
        ) const;

        //- Return normalized processor-weight fractions for the bins along
        //- one geometric axis
        scalarList axisWeightFractions(const direction cmpt) const;


public:

    // Constructors

        //- Construct with number of x/y/z division
        //- (no coefficients or constraints)
        explicit geomDecomp(const Vector<label>& divisions);

        //- Construct for derived type name and decomposition dictionary
        geomDecomp
        (
            const word& derivedType,
            const dictionary& decompDict,
            int select = selectionType::DEFAULT
        );

        //- Construct for derived type name, decomposition dictionary
        //- and region name
        geomDecomp
        (
            const word& derivedType,
            const dictionary& decompDict,
            const word& regionName,
            int select = selectionType::DEFAULT
        );


    // Member Functions

        //- Purely geometric methods
        virtual bool geometric() const { return true; }


    // Domain Weighting

        //- Has processor/domain weights?
        bool hasWeights() const noexcept { return !procWeights_.empty(); }

        //- The target load per domain (non-normalized)
        const UList<scalar>& weights() const noexcept { return procWeights_; }

        //- Remove processor/domain weighting
        void clearWeights() { procWeights_.clear(); }

        //- Set/reset the processor/domain weighting
        void setWeights(const labelUList& weights);

        //- Set/reset the processor/domain selection (weighting)
        void setWeights(const bitSet& selection);

        //- Set/reset the processor/domain selection (weighting)
        void setWeights(const boolUList& selection);


    // Member Functions

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights = scalarField::null()
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
