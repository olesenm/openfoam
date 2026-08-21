/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::treeDecomp

Description
    Geometrical domain decomposition using a binary kd-tree.

    The default uses an axis-aligned kd-tree,
    but an optional coordinate system transform can also be specified.

    Method coefficients:
    \table
        Property  | Description                             | Required | Default
        rotation  | cartesian coordinate rotation           | no  |
        processorWeights | Relative target load per domain  | no  | uniform
    \endtable

    At each recursion level, the current group of points is split along the
    coordinate direction with the largest span. The split index is choosen
    to yield a balanced count on each side of the split.

    When processor weights are specified, the bin split and the point split
    fraction are based on the relative weight of the domains on each side.

Note
  - This decomposition method is primarily useful for distributing
    surface meshes. For volume meshes, the processor connectivity will
    not be good at all.

SourceFiles
    treeDecomp.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_treeDecomp_H
#define Foam_treeDecomp_H

#include "decompositionMethod.H"
#include "cartesianCS.H"

namespace Foam
{

// Forward Declarations
class bitSet;

/*---------------------------------------------------------------------------*\
                         Class treeDecomp Declaration
\*---------------------------------------------------------------------------*/

class treeDecomp
:
    public decompositionMethod
{
    // Private Data

        //- Local coordinate rotation (optional)
        coordSystem::cartesian csys_;

        //- Relative target load per domain
        scalarList procWeights_;


    // Private Member Functions

        //- Recursively split the points addressed by the (index) span
        //- among the domains addressed by bins, always splitting along
        //- the coordinate direction with the largest span.
        //
        // Without weights (empty list), the point-count split is
        // proportional to the number of domains assigned to each side so
        // that odd splits are handled correctly (eg, 9 -> 4/9 front,
        // 5/9 back).
        //
        // With weights (indexed by domain id, ie, by the values addressed
        // by bins), the bin split and the point-count split are instead
        // based on the relative (summed) weight assigned to each side.
        // A side with zero summed weight may legitimately receive zero
        // points.
        static void splitPoints
        (
            const UList<point>& points,
            labelUList& indices,  //!< [in,out] current indices into points
            const labelRange& span,  //!< span of indices to split
            const labelRange& bins,  //!< range of target domain ids
            //! relative weight per domain id or empty (ie, uniform)
            const UList<scalar>& procWeights,
            labelUList& finalDecomp  //!< [out] the decomposition
        );

        //- Decompose on the master process only.
        //- Points already gathered, rotated etc
        static labelList treeDecompose_impl
        (
            const label numDomains,
            const UList<point>& points,
            //! relative weight per domain id or empty (ie, uniform)
            const UList<scalar>& procWeights = UList<scalar>::null()
        );

        //- Read input values and set the rotation matrix
        void readCoeffs(const dictionary& coeffs);

        //- Decompose the point field by gathering points (in parallel),
        //- applying any coordinate transformations and finally calling
        //- treeDecompose_impl() to do the real work.
        labelList treeDecompose(const pointField& points) const;

public:

    // Generated Methods

        //- No copy construct
        treeDecomp(const treeDecomp&) = delete;

        //- No copy assignment
        void operator=(const treeDecomp&) = delete;


    //- Runtime type information
    TypeName("tree");


    // Constructors

        //- Construct with number of domains,
        //- no coefficients or constraints
        explicit treeDecomp(label numDomains);

        //- Construct with number of domains and specified domain weighting
        treeDecomp(label numDomains, const labelUList& weights);

        //- Construct with number of domains and specified domain selection
        treeDecomp(label numDomains, const boolUList& selection);

        //- Construct with number of domains and specified domain selection
        treeDecomp(label numDomains, const bitSet& selection);

        //- Construct with number of domains and coordinate transformation,
        //- no coefficients or constraints
        treeDecomp(label numDomains, const coordinateSystem& csys);

        //- Construct for decomposition dictionary and optional region name
        explicit treeDecomp
        (
            const dictionary& decompDict,
            const word& regionName = "",
            int select = selectionType::DEFAULT
        );


    //- Destructor
    virtual ~treeDecomp() = default;


    // Member Functions

        //- The target load per domain (non-normalized)
        const UList<scalar>& weights() const noexcept { return procWeights_; }

        //- Purely geometric methods
        virtual bool geometric() const { return true; }

        //- Sends all points to the master for decomposition
        virtual bool parallelAware() const { return true; }


    // Member Functions

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const pointField& points,
            //! Weighting is ignored (ie, uniform)
            const scalarField& pointWeights = scalarField::null()
        ) const
        {
            return treeDecompose(points);
        }

        //- Decompose with uniform weights
        virtual labelList decompose
        (
            const polyMesh& mesh,  //!< unused
            const pointField& points,
            //! Weighting is ignored (ie, uniform)
            const scalarField& pointWeights = scalarField::null()
        ) const
        {
            return treeDecompose(points);
        }

        //- Decompose with uniform weights (connectivity ignored)
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,  //!< unused
            const pointField& cc,
            //! Weighting is ignored (ie, uniform)
            const scalarField& cWeights = scalarField::null()
        ) const
        {
            return treeDecompose(cc);
        }

        //- Decompose with uniform weights (connectivity ignored)
        virtual labelList decompose
        (
            const labelListList& globalCellCells,  //!< unused
            const pointField& cc,
            //! Weighting is ignored (ie, uniform)
            const scalarField& cWeights = scalarField::null()
        ) const
        {
            return treeDecompose(cc);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
