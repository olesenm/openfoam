/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenFOAM Foundation
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::smoothTriSurfaceMesh

Description
    Variant of triSurfaceMesh that interpolates face-normals to obtain
    point normals.

    Note: when constructing from dictionary has optional parameters:
        - scale     : scaling factor.
        - tolerance : relative tolerance for doing intersections
                      (see triangle::intersection)
        - minQuality: discard triangles with low quality when getting normal
        - featureAngle: limit interpolation to current 'region' only.
                        Region limited by feature edges.

SourceFiles
    smoothTriSurfaceMesh.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_smoothTriSurfaceMesh_H
#define Foam_smoothTriSurfaceMesh_H

#include "triSurfaceMesh.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class smoothTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class smoothTriSurfaceMesh
:
    public triSurfaceMesh
{
    // Private Data

        //- Edges with discontinuity
        bitSet isBorderEdge_;

        //- Points on border edges
        bitSet isPointOnBorderEdge_;


    // Private Member Functions

        //- Mark feature edges
        void calcFeatureEdges(const scalar featureAngle);

        //- Calculate local point normal
        vector pointNormal(label startFacei, label localPointi) const;


public:

    // Generated Methods

        //- No copy construct
        smoothTriSurfaceMesh(const smoothTriSurfaceMesh&) = delete;

        //- No copy assignment
        void operator=(const smoothTriSurfaceMesh&) = delete;


    //- Runtime type information
    TypeName("smoothTriSurfaceMesh");


    // Constructors

        //- Construct read, specify feature angle
        smoothTriSurfaceMesh(const IOobject& io, scalar featureAngle);

        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        smoothTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~smoothTriSurfaceMesh() = default;


    // Member Functions

        //- From a set of points and indices get the normal
        virtual void getNormal
        (
            const List<pointIndexHit>&,
            vectorField& normal
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
