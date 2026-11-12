/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::searchableCylinder

Description
    Searching on a cylinder.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required | Default
        type        | cylinder                          | selector |
        point1      | coordinate of endpoint            | yes   |
        point2      | coordinate of endpoint            | yes   |
        radius      | cylinder radius                   | yes   |
    \endtable

Note
    Longer type name : \c searchableCylinder

SourceFiles
    searchableCylinder.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_searchableCylinder_H
#define Foam_searchableCylinder_H

#include "treeBoundBox.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class searchableCylinder Declaration
\*---------------------------------------------------------------------------*/

class searchableCylinder
:
    public searchableSurface
{
private:

    // Private Member Data

        //- The 'left' point
        const point point1_;

        //- The 'right' point
        const point point2_;

        //- Length of vector point2-point1
        const scalar magDir_;

        //- Normalised vector point2-point1
        const vector unitDir_;

        //- The radius
        const scalar radius_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Find nearest point on cylinder.
        void findNearest_impl
        (
            const point& sample,
            const scalar nearestDistSqr,
            pointIndexHit& info
        ) const;

        scalar radius2(const point& pt) const;

        //- Find intersection with cylinder
        void findLineAll_impl
        (
            const point& start,
            const point& end,
            pointIndexHit& near,
            pointIndexHit& far
        ) const;

        //- Return the boundBox of the cylinder
        boundBox calcBounds() const;

public:

    // Generated Methods

        //- No copy construct
        searchableCylinder(const searchableCylinder&) = delete;

        //- No copy assignment
        void operator=(const searchableCylinder&) = delete;


    //- Runtime type information
    TypeName("searchableCylinder");


    // Constructors

        //- Construct from components
        searchableCylinder
        (
            const IOobject& io,
            const point& point1,
            const point& point2,
            const scalar radius
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableCylinder
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableCylinder() = default;


    // Member Functions

        //- Names of regions
        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return true;
        }

        //- What is type of points outside bounds
        virtual volumeType outsideVolumeType() const
        {
            return volumeType::OUTSIDE;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- The number of the points that define the surface.
        virtual label nPoints() const
        {
            return 2;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const
        {
            NotImplemented;
            return false;
        }


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point.
            //  Unknown if cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField& points,
                List<volumeType>& volType
            ) const;

            //- WIP.
            //- Query for surface information.
            //- \returns the number of surface elements queried.
            //- \returns 0 for an empty surface or if all query parameters
            //- were nullptr.
            virtual label measure
            (
                //- [out] surface element centres
                List<vector>* centres,
                //- [out] surface element bounding sphere radius(squared)
                List<scalar>* boundingRadiusSqr,
                //- [out] surface element region ids
                List<label>* regions,
                //- [out] surface element unit normals
                List<vector>* normals
            ) const;



        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
