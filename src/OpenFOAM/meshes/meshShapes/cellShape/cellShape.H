/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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
    Foam::cellShape

Description
    An analytical geometric cellShape.

    The optional collapse functionality changes the cellModel to the
    correct type after removing any duplicate points.

SourceFiles
    cellShapeI.H
    cellShape.cxx
    cellShapeIO.cxx
    cellShapeEqual.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellShape_H
#define Foam_cellShape_H

#include "pointField.H"
#include "labelList.H"
#include "cellModel.H"
#include "autoPtr.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cell;
class cellShape;
class hexCell;
class tetCell;

Istream& operator>>(Istream& is, cellShape& s);
Ostream& operator<<(Ostream& os, const cellShape& s);

template<>
Ostream& operator<<(Ostream& os, const InfoProxy<cellShape>& ip);

// Common list types
typedef List<cellShape> cellShapeList;  //!< List of cellShape

/*---------------------------------------------------------------------------*\
                          Class cellShape Declaration
\*---------------------------------------------------------------------------*/

class cellShape
:
    public labelList
{
    // Private Data

        //- Access to the cellShape's model
        const cellModel *m;


public:

    // Constructors

        //- Default construct. Empty shape, no cell model.
        inline constexpr cellShape() noexcept;

        //- Move construct from components
        inline cellShape
        (
            const cellModel& model,
            labelList&& labels,
            const bool doCollapse = false
        );

        //- Copy construct from components
        inline cellShape
        (
            const cellModel& model,
            const labelUList& labels,
            const bool doCollapse = false
        );

        //- Copy construct from components
        template<unsigned N>
        inline cellShape
        (
            const cellModel& model,
            const FixedList<label, N>& labels,
            const bool doCollapse = false
        );

        //- Move construct from components, lookup cellModel by type
        inline cellShape
        (
            const cellModel::modelType model,
            labelList&& labels,
            const bool doCollapse = false
        );

        //- Copy construct from components, lookup cellModel by type
        inline cellShape
        (
            const cellModel::modelType model,
            const labelUList& labels,
            const bool doCollapse = false
        );

        //- Copy construct from components, lookup cellModel by type
        template<unsigned N>
        inline cellShape
        (
            const cellModel::modelType model,
            const FixedList<label, N>& labels,
            const bool doCollapse = false
        );

        //- Copy construct from components, lookup cellModel by name
        inline cellShape
        (
            const word& modelName,
            const labelUList& labels,
            const bool doCollapse = false
        );

        //- Copy construct from hex cell labels
        explicit cellShape(const hexCell& labels, bool doCollapse = false);

        //- Copy construct from tet cell labels
        explicit cellShape(const tetCell& labels, bool doCollapse = false);

        //- Construct from Istream
        inline explicit cellShape(Istream& is);

        //- Clone
        autoPtr<cellShape> clone() const
        {
            return autoPtr<cellShape>::New(*this);
        }


    // Member Functions

        //- Model reference
        inline const cellModel& model() const;

        //- Number of points
        inline label nPoints() const noexcept;

        //- Number of edges
        inline label nEdges() const noexcept;

        //- Number of faces
        inline label nFaces() const noexcept;

        //- The points corresponding to this shape
        inline pointField points(const UList<point>& meshPoints) const;

        //- Mesh face labels of this cell (in order of model)
        inline labelList meshFaces
        (
            const UList<Foam::face>& allFaces,
            const cell& cFaces
        ) const;

        //- Mesh edge labels of this cell (in order of model)
        inline labelList meshEdges
        (
            const UList<Foam::edge>& allEdges,
            const labelUList& cEdges
        ) const;

        //- The face for the specified model face
        inline Foam::face face(const label modelFacei) const;

        //- Faces of this cell
        inline faceList faces() const;

        //- Collapsed faces of this cell
        inline faceList collapsedFaces() const;

        //- The edge for the specified model edge
        inline Foam::edge edge(const label modelEdgei) const;

        //- Edges of this shape
        inline edgeList edges() const;

        //- Centroid of the cell
        inline point centre(const UList<point>& points) const;

        //- Scalar magnitude
        inline scalar mag(const UList<point>& points) const;

        //- Return the smallest vertex label used by the shape
        //- or -1 if the shape is empty.
        inline label min() const;

        //- Return the largest vertex label used by the shape
        //- or -1 if the shape is empty.
        inline label max() const;

        //- Reset from components, moving label contents
        inline void reset
        (
            const cellModel& model,
            labelList&& labels,
            bool doCollapse = false
        );

        //- Reset from components, copying label contents
        inline void reset
        (
            const cellModel& model,
            const labelUList& labels,
            bool doCollapse = false
        );

        //- Reset from components, copying label contents
        template<unsigned N>
        inline void reset
        (
            const cellModel& model,
            const FixedList<label, N>& labels,
            bool doCollapse = false
        );

        //- Reset from components, moving label contents
        inline void reset
        (
            cellModel::modelType model,
            labelList&& labels,
            bool doCollapse = false
        );

        //- Reset from components, copying label contents
        inline void reset
        (
            const cellModel::modelType model,
            const labelUList& labels,
            bool doCollapse = false
        );

        //- Reset from components, copying label contents
        template<unsigned N>
        inline void reset
        (
            const cellModel::modelType model,
            const FixedList<label, N>& labels,
            bool doCollapse = false
        );

        //- Reset from hex cell labels
        void reset(const hexCell& labels, bool doCollapse = false);

        //- Reset from tet cell labels
        void reset(const tetCell& labels, bool doCollapse = false);

        //- Collapse shape to correct one after removing duplicate vertices
        void collapse();

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<cellShape> info() const noexcept { return *this; }


    // Comparison

        //- Compare cellShape vertices
        //  \return
        //  -  0: different
        //  - +1: identical values and order used
        //  - -1: identical values, but in different order
        static int compare(const cellShape& a, const cellShape& b);


    // Member Operators

        //- Increment (offset) vertices by given amount
        inline void operator+=(const label vertexOffset);

        //- Compare for identical values (ignoring order)
        bool operator==(const cellShape& other) const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, cellShape& s);
        friend Ostream& operator<<(Ostream& os, const cellShape& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellShapeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
