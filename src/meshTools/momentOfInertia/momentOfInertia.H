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
    Foam::momentOfInertia

Description
    Calculates the inertia tensor and principal axes and moments of a
    polyhedra/cells/triSurfaces.  Inertia can either be of the solid body or
    of a thin shell.

    References:
    \verbatim
        "Fast and Accurate Computation of Polyhedral Mass Properties",
        Brian Mirtich,
        Journal of Graphics Tools, Volume 1, Number 2, 1996.
    \endverbatim

    https://people.eecs.berkeley.edu/~jfc/mirtich/massProps.html
    https://people.eecs.berkeley.edu/~jfc/mirtich/papers/volInt.ps.gz
    https://people.eecs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar

SourceFiles
    momentOfInertia.H

\*---------------------------------------------------------------------------*/

#ifndef momentOfInertia_H
#define momentOfInertia_H

#include "triFaceList.H"
#include "triSurface.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class momentOfInertia Declaration
\*---------------------------------------------------------------------------*/

class momentOfInertia
{

public:

    static void massPropertiesSolid
    (
        const pointField& pts,
        const triFaceList& triFaces,
        scalar density,
        scalar& mass,
        vector& cM,
        symmTensor& J
    );

    //- Tensor overload of the correct symmTensor overload
    //
    //  \deprecated(2026-11) - use symmTensor overload
    FOAM_DEPRECATED_FOR(2026-11, "symmTensor overload")
    static void massPropertiesSolid
    (
        const pointField& pts,
        const triFaceList& triFaces,
        scalar density,
        scalar& mass,
        vector& cM,
        tensor& J
    )
    {
        symmTensor symmJ(Foam::zero{});
        massPropertiesSolid(pts, triFaces, density, mass, cM, symmJ);
        J = tensor(symmJ);
    }

    static void massPropertiesShell
    (
        const pointField& pts,
        const triFaceList& triFaces,
        scalar density,
        scalar& mass,
        vector& cM,
        symmTensor& J,
        bool doReduce = false
    );

    //- Tensor overload of the correct symmTensor overload
    //
    //  \deprecated(2026-11) - use symmTensor overload
    FOAM_DEPRECATED_FOR(2026-11, "symmTensor overload")
    static void massPropertiesShell
    (
        const pointField& pts,
        const triFaceList& triFaces,
        scalar density,
        scalar& mass,
        vector& cM,
        tensor& J,
        bool doReduce = false
    )
    {
        symmTensor symmJ(Foam::zero{});
        massPropertiesShell
        (
            pts,
            triFaces,
            density,
            mass,
            cM,
            symmJ,
            doReduce
        );
        J = tensor(symmJ);
    }

    static void massPropertiesSolid
    (
        const triSurface& surf,
        scalar density,
        scalar& mass,
        vector& cM,
        symmTensor& J
    );

    //- Tensor overload of the correct symmTensor overload
    //
    //  \deprecated(2026-11) - use symmTensor overload
    FOAM_DEPRECATED_FOR(2026-11, "symmTensor overload")
    static void massPropertiesSolid
    (
        const triSurface& surf,
        scalar density,
        scalar& mass,
        vector& cM,
        tensor& J
    )
    {
        symmTensor symmJ(Foam::zero{});
        massPropertiesSolid(surf, density, mass, cM, symmJ);
        J = tensor(symmJ);
    }

    static void massPropertiesShell
    (
        const triSurface& surf,
        scalar density,
        scalar& mass,
        vector& cM,
        symmTensor& J,
        bool doReduce = false
    );

    //- Tensor overload of the correct symmTensor overload
    //
    //  \deprecated(2026-11) - use symmTensor overload
    FOAM_DEPRECATED_FOR(2026-11, "symmTensor overload")
    static void massPropertiesShell
    (
        const triSurface& surf,
        scalar density,
        scalar& mass,
        vector& cM,
        tensor& J,
        bool doReduce = false
    )
    {
        symmTensor symmJ(Foam::zero{});
        massPropertiesShell(surf, density, mass, cM, symmJ, doReduce);
        J = tensor(symmJ);
    }

    static void massPropertiesPatch
    (
        const polyPatch& pp,
        scalar density,
        scalar& mass,
        vector& cM,
        symmTensor& J,
        bool doReduce = false
    );

    //- Tensor overload of the correct symmTensor overload
    //
    //  \deprecated(2026-11) - use symmTensor overload
    FOAM_DEPRECATED_FOR(2026-11, "symmTensor overload")
    static void massPropertiesPatch
    (
        const polyPatch& pp,
        scalar density,
        scalar& mass,
        vector& cM,
        tensor& J,
        bool doReduce = false
    )
    {
        symmTensor symmJ(Foam::zero{});
        massPropertiesPatch(pp, density, mass, cM, symmJ, doReduce);
        J = tensor(symmJ);
    }

    static symmTensor applyParallelAxisTheorem
    (
        scalar mass,
        const vector& cM,
        const symmTensor& J,
        const vector& refPt
    );

    //- Tensor overload of the correct symmTensor overload
    //
    //  \deprecated(2026-11) - use symmTensor overload
    FOAM_DEPRECATED_FOR(2026-11, "symmTensor overload")
    static tensor applyParallelAxisTheorem
    (
        scalar mass,
        const vector& cM,
        const tensor& J,
        const vector& refPt
    )
    {
        const vector d(refPt - cM);
        return J + mass*((d & d)*I - d*d);
    }

    static void meshInertia
    (
        const polyMesh& mesh,
        List<symmTensor>& result
    );

    static void meshInertia
    (
        const polyMesh& mesh,
        List<tensor>& result
    );

    template<class TensorType = symmTensor>
    static tmp<Field<TensorType>> meshInertia
    (
        const polyMesh& mesh
    )
    {
        auto tresult = tmp<Field<TensorType>>::New();
        meshInertia(mesh, tresult.ref());
        return tresult;
    }

    static symmTensor meshInertia
    (
        const polyMesh& mesh,
        label celli
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
