/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenFOAM Foundation
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "smoothTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
   defineTypeNameAndDebug(smoothTriSurfaceMesh, 0);
   addToRunTimeSelectionTable(searchableSurface, smoothTriSurfaceMesh, dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::smoothTriSurfaceMesh::calcFeatureEdges(const scalar featureAngle)
{
    if (featureAngle <= -180)
    {
        return;
    }
    const scalar cosAngle = Foam::cos(Foam::degToRad(featureAngle));


    const triSurface& s = *this;
    const auto& pts = s.points();

    const auto& edgeFaces = s.edgeFaces();
    const edgeList& edges = s.edges();

    isBorderEdge_.resize_fill(s.nEdges(), false);
    isPointOnBorderEdge_.resize_fill(s.nPoints(), false);

    forAll(edgeFaces, edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];

        if (eFaces.size() > 2)
        {
            isBorderEdge_.set(edgei);
        }
        else if (eFaces.size() == 2)
        {
            const vector n0(s[eFaces[0]].unitNormal(pts));
            const vector n1(s[eFaces[1]].unitNormal(pts));

            if ((n0&n1) < cosAngle)
            {
                isBorderEdge_.set(edgei);
            }
        }
    }

    for (label edgei : isBorderEdge_)
    {
        const auto& e = edges[edgei];

        isPointOnBorderEdge_.set(e.first());
        isPointOnBorderEdge_.set(e.second());
    }
}


Foam::vector Foam::smoothTriSurfaceMesh::pointNormal
(
    label startFacei,
    label localPointi
) const
{
    const triSurface& s = *this;
    const auto& pts = s.points();

    if (!isPointOnBorderEdge_.test(localPointi))
    {
        vector pn(Zero);
        for (label facei : s.pointFaces()[localPointi])
        {
            pn += s[facei].unitNormal(pts);
        }
        return pn.normalise();
    }


    // Calculate the local point normal on the face. This routine
    // only gets called if the point is on a border edge so we can
    // walk and always hit a border edge.

    const edgeList& edges = s.edges();

    // Get the two edges on the point
    FixedList<label, 2> edgesOnPoint(-1);
    {
        const auto& fEdges = s.faceEdges()[startFacei];

        forAll(fEdges, i)
        {
            const edge& e = edges[fEdges[i]];
            if (e.otherVertex(localPointi) < 0)
            {
                edgesOnPoint[0] = fEdges[fEdges.rcIndex(i)];
                edgesOnPoint[1] = fEdges[fEdges.fcIndex(i)];
                break;
            }
        }
    }

    // Initialise normal
    vector n(s[startFacei].unitNormal(pts));

    // Walk in both directions
    for (const label startEdgei : edgesOnPoint)
    {
        if (startEdgei < 0)
        {
            WarningInFunction
                << "Missing edge for face " << startFacei
                << " appears to be degenerate" << nl;
            continue;
        }

        label facei = startFacei;
        label edgei = startEdgei;

        while (!isBorderEdge_.test(edgei))
        {
            // Cross edge to next face
            const auto& eFaces = s.edgeFaces()[edgei];
            if (eFaces.size() != 2)
            {
                break;
            }

            for (label efacei : eFaces)
            {
                if (facei != efacei)
                {
                    facei = efacei;
                    break;
                }
            }

            // No progress - stop
            if (facei == startFacei)
            {
                break;
            }

            n += s[facei].unitNormal(pts);

            // Cross face to next edge
            for (label ei : s.faceEdges()[facei])
            {
                if (ei != edgei && edges[ei].otherVertex(localPointi) >= 0)
                {
                    edgei = ei;
                    break;
                }
            }
        }
    }

    return n.normalise();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    scalar featureAngle
)
:
    triSurfaceMesh(io)
{
    calcFeatureEdges(featureAngle);
}


Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, dict)
{
    if (scalar featureAngle; dict.readIfPresent("featureAngle", featureAngle))
    {
        calcFeatureEdges(featureAngle);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoothTriSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    const triSurface& s = *this;
    const auto& pts = s.points();

    DebugPout
        << "smoothTriSurfaceMesh::getNormal :"
        << " surface:" << this->searchableSurface::name()
        << " tris:" << s.size()
        << " feat edges:" << isBorderEdge_.count()
        << endl;

    normal.resize_nocopy(info.size());

    forAll(info, i)
    {
        if (const auto& p = info[i].point(); info[i].hit())
        {
            label facei = info[i].index();

            // Get local coordinates in triangle
            barycentric2D coordinates
            (
                s[facei].tri(pts).pointToBarycentric(p)
            );

            // Average point normals
            const triFace& localTri = s.localFaces()[facei];

            vector n(Zero);
            forAll(localTri, fp)
            {
                n += coordinates[fp]*pointNormal(facei, localTri[fp]);
            }
            normal[i] = n.normalise();
        }
        else
        {
            // Set to what?
            normal[i] = Foam::zero{};
        }
    }
}


Foam::label Foam::smoothTriSurfaceMesh::measure
(
    List<vector>* centres,
    List<scalar>* boundingRadiusSqr,
    List<label>* regions,
    List<vector>* normals
) const
{
    const triSurface& s = *this;

    const label output_size = s.size();

    const label parentCount = triSurfaceMesh::measure
    (
        centres,
        boundingRadiusSqr,
        regions,
        nullptr  // defer normals to below
    );

    bool usedParams = (parentCount > 0);

    // Extract smooth normals.
    //
    // Since the sample point is the triangle centre, just need to equally
    // weight the point normals
    if (normals)
    {
        usedParams = true;
        auto& normal = *normals;
        normal.resize_nocopy(output_size);

        forAll(normal, facei)
        {
            auto& n = normal[facei];
            const triFace& localTri = s.localFaces()[facei];

            n = pointNormal(facei, localTri.a());
            n += pointNormal(facei, localTri.b());
            n += pointNormal(facei, localTri.c());

            n.normalise();  // inplace normalise
        }
    }

    return (usedParams ? output_size : 0);
}


// ************************************************************************* //
