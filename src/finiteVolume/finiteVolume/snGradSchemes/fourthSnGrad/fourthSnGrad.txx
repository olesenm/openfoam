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

\*---------------------------------------------------------------------------*/

#include "fourthSnGrad.H"
#include "correctedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<SurfaceField<Type>> fourthSnGrad<Type>::correction
(
    const VolumeField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    auto tcorr = tmp<SurfaceField<Type>>::New
    (
        IOobject
        (
            "snGradCorr(" + vf.name() + ')',
            vf.instance(),
            mesh,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE
        ),
        mesh,
        vf.dimensions()*mesh.deltaCoeffs().dimensions()
    );
    auto& corr = tcorr.ref();
    corr.setOriented();

    // Work component-wise to support all tensor ranks.
    // For each scalar component, compute:
    //   (1/6)*(corrected_snGrad - face-normal Gauss gradient)
    //   + corrected non-orthogonal correction

    correctedSnGrad<scalar> corrScheme(mesh);

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        tmp<volScalarField> tvsf = vf.component(cmpt);
        const volScalarField& vsf = tvsf.cref();

        corr.replace
        (
            cmpt,
            (1.0/6.0)
           *(
                corrScheme.snGrad(vsf)
              - (
                    mesh.Sf()
                  & linearInterpolate(gaussGrad<scalar>(mesh).grad(vsf))
                )
               /mesh.magSf()
            )
          + corrScheme.correction(vsf)
        );
    }

    return tcorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
