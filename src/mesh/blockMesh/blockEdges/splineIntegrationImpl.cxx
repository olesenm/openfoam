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

#ifndef Foam_splineIntegration_H
#define Foam_splineIntegration_H

#include "label.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace splineIntegration
{

namespace detail
{

template<class Function>
inline scalar gaussKronrod15
(
    const Function& f,
    const scalar a,
    const scalar b,
    scalar& error
)
{
    static const FixedList<scalar, 8> xgk
    ({
        0.9914553711208126392069,
        0.9491079123427585245262,
        0.8648644233597690727897,
        0.7415311855993944398639,
        0.5860872354676911302941,
        0.4058451513773971669066,
        0.2077849550078984676007,
        0
    });

    static const FixedList<scalar, 8> wgk
    ({
        0.02293532201052922496373,
        0.06309209262997855329070,
        0.1047900103222501838399,
        0.1406532597155259187452,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.2044329400752988924142,
        0.2094821410847278280130
    });

    static const FixedList<scalar, 4> wg
    ({
        0.1294849661688696932706,
        0.2797053914892766679015,
        0.3818300505051189449504,
        0.4179591836734693877551
    });

    const scalar centre = 0.5*(a + b);
    const scalar halfLen = 0.5*(b - a);

    const scalar fc = f(centre);
    scalar kronrodSum = wgk[7]*fc;
    scalar gaussSum = wg[3]*fc;

    for (int i = 0; i < 7; ++i)
    {
        const scalar abscissa = halfLen*xgk[i];
        const scalar fsum = f(centre - abscissa) + f(centre + abscissa);

        kronrodSum += wgk[i]*fsum;

        if (i == 1)
        {
            gaussSum += wg[0]*fsum;
        }
        else if (i == 3)
        {
            gaussSum += wg[1]*fsum;
        }
        else if (i == 5)
        {
            gaussSum += wg[2]*fsum;
        }
    }

    const scalar kronrodLen = halfLen*kronrodSum;
    const scalar gaussLen = halfLen*gaussSum;

    error = Foam::mag(kronrodLen - gaussLen);

    return kronrodLen;
}


template<class Function>
inline scalar adaptiveGaussKronrod15
(
    const Function& f,
    const scalar a,
    const scalar b,
    const scalar absTol,
    const scalar relTol,
    const label depth,
    const label maxDepth
)
{
    scalar error = 0;
    const scalar len = gaussKronrod15(f, a, b, error);
    const scalar tolerance = Foam::max(absTol, relTol*Foam::mag(len));

    if (error <= tolerance || depth >= maxDepth)
    {
        return len;
    }

    const scalar mid = 0.5*(a + b);

    if (mid <= a || mid >= b)
    {
        return len;
    }

    return
        adaptiveGaussKronrod15
        (
            f,
            a,
            mid,
            0.5*absTol,
            relTol,
            depth + 1,
            maxDepth
        )
      + adaptiveGaussKronrod15
        (
            f,
            mid,
            b,
            0.5*absTol,
            relTol,
            depth + 1,
            maxDepth
        );
}

} // End namespace detail


template<class Function>
scalar adaptiveGaussKronrod15(const Function& f)
{
    const scalar absTol = Foam::max(SMALL, scalar(1e-12));
    const scalar relTol = Foam::max(SMALL, scalar(1e-8));
    constexpr label maxDepth = 12;

    return detail::adaptiveGaussKronrod15
    (
        f,
        scalar(0),
        scalar(1),
        absTol,
        relTol,
        0,
        maxDepth
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace splineIntegration
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
