/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "geomDecomp.H"
#include "bitSet.H"
#include "specifiedRotation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::geomDecomp::setOrder()
{
    const word order(coeffsDict_.getOrDefault<word>("order", ""));

    if (order.empty())
    {
        return;
    }
    else if (order.size() != 3)
    {
        FatalIOErrorInFunction(decompDict_)
            << "Number of characters in order (" << order << ") != 3"
            << exit(FatalIOError);
    }

    for (int i = 0; i < 3; ++i)
    {
        // Change [x-z] -> [0-2]

        switch (order[i])
        {
            case 'x': order_[i] = 0; break;
            case 'y': order_[i] = 1; break;
            case 'z': order_[i] = 2; break;

            default:
                FatalIOErrorInFunction(decompDict_)
                    << "Illegal decomposition order " << order << nl
                    << "It should only contain x, y or z"
                    << exit(FatalIOError);
                break;
        }
    }
}


void Foam::geomDecomp::readCoeffs()
{
    coeffsDict_.readIfPresent("delta", delta_);

    coeffsDict_.readEntry("n", n_);

    if (nDomains_ != n_.x()*n_.y()*n_.z())
    {
        // Verify that the input makes sense
        FatalIOErrorInFunction(coeffsDict_)
            << "Wrong number of domain divisions in geomDecomp:" << nl
            << "Number of domains    : " << nDomains_ << nl
            << "Wanted decomposition : " << n_
            << exit(FatalIOError);
    }
    setOrder();

    if
    (
        coeffsDict_.readIfPresent("processorWeights", procWeights_)
     && (procWeights_.size() > 1)
    )
    {
        if (procWeights_.size() < nDomains_)
        {
            FatalIOErrorInFunction(coeffsDict_)
                << "processorWeights (" << procWeights_.size()
                << ") less than numberOfSubdomains (" << nDomains_
                << ")\n"
                << exit(FatalIOError);
        }
        else if (procWeights_.size() > nDomains_)
        {
            WarningInFunction
                << "processorWeights (" << procWeights_.size()
                << ") greater than numberOfSubdomains (" << nDomains_
                << ") ... truncating" << endl;

            procWeights_.resize(nDomains_, 0);
        }

        sanitizeWeights(nDomains_, procWeights_);
    }
    else
    {
        // Also treat size 1 like uniform weighted
        procWeights_.clear();
    }

    // Optional (cartesian) coordinate system transform
    auto csysPtr = coordinateSystem::NewIfPresent(coeffsDict_, "transform");

    if (csysPtr)
    {
        csys_ = *csysPtr;
    }
    else if (equal(delta_, 0))
    {
        csys_.clear();  // Reset to identity
    }
    else
    {
        const scalar d = 1 - 0.5*delta_*delta_;
        const scalar d2 = sqr(d);

        const scalar a = delta_;
        const scalar a2 = sqr(a);

        // Direction (forward/reverse) doesn't matter much
        tensor rot
        (
            d2,         -a*d,         a,
            a*d - a2*d,  a*a2 + d2,  -2*a*d,
            a*d2 + a2,   a*d - a2*d,  d2 - a2
        );

        // origin=0
        csys_ = coordinateSystem(coordinateRotations::specified(rot));
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::geomDecomp::adjustPoints
(
    const pointField& points
) const
{
    return csys_.localPosition(points);
}


void Foam::geomDecomp::checkDecompositionDirections
(
    const Vector<label>& meshDirs
) const
{
    for (direction dir = 0; dir < Vector<label>::nComponents; ++dir)
    {
        if (n_[dir] > 1 && meshDirs[dir] == -1)
        {
            WarningInFunction
                << "Trying to decompose a 1/2D mesh"
                << " into " << n_[dir]
                << " parts in direction "
                << Vector<label>::componentNames[dir]
                << endl;
        }
    }
}


Foam::scalarList Foam::geomDecomp::binWeightFractions
(
    const label procOffset,
    const label stride,
    const label nBins
) const
{
    // procOffset: first proc ID in current subtree
    // stride: spacing between proc IDs while advancing current decomp digit
    // nBins: num of child bins created at this recursion level

    if (stride <= 0 || nBins <= 0)
    {
        return scalarList();
    }
    else if (procWeights_.size() != nDomains_)
    {
        return scalarList(nBins, scalar(1));
    }

    scalarList fractions(nBins, Zero);

    // group: Number of final processors represented by each current bin
    const label group = nDomains_/(stride*nBins);

    scalar total = 0;
    for (label b = 0; b < nBins; ++b)
    {
        scalar w = 0;
        for (label m = 0; m < group; ++m)
        {
            const label proci = procOffset + (b + m*nBins)*stride;
            w += procWeights_[proci];
        }
        fractions[b] = w;
        total += w;
    }

    if (total > SMALL)
    {
        // Normalize the fraction
        for (auto& fraction : fractions)
        {
            fraction /= total;
        }
    }
    else
    {
        fractions = scalar(1)/scalar(nBins);
    }

    return fractions;
}


Foam::scalarList Foam::geomDecomp::axisWeightFractions
(
    const direction cmpt
) const
{
    const label nBins = n_[cmpt];

    if (nBins <= 0)
    {
        return scalarList();
    }
    else if (procWeights_.size() != nDomains_)
    {
        return scalarList(nBins, scalar(1));
    }

    scalarList fractions(nBins, Zero);
    scalar total = 0;

    forAll(procWeights_, proci)
    {
        const label x = proci % n_.x();
        const label y = (proci/n_.x()) % n_.y();
        const label z = proci/(n_.x()*n_.y());

        label bin = x;
        if (cmpt == vector::Y)
        {
            bin = y;
        }
        else if (cmpt == vector::Z)
        {
            bin = z;
        }

        fractions[bin] += procWeights_[proci];
        total += procWeights_[proci];
    }

    if (total > SMALL)
    {
        // Normalize the fraction
        for (auto& fraction : fractions)
        {
            fraction /= total;
        }
    }
    else
    {
        fractions = scalar(1)/scalar(nBins);
    }

    return fractions;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geomDecomp::geomDecomp(const Vector<label>& divisions)
:
    decompositionMethod(divisions.x()*divisions.y()*divisions.z()),
    delta_(0.001),
    csys_(),
    n_(divisions),
    order_(0,1,2),
    coeffsDict_(dictionary::null)
{}


Foam::geomDecomp::geomDecomp
(
    const word& derivedType,
    const dictionary& decompDict,
    int select
)
:
    decompositionMethod(decompDict),
    delta_(0.001),
    csys_(),
    n_(1,1,1),
    order_(0,1,2),
    coeffsDict_(findCoeffsDict(derivedType + "Coeffs", select))
{
    readCoeffs();
}


Foam::geomDecomp::geomDecomp
(
    const word& derivedType,
    const dictionary& decompDict,
    const word& regionName,
    int select
)
:
    decompositionMethod(decompDict, regionName),
    delta_(0.001),
    csys_(),
    n_(1,1,1),
    order_(0,1,2),
    coeffsDict_(findCoeffsDict(derivedType + "Coeffs", select))
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::geomDecomp::setWeights(const labelUList& weights)
{
    procWeights_.clear();
    copyWeights(nDomains_, weights, procWeights_);
}


void Foam::geomDecomp::setWeights(const bitSet& selection)
{
    procWeights_.clear();
    selectWeights(nDomains_, selection, procWeights_);
}


void Foam::geomDecomp::setWeights(const boolUList& selection)
{
    procWeights_.clear();
    selectWeights(nDomains_, selection, procWeights_);
}


// ************************************************************************* //
