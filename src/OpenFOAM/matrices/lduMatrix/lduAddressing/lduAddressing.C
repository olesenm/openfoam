/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "lduAddressing.H"
#include "scalarField.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduAddressing::calcPatchSort() const
{
    if (patchSortAddr_.size() == nPatches())
    {
        FatalErrorIn("lduAddressing::calcPatchSort() const")
            << "patch sort already calculated"
            << abort(FatalError);
    }

    patchSortAddr_.setSize(nPatches());
    patchSortCells_.setSize(nPatches());

    for(label i=0; i<nPatches(); i++)
    {
        if (!patchAvailable(i))
        {
            continue;
        }

        const labelList& nbr = patchAddr(i);
        labelList* sortPtr_ = new labelList(nbr.size(), -1);
        patchSortAddr_.set(i, sortPtr_);

        labelList& lst = *sortPtr_;

        labelList nbrTmp(nbr);

        const labelRange indices(0, label(nbr.size()));

        apply_policy(static_cast<std::size_t>(nbr.size()), [&](auto policy)
        {
            std::copy
            (
                policy,
                indices.begin(),
                indices.end(),
                lst.begin()
            );
        });

        apply_policy(static_cast<std::size_t>(lst.size()), [&](auto policy)
        {
            std::stable_sort
            (
                policy,
                lst.begin(),
                lst.end(),
                [&](label a, label b)
                {
                    return nbrTmp[a] < nbrTmp[b];
                }
            );
        });

        labelList* cellsSortPtr = new labelList(nbr.size());
        patchSortCells_.set(i, cellsSortPtr);
        labelList& cellsSort = *cellsSortPtr;

        apply_policy(static_cast<std::size_t>(lst.size()), [&](auto policy)
        {
            std::transform
            (
                policy,
                lst.begin(),
                lst.end(),
                cellsSort.begin(),
                [&](label idx)
                {
                    return nbr[idx];
                }
            );
        });

        auto new_end = std::unique(cellsSort.begin(), cellsSort.end());
        cellsSort.resize(std::distance(cellsSort.begin(), new_end));
    }
}


void Foam::lduAddressing::calcPatchSortStart() const
{
    if (patchSortStartAddr_.size() == nPatches())
    {
        // Already calculated
        return;
    }

    patchSortStartAddr_.setSize(nPatches());

    for (label i=0; i<nPatches(); i++)
    {
        if (!patchAvailable(i))
        {
            continue;
        }

        const labelList& nbr = patchAddr(i);
        // assuming this returns std::vector<label>
        const labelList& lsrt = patchSortAddr(i);

        // 1. Permutation: fill nbrSort with nbr[idx] based on lsrt
        std::vector<label> nbrSort(nbr.size());

        apply_policy(static_cast<std::size_t>(lsrt.size()), [&](auto policy)
        {
            std::transform
            (
                policy,
                lsrt.begin(),
                lsrt.end(),
                nbrSort.begin(),
                [&](label idx)
                {
                    return nbr[idx];
                }
            );
        });

        std::vector<size_t> tmpSum;
        if (!nbrSort.empty()) 
        {
            size_t count = 1;
            for (size_t i=1; i<nbrSort.size(); i++)
            {
                if (nbrSort[i] == nbrSort[i-1])
                {
                    ++count;
                }
                else
                {
                    tmpSum.push_back(count);
                    count = 1;
                }
            }
            tmpSum.push_back(count); // last group
        }
  
        // (nbrsize+1)
        labelList* patchSortStartAddr = new labelList(tmpSum.size() + 1);
        patchSortStartAddr_.set(i, patchSortStartAddr);
        labelList& lsrtStart = *patchSortStartAddr;

        apply_policy(static_cast<std::size_t>(tmpSum.size()), [&](auto policy)
        {
            std::exclusive_scan
            (
                policy,
                tmpSum.begin(),
                tmpSum.end(),
                lsrtStart.begin(),
                label(0)
            );
        });

        lsrtStart[tmpSum.size()] = lsrt.size(); // last element is total size
        // lsrtStart now contains starting indices
        // (exclusive scan of group sizes)
    }
}

void Foam::lduAddressing::calcLosort() const
{
    if (losortPtr_)
    {
        FatalErrorInFunction
            << "losort already calculated"
            << abort(FatalError);
    }

    // Scan the neighbour list to find out how many times the cell
    // appears as a neighbour of the face. Done this way to avoid guessing
    // and resizing list
    labelList nNbrOfFace(size(), Foam::zero{});

    const labelUList& nbr = upperAddr();

    forAll(nbr, nbrI)
    {
        nNbrOfFace[nbr[nbrI]]++;
    }

    // Create temporary neighbour addressing
    labelListList cellNbrFaces(size());

    forAll(cellNbrFaces, celli)
    {
        cellNbrFaces[celli].setSize(nNbrOfFace[celli]);
    }

    // Reset the list of number of neighbours to zero
    nNbrOfFace = 0;

    // Scatter the neighbour faces
    forAll(nbr, nbrI)
    {
        cellNbrFaces[nbr[nbrI]][nNbrOfFace[nbr[nbrI]]] = nbrI;

        nNbrOfFace[nbr[nbrI]]++;
    }

    // Gather the neighbours into the losort array
    losortPtr_ = std::make_unique<labelList>(nbr.size(), -1);
    auto& lst = *losortPtr_;

    // Set counter for losort
    label lstI = 0;

    forAll(cellNbrFaces, celli)
    {
        const labelUList& curNbr = cellNbrFaces[celli];

        forAll(curNbr, curNbrI)
        {
            lst[lstI] = curNbr[curNbrI];
            lstI++;
        }
    }
}


void Foam::lduAddressing::calcOwnerStart() const
{
    if (ownerStartPtr_)
    {
        FatalErrorInFunction
            << "owner start already calculated"
            << abort(FatalError);
    }

    const labelUList& own = lowerAddr();

    ownerStartPtr_ = std::make_unique<labelList>(size() + 1, own.size());
    auto& ownStart = *ownerStartPtr_;

    // Set up first lookup by hand
    ownStart[0] = 0;
    label nOwnStart = 0;
    label i = 1;

    forAll(own, facei)
    {
        label curOwn = own[facei];

        if (curOwn > nOwnStart)
        {
            while (i <= curOwn)
            {
                ownStart[i++] = facei;
            }

            nOwnStart = curOwn;
        }
    }
}


void Foam::lduAddressing::calcLosortStart() const
{
    if (losortStartPtr_)
    {
        FatalErrorInFunction
            << "losort start already calculated"
            << abort(FatalError);
    }

    const labelUList& nbr = upperAddr();
    losortStartPtr_ = std::make_unique<labelList>(size() + 1, nbr.size());
    auto& lsrtStart = *losortStartPtr_;

    const labelUList& lsrt = losortAddr();

    // Set up first lookup by hand
    lsrtStart[0] = 0;
    label nLsrtStart = 0;
    label i = 0;

    forAll(lsrt, facei)
    {
        // Get neighbour
        const label curNbr = nbr[lsrt[facei]];

        if (curNbr > nLsrtStart)
        {
            while (i <= curNbr)
            {
                lsrtStart[i++] = facei;
            }

            nLsrtStart = curNbr;
        }
    }

    // Set up last lookup by hand
    lsrtStart[size()] = nbr.size();
}


void Foam::lduAddressing::calcLoCSR() const
{
    if (lowerCSRAddrPtr_)
    {
        FatalErrorInFunction
            << "lowerCSRAddr already calculated"
            << abort(FatalError);
    }

    lowerCSRAddrPtr_ = std::make_unique<labelList>(lowerAddr().size());
    map(lowerAddr(), *lowerCSRAddrPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::lduAddressing::losortAddr() const
{
    if (!losortPtr_)
    {
        calcLosort();
    }

    return *losortPtr_;
}


const Foam::labelUList& Foam::lduAddressing::ownerStartAddr() const
{
    if (!ownerStartPtr_)
    {
        calcOwnerStart();
    }

    return *ownerStartPtr_;
}


const Foam::labelUList& Foam::lduAddressing::losortStartAddr() const
{
    if (!losortStartPtr_)
    {
        calcLosortStart();
    }

    return *losortStartPtr_;
}

const Foam::labelUList& Foam::lduAddressing::patchSortCells(const label i) const
{
    if (patchSortCells_.size() != nPatches())
    {
        calcPatchSort();
    }
    return patchSortCells_[i];
}


const Foam::labelUList&
Foam::lduAddressing::patchSortAddr(const label i) const
{
    if (patchSortAddr_.size() != nPatches())
    {
        calcPatchSort();
    }

    return patchSortAddr_[i];
}


const Foam::labelUList&
Foam::lduAddressing::patchSortStartAddr(const label i) const
{
    if (patchSortStartAddr_.size() != nPatches())
    {
        calcPatchSortStart();
    }

    return patchSortStartAddr_[i];
}


const Foam::labelUList& Foam::lduAddressing::lowerCSRAddr() const
{
    if (!lowerCSRAddrPtr_)
    {
        calcLoCSR();
    }

    return *lowerCSRAddrPtr_;
}


void Foam::lduAddressing::clearOut()
{
    losortPtr_.reset(nullptr);
    ownerStartPtr_.reset(nullptr);
    losortStartPtr_.reset(nullptr);
    lowerCSRAddrPtr_.reset(nullptr);
    patchSortCells_.clear();
    patchSortAddr_.clear();
    patchSortStartAddr_.clear();
}


Foam::label Foam::lduAddressing::triIndex(const label a, const label b) const
{
    label own = min(a, b);

    label nbr = max(a, b);

    label startLabel = ownerStartAddr()[own];

    label endLabel = ownerStartAddr()[own + 1];

    const labelUList& neighbour = upperAddr();

    for (label i=startLabel; i<endLabel; i++)
    {
        if (neighbour[i] == nbr)
        {
            return i;
        }
    }

    // If neighbour has not been found, something has gone seriously
    // wrong with the addressing mechanism
    FatalErrorInFunction
        << "neighbour " << nbr << " not found for owner " << own << ". "
        << "Problem with addressing"
        << abort(FatalError);

    return -1;
}


Foam::Tuple2<Foam::label, Foam::scalar> Foam::lduAddressing::band() const
{
    const labelUList& owner = lowerAddr();
    const labelUList& neighbour = upperAddr();

    labelList cellBandwidth(size(), Foam::zero{});

    forAll(neighbour, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // Note: mag not necessary for correct (upper-triangular) ordering.
        label diff = nei-own;
        cellBandwidth[nei] = max(cellBandwidth[nei], diff);
    }

    label bandwidth = max(cellBandwidth);

    // Do not use field algebra because of conversion label to scalar
    scalar profile = 0.0;
    forAll(cellBandwidth, celli)
    {
        profile += 1.0*cellBandwidth[celli];
    }

    return Tuple2<label, scalar>(bandwidth, profile);
}


// ************************************************************************* //
