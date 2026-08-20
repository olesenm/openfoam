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

#include "boundBoxTree.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Visitor>
void Foam::boundBoxTree::visit
(
    label nodei,
    Visitor& visitor
) const
{
    const node& nd = nodes_[nodei];
    const auto& bb = nd.bb();

    if (nd.is_leaf())
    {
        // Filter on (nd.size() > 0) ?
        // Probably should not have a leaf without content anyhow

        visitor.leaf(bb, nd.slice(indices_));
    }
    else  // is_interior
    {
        if (visitor.interior(bb))
        {
            // Do left first, right second
            visit(nd.left(), visitor);
            visit(nd.right(), visitor);
        }
    }
}


template<class Visitor>
void Foam::boundBoxTree::visit_leaf
(
    label nodei,
    Visitor& visitor
) const
{
    const node& nd = nodes_[nodei];

    if (nd.is_leaf())
    {
        visitor(nd);
    }
    else  // is_interior
    {
        // Do left first, right second
        visit_leaf(nd.left(), visitor);
        visit_leaf(nd.right(), visitor);
    }
}


template<class Predicate>
Foam::label Foam::boundBoxTree::find_first_of
(
    label nodei,
    const Predicate& pred
) const
{
    const node& nd = nodes_[nodei];
    const auto& bb = nd.bb();

    if (!pred(bb))
    {
        return -1;
    }

    if (nd.is_leaf())
    {
        for (auto i : nd.slice(indices_))
        {
            if (pred(leafBoxes_[i]))
            {
                return i;
            }
        }
    }
    else  // is_interior
    {
        // Do left first, right second
        if (label i = find_first_of(nd.left(), pred); i >= 0)
        {
            return i;
        }
        if (label i = find_first_of(nd.right(), pred); i >= 0)
        {
            return i;
        }
    }

    return -1;
}


template<class Predicate>
Foam::label Foam::boundBoxTree::count_if
(
    label nodei,
    const Predicate& pred
) const
{
    const node& nd = nodes_[nodei];
    const auto& bb = nd.bb();

    if (!pred(bb))
    {
        return 0;
    }

    if (nd.is_leaf())
    {
        label count = 0;

        for (auto i : nd.slice(indices_))
        {
            if (pred(leafBoxes_[i]))
            {
                ++count;
            }
        }

        return count;
    }
    else  // is_interior
    {
        return
        (
            count_if(nd.left(),  pred)
          + count_if(nd.right(), pred)
        );
    }
}


template<class Predicate>
void Foam::boundBoxTree::find_all
(
    label nodei,
    const Predicate& pred,
    DynamicList<label>& hits
) const
{
    const node& nd = nodes_[nodei];
    const auto& bb = nd.bb();

    if (!pred(bb))
    {
        return;
    }

    if (nd.is_leaf())
    {
        for (auto i : nd.slice(indices_))
        {
            if (pred(leafBoxes_[i]))
            {
                hits.push_back(i);
            }
        }
    }
    else  // is_interior
    {
        // Do left first, right second
        find_all(nd.left(),  pred, hits);
        find_all(nd.right(), pred, hits);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Predicate>
Foam::label Foam::boundBoxTree::find_first_of(const Predicate& pred) const
{
    if (!nodes_.empty())
    {
        return find_first_of(root(), pred);
    }

    return -1;
}


template<class Predicate>
bool Foam::boundBoxTree::any_of(const Predicate& pred) const
{
    return (find_first_of(pred) >= 0);
}


template<class Predicate>
Foam::label Foam::boundBoxTree::count_if
(
    const Predicate& pred
) const
{
    if (!nodes_.empty())
    {
        return count_if(root(), pred);
    }

    return 0;
}


template<class Predicate>
bool Foam::boundBoxTree::find_all
(
    const Predicate& pred,
    DynamicList<label>& hits
) const
{
    hits.clear();
    if (!nodes_.empty())
    {
        find_all(root(), pred, hits);
    }
    return !hits.empty();
}


template<class Visitor>
void Foam::boundBoxTree::visit(Visitor& visitor) const
{
    if (!nodes_.empty())
    {
        visit(root(), visitor);
    }
}


template<class Visitor>
void Foam::boundBoxTree::visit_leaf
(
    Visitor&& visitor
) const
{
    if (!nodes_.empty())
    {
        visit_leaf(root(), visitor);
    }
}


// ************************************************************************* //
