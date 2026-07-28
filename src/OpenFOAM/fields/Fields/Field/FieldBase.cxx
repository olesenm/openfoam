/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "Field.H"
#include "debug.H"
#include "dictionary.H"
#include "error.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::FieldBase::typeName("Field");

bool Foam::FieldBase::allowConstructFromLargerSize = false;

bool Foam::FieldBase::unifiedGeometricField
(
    Foam::debug::optimisationSwitch("unifiedGeometricField", 0)
);
registerOptSwitch
(
    "unifiedGeometricField",
    bool,
    Foam::FieldBase::unifiedGeometricField
);


int Foam::FieldBase::localBoundaryConsistency_
(
    Foam::debug::optimisationSwitch("localBoundaryConsistency", 1)
);
registerOptSwitch
(
    "localConsistency",
    int,
    Foam::FieldBase::localBoundaryConsistency_
);


Foam::scalar Foam::FieldBase::localBoundaryTolerance_
(
    Foam::debug::floatOptimisationSwitch
    (
        "localBoundaryConsistency::tolerance", 0
    )
);
registerOptSwitch
(
    "localBoundaryConsistency::tolerance",
    Foam::scalar,
    Foam::FieldBase::localBoundaryTolerance_
);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// This is a really ugly solution, but no obvious simpler method

void Foam::FieldBase::warnLocalBoundaryConsistencyCompat
(
    const dictionary& dict
)
{
    // New: localBoundaryConsistency
    // Old:
    // - localConsistency
    // - point(Vector)Field::Boundary::localConsistency
    // - point(Spherical|Symm)?TensorField::Boundary::localConsistency
    // - (surface|vol)(Scalar|Vector)Field::Boundary::localConsistency
    // - (surface|vol)(Spherical|Symm)?TensorField::Boundary::localConsistency

    // New: localBoundaryConsistency::tolerance
    // Old:
    // - (vol|surface)(Scalar|Vector)Field::Boundary::tolerance
    // - (vol|surface)(Spherical|Symm)?TensorField::Boundary::tolerance

    {
        constexpr int version(2412);

        const word flagName("localBoundaryConsistency");
        const word tolName("localBoundaryConsistency::tolerance");

        // Warn: "using " kw -> flagName;
        const auto emitWarning =
            [=](const std::string& oldName, const std::string& newName)
            {
                if (error::master())
                {
                    std::cerr
                        << "--> FOAM IOWarning :" << nl
                        << "    Found [v" << version << "] '"
                        << oldName.c_str() << "' entry instead of '"
                        << newName.c_str() << "'" << nl
                        << std::endl;
                    error::warnAboutAge("keyword", version);
                }
            };

        for (const entry& e : dict)
        {
            const auto& kw = e.keyword();

            if (kw.contains("Field::Boundary::"))
            {
                if (kw.ends_with("Field::Boundary::localConsistency"))
                {
                    emitWarning(kw, flagName);
                }
                else if (kw.ends_with("Field::Boundary::tolerance"))
                {
                    emitWarning(kw, tolName);
                }
            }
            else if (kw == "localConsistency")
            {
                emitWarning(kw, flagName);
            }
        }
    }
}


// ************************************************************************* //
