/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tanuj Ravi
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

Note
    Some ideas taken from https://github.com/tanujravi/numpyToFoam

\*---------------------------------------------------------------------------*/

#include "foamNumpyCore.H"
#include <sstream>

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

// Header magic (6 chars)
constexpr char npy_magic[] = "\x93NUMPY";

// The size of the magic + version + header_len (version 1.0)
constexpr std::size_t preamble_1_0 = 10;

// The size of the magic + version + header_len (version 2.0)
// constexpr std::size_t preamble_2_0 = 12;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

std::string Foam::fileFormats::numpyCore::headerDict
(
    const std::string& dtype,
    const std::string& shape,
    bool fortranOrder
)
{
    std::string header;
    {
        std::ostringstream dict;
        dict
            << "{'descr': '" << dtype
            << "', 'fortran_order': " << (fortranOrder ? "True" : "False")
            << ", 'shape': " << shape << '}';

        header = dict.str();
    }

    const std::size_t total = (preamble_1_0 + header.size() + 1);

    // Require 16 byte alignment
    constexpr std::size_t alignment = 16;

    if (auto remainer = (total % alignment); remainer)
    {
        // Pad with spaces
        header.append(alignment - remainer, ' ');
    }

    header.push_back('\n');

    return header;
}


std::size_t Foam::fileFormats::numpyCore::writeHeader
(
    std::ostream& os,
    const std::string& dtype,
    const std::string& shape,
    bool fortranOrder
)
{
    const std::string header
    (
        headerDict(dtype, shape, fortranOrder)
    );

    os.write(npy_magic, 6);
    os.put(char(0x01));  // major
    os.put(char(0x00));  // minor

    // version 1.0 : Header length as little endian uint16
    {
        auto hlen = static_cast<std::uint16_t>(header.size());
        os.put(static_cast<char>(hlen & 0xFF));
        os.put(static_cast<char>((hlen >> 8) & 0xFF));
    }

    // version 2.0 : Header length as little endian uint32
    // {
    //     auto hlen = static_cast<std::uint32_t>(header.size());
    //     os.put(static_cast<char>(hlen & 0xFF));
    //     os.put(static_cast<char>((hlen >> 8) & 0xFF));
    //     os.put(static_cast<char>((hlen >> 16) & 0xFF));
    //     os.put(static_cast<char>((hlen >> 24) & 0xFF));
    // }

    os.write(header.data(), header.size());

    return (preamble_1_0 + header.size());
}


// ************************************************************************* //
