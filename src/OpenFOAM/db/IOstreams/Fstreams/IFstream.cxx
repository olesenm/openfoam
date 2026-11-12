/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
    Copyright (C) 2026 Mark Olesen
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

#include "IFstream.H"
#include "OSspecific.H"  // For isFile(), fileSize()

#include <algorithm>
#include <limits>
#include <type_traits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IFstream, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace
{

// Regular or chunk-wise reading of file contents into
// a std::string, std::vector<char> or DynamicList<char>.
//
// The reported file size is only used as a size hint.

template<class CharBufferType>
CharBufferType readContents_impl(Foam::IFstream& ifs)
{
    // Report the public API name instead of the template's mangled signature
    [[maybe_unused]]
    constexpr const char* funcName = "IFstream::readContents";

    constexpr bool is_string_buffer =
        std::is_base_of_v<std::string, CharBufferType>;

    constexpr bool is_vector_buffer =
        std::is_base_of_v<std::vector<char>, CharBufferType>;

    // DynamicList if it is not std::string or std::vector
    constexpr bool is_dynlist_buffer = (!is_string_buffer && !is_vector_buffer);

    // Reading in 4kB increments
    constexpr uint64_t base_chunk = uint64_t(4096);

    // Limits when using Foam::List
    [[maybe_unused]]
    constexpr uint64_t maxLabelLimit =
        uint64_t(std::numeric_limits<Foam::label>::max());


    const std::streamsize inputSize = ifs.fileSize();

    CharBufferType buffer;

    if (inputSize > 0)
    {
        auto& iss = ifs.stdStream();

        const bool is_compressed =
            (Foam::IOstreamOption::COMPRESSED == ifs.compression());

        if constexpr (is_dynlist_buffer)
        {
            if (!is_compressed && uint64_t(inputSize) > maxLabelLimit)
            {
                // FatalIOErrorInFunction
                FatalIOErrorIn(funcName, ifs)
                    << "File size " << inputSize
                    << " exceeds the Foam::List<char> "
                    << int(8*sizeof(Foam::label)) << "-bit size limit"
                    << Foam::exit(Foam::FatalIOError);
            }
        }

        // For compressed files, no idea how large the result will be.
        // So read chunk-wise.
        //
        // Using the compressed size for the chunk size:
        // 50% compression = 2 iterations
        // 66% compression = 3 iterations
        // ...
        //
        // For uncompressed files:
        // use (fileSize+1) to detect EOF in a single pass,
        // unless this exceeds the representable size of a Foam::List

        const uint64_t firstSize =
        (
            is_compressed
          ? std::max(base_chunk, 2*uint64_t(inputSize))
          :
            // uncompressed:
            (
                // Foam::List needs a separate EOF probe if the file size
                // is already at the List::max_size()
                (is_dynlist_buffer && uint64_t(inputSize) == maxLabelLimit)
              ? uint64_t(inputSize)
              : (uint64_t(inputSize) + 1u)
            )
        );

        // Subsequent growth is linear (never doubling) to limit over-allocation
        const uint64_t chunkSize =
        (
            is_compressed
          ? firstSize
          : std::max(base_chunk, uint64_t(inputSize)/8)
        );

        uint64_t total = 0;

        // Loop until nothing more can be read
        for (uint64_t want = firstSize; true; want = chunkSize)
        {
            const uint64_t nextTotal = total + want;

            if constexpr (is_dynlist_buffer)
            {
                if (nextTotal > maxLabelLimit)
                {
                    FatalIOErrorIn(funcName, ifs)
                        << "File contents exceed the Foam::List<char> "
                        << int(8*sizeof(Foam::label)) << "-bit size limit"
                        << Foam::exit(Foam::FatalIOError);
                }
            }

            // Exact capacity, then size
            if constexpr (is_string_buffer || is_vector_buffer)
            {
                buffer.reserve(size_t(nextTotal));
                buffer.resize(size_t(nextTotal));
            }
            else
            {
                // Foam::DynamicList
                buffer.setCapacity(Foam::label(nextTotal));
                buffer.resize(Foam::label(nextTotal));
            }

            ifs.readRaw(buffer.data() + total, std::streamsize(want));
            const std::streamsize nread = iss.gcount();

            if
            (
                nread <= 0
             || nread == std::numeric_limits<std::streamsize>::max()
            )
            {
                break;  // Failed or EOF: treat as normal 'done'
            }

            total += uint64_t(nread);

            if constexpr (is_dynlist_buffer)
            {
                if (total == maxLabelLimit)
                {
                    break;
                }
            }

            if (uint64_t(nread) < want)
            {
                break;  // Short read: EOF
            }
        }

        if constexpr (is_dynlist_buffer)
        {
            if (total == maxLabelLimit)
            {
                char eofProbe;
                ifs.readRaw(&eofProbe, 1);

                if (iss.gcount() > 0)
                {
                    // FatalIOErrorInFunction
                    FatalIOErrorIn(funcName, ifs)
                        << "File contents exceed the Foam::List<char> "
                        << int(8*sizeof(Foam::label)) << "-bit size limit"
                        << Foam::exit(Foam::FatalIOError);
                }
            }
        }

        if constexpr (is_dynlist_buffer)
        {
            // DynamicList
            buffer.resize(Foam::label(total));
        }
        else
        {
            buffer.resize(size_t(total));
        }
    }


    // Future (experimental) idea.
    // Trim out excessive capacity, but with padding to a 8-byte word
    // boundary in case we want to quickly resize and send as 8-byte
    // aggregates for large send/recv.

    // constexpr bool shrink_somewhat = true;
    //
    // if constexpr (shrink_somewhat)
    // {
    //     constexpr int64_t smaller_base_chunk = int64_t(base_chunk/2);
    //
    //     const auto realLen = int64_t(buffer.size());
    //
    //     if
    //     (
    //         int64_t(buffer.capacity()) - smaller_base_chunk
    //       > int64_t(realLen)
    //     )
    //     {
    //         int padding = (realLen % 8);
    //         if (!padding) padding = 8;
    //
    //         buffer.resize(realLen + padding, char(0));
    //         buffer.shrink_to_fit();
    //         buffer.resize(realLen);
    //     }
    // }

    return buffer;
}

}  // End anonymous namespace


Foam::DynamicList<char>
Foam::IFstream::readContents_list(IFstream& ifs)
{
    return readContents_impl<DynamicList<char>>(ifs);
}


Foam::DynamicList<char>
Foam::IFstream::readContents_list(const fileName& pathname)
{
    if (!pathname.empty())
    {
        IFstream ifs(pathname, IOstreamOption::BINARY);

        if (ifs.good())
        {
            return readContents(ifs);
        }
    }

    return DynamicList<char>();
}


std::string
Foam::IFstream::readContents_string(IFstream& ifs)
{
    return readContents_impl<std::string>(ifs);
}


std::string
Foam::IFstream::readContents_string(const fileName& pathname)
{
    if (!pathname.empty())
    {
        IFstream ifs(pathname, IOstreamOption::BINARY);

        if (ifs.good())
        {
            return readContents_string(ifs);
        }
    }

    return std::string();
}


std::vector<char>
Foam::IFstream::readContents_vector(IFstream& ifs)
{
    return readContents_impl<std::vector<char>>(ifs);
}


std::vector<char>
Foam::IFstream::readContents_vector(const fileName& pathname)
{
    if (!pathname.empty())
    {
        IFstream ifs(pathname, IOstreamOption::BINARY);

        if (ifs.good())
        {
            return readContents_vector(ifs);
        }
    }

    return std::vector<char>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IFstream::IFstream
(
    const fileName& pathname,
    IOstreamOption streamOpt
)
:
    Foam::ifstreamPointer(pathname, streamOpt),
    ISstream(*(ifstreamPointer::get()), pathname, streamOpt)
{
    IOstreamOption::compression(ifstreamPointer::whichCompression());

    setClosed();

    setState(ifstreamPointer::get()->rdstate());

    if (good())
    {
        setOpened();
    }
    else
    {
        setBad();
    }

    lineNumber_ = 1;

    if (debug)
    {
        if (pathname.empty())
        {
            InfoInFunction
                << "Cannot open empty file name"
                << Foam::endl;
        }
        else if (IOstreamOption::COMPRESSED == IOstreamOption::compression())
        {
            InfoInFunction
                << "Decompressing " << (this->name() + ".gz") << Foam::endl;
        }

        if (!opened())
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for input\n" << info() << Foam::endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::streamsize Foam::IFstream::fileSize() const
{
    const std::istream* ptr = ifstreamPointer::get();

    if (!ptr || this->name().empty())
    {
        return std::streamsize(-1);
    }

    off_t fileLen = -1;

    if (IOstreamOption::COMPRESSED == ifstreamPointer::whichCompression())
    {
        fileLen = Foam::fileSize(this->name() + ".gz");
    }
    else
    {
        // TBD: special handing for wrapped icharstream
        // if
        // (
        //     const Foam::icharstream* charstr
        //   = dynamic_cast<const Foam::icharstream*>(ptr)>(ptr)
        // )
        // {
        //     return charstr->capacity();
        // }

        fileLen = Foam::fileSize(this->name());
    }

    if (fileLen >= 0)
    {
        return std::streamsize(fileLen);
    }

    return std::streamsize(-1);
}


std::istream& Foam::IFstream::stdStream()
{
    std::istream* ptr = ifstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


const std::istream& Foam::IFstream::stdStream() const
{
    const std::istream* ptr = ifstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


void Foam::IFstream::rewind()
{
    Istream::rewind();  // Drop any putback
    lineNumber_ = 1;    // Reset line number

    if (IOstreamOption::COMPRESSED == ifstreamPointer::whichCompression())
    {
        ifstreamPointer::reopen_gz(this->name());
        setState(ifstreamPointer::get()->rdstate());
    }
    else
    {
        ISstream::rewind();
    }
}


void Foam::IFstream::print(Ostream& os) const
{
    os  << "IFstream: ";
    ISstream::print(os);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::IFstream& Foam::IFstream::operator()() const
{
    if (!good())
    {
        // Also checks .gz file
        if (Foam::isFile(this->name(), true))
        {
            check(FUNCTION_NAME);
            FatalIOError.exit();
        }
        else
        {
            FatalIOErrorInFunction(*this)
                << "File " << this->name() << " does not exist"
                << exit(FatalIOError);
        }
    }

    return const_cast<IFstream&>(*this);
}


// ************************************************************************* //
