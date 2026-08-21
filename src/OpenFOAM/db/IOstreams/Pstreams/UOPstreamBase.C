/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2024 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "UOPstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Return the position with specified byte alignment (eg, word boundary)
template<typename IntType>
inline IntType byteAlign(IntType pos, size_t align) noexcept
{
    return
    (
        (align > 1)
      ? (align + ((pos - 1) & ~(align - 1)))
      : pos
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UOPstreamBase::alignp(size_t alignment)
{
    // Aligned position
    const auto pos = byteAlign(sendBuf_.size(), alignment);

    // Move to the aligned output position. Fill any gap with nul char.
    sendBuf_.resize(pos, '\0');
}


inline void Foam::UOPstreamBase::prepareBuffer
(
    const size_t count,
    const size_t align
)
{
    if (!count)
    {
        return;
    }

    // Align for the next output position
    const auto pos = byteAlign(sendBuf_.size(), align);

    const int64_t newLen = int64_t(pos) + int64_t(count);

    // The max buffer length:
    // - limited by OpenFOAM label size, or by MPI <int> limit (prior to MPI-4)

    constexpr label minLen = 1024;

    constexpr int64_t maxLen =
    (
        (sizeof(UList<char>::max_size()) < sizeof(int64_t))
      ? int64_t(UList<char>::max_size())            // OpenFOAM (label)
      : int64_t(std::numeric_limits<int>::max())    // MPI (int)
    );

    if (FOAM_UNLIKELY(newLen > maxLen))
    {
        FatalErrorInFunction
            << "Buffer requirement exceeds <int32> range by "
            << (newLen - maxLen) << " values.\n";

        if constexpr (sizeof(UList<char>::max_size()) < sizeof(int64_t))
        {
            // Requested amount exceeds OpenFOAM (int32) label limit
            FatalError
                << "Using an <int64> build will not likely solve the problem"
                   " since MPI itself still has <int> limits\n";
        }
        else
        {
            // Requested amount exceeds MPI (int) limit.
            // Need to wait for MPI-4 (big sizes)
            FatalError
                << "This is still an MPI limit (as of AUG-2026)\n";
        }
        FatalError<< Foam::abort(FatalError);
    }


    // The int32 versions need to be careful when resizing a buffer of bytes:
    //
    // Currently DynamicList::reserve() uses doubling without additional
    // guards. Use a slower (1.5) growth when the buffer becomes larger.

    if constexpr (sizeof(UList<char>::max_size()) < sizeof(int64_t))
    {
        // OpenFOAM <label> limit:
        constexpr int64_t threshold = int64_t(0.4*maxLen);

        const int64_t oldCapacity(sendBuf_.capacity());
        if ((newLen > oldCapacity) && (oldCapacity > threshold))
        {
            // Slower (3/2 = 1.5x) growth to avoid doubling overflow,
            // clamped to UList::max_size().
            // No min-size clamping needed for this code path.

            int64_t newCapacity =
            (
                Foam::ListPolicy::reserve_size<0,3,2>(newLen, oldCapacity)
            );
            newCapacity = std::min(newCapacity, maxLen);

            // Impose this capacity directly
            sendBuf_.reserve_exact(label(newCapacity));
        }
        else
        {
            // Extend buffer (as required)
            sendBuf_.reserve(Foam::max(minLen, label(newLen)));
        }
    }
    else
    {
        // MPI <int> limit:

        // Extend buffer (as required)
        sendBuf_.reserve(Foam::max(minLen, label(newLen)));
    }

    // Move to the aligned output position. Fill any gap with nul char.
    sendBuf_.resize(pos, '\0');
}


template<class T>
inline void Foam::UOPstreamBase::writeToBuffer(const T& val)
{
    writeToBuffer(&val, sizeof(T), sizeof(T));
}


inline void Foam::UOPstreamBase::writeToBuffer
(
    const void* data,
    const size_t count,
    const size_t align
)
{
    if (!count)
    {
        // Nothing to do, don't even do any alignment
        return;
    }

    prepareBuffer(count, align);

    // The aligned output position
    const auto pos = sendBuf_.size();

    if (!data)
    {
        // No data to write. Extend and fill with block of zeros
        sendBuf_.resize(pos + count, '\0');
        return;
    }

    // Extend the addressable range for direct pointer access
    sendBuf_.resize(pos + count);

    char* const __restrict__ buf = (sendBuf_.data() + pos);
    const char* const __restrict__ input = reinterpret_cast<const char*>(data);

    std::copy(input, input + count, buf);
}


void Foam::UOPstreamBase::writeToBuffer_at
(
    // The target output position
    int64_t pos,
    const void* data,
    const size_t count,
    const size_t align
)
{
    if (!count)
    {
        // Nothing to do
        return;
    }

    if (pos < 0)
    {
        prepareBuffer(count, align);
        pos = sendBuf_.size();
    }
    else
    {
        // Align the position to the specified boundary
        pos = byteAlign(pos, align);
    }

    // The end position after the write
    const int64_t end_pos = pos + static_cast<int64_t>(count);

    // Need additional size check, as per prepareBuffer? (issue #3152)

    // Extend buffer if needed to accommodate the write at pos
    if (sendBuf_.size() < end_pos)
    {
        sendBuf_.resize(end_pos, '\0');
    }

    if (!data)
    {
        // No data to write. Fill with zeros at the target position
        std::fill(sendBuf_.data() + pos, sendBuf_.data() + end_pos, '\0');
        return;
    }

    // Write data directly at the target position
    char* const __restrict__ buf = (sendBuf_.data() + pos);
    const char* const __restrict__ input = reinterpret_cast<const char*>(data);

    std::copy(input, input + count, buf);
}


inline void Foam::UOPstreamBase::putChar(const char c)
{
    if (!sendBuf_.capacity())
    {
        sendBuf_.setCapacity(1024);
    }
    sendBuf_.push_back(c);
}


inline void Foam::UOPstreamBase::putString
(
    const char* str,
    const size_t len
)
{
    writeToBuffer(len);
    writeToBuffer(str, len, 1);  // no-op when len == 0
}


inline void Foam::UOPstreamBase::putString(const std::string& str)
{
    putString(str.data(), str.size());
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UOPstreamBase::UOPstreamBase
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    DynamicList<char>& sendBuf,
    const int tag,
    const int communicator,
    const bool sendAtDestruct,
    IOstreamOption::streamFormat fmt
)
:
    UPstream(commsType),
    Ostream(fmt),
    toProcNo_(toProcNo),
    tag_(tag),
    comm_(communicator),
    sendAtDestruct_(sendAtDestruct),
    sendBuf_(sendBuf)
{
    setOpened();
    setGood();
}


Foam::UOPstreamBase::UOPstreamBase(const int toProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType()),
    Ostream(buffers.format()),
    toProcNo_(toProcNo),
    tag_(buffers.tag()),
    comm_(buffers.comm()),
    sendAtDestruct_(buffers.commsType() != UPstream::commsTypes::nonBlocking),
    sendBuf_(buffers.accessSendBuffer(toProcNo))
{
    setOpened();
    setGood();
}


Foam::UOPstreamBase::UOPstreamBase
(
    DynamicList<char>& sendBuf,
    IOstreamOption::streamFormat fmt
)
:
    UPstream(UPstream::commsTypes::nonBlocking), // placeholder
    Ostream(fmt),
    toProcNo_(UPstream::masterNo()),        // placeholder
    tag_(UPstream::msgType()),              // placeholder
    comm_(UPstream::commSelf()),            // placeholder
    sendAtDestruct_(false),   // Never sendAtDestruct!!
    sendBuf_(sendBuf)
{
    sendBuf_.clear();  // Overwrite into buffer
    setOpened();
    setGood();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UOPstreamBase::~UOPstreamBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::UOPstreamBase::write(const token& tok)
{
    // Direct token handling only for some types

    switch (tok.type())
    {
        case token::tokenType::FLAG :
        {
            putChar(token::tokenType::FLAG);
            putChar(tok.flagToken());

            return true;
        }

        // The word-variants
        case token::tokenType::WORD :
        case token::tokenType::DIRECTIVE :
        {
            putChar(tok.type());
            putString(tok.wordToken());

            return true;
        }

        // The string-variants
        case token::tokenType::STRING :
        case token::tokenType::EXPRESSION :
        case token::tokenType::VARIABLE :
        case token::tokenType::VERBATIM :
        case token::tokenType::CHAR_DATA :
        {
            putChar(tok.type());
            putString(tok.stringToken());

            return true;
        }

        default:
            break;
    }

    return false;
}


Foam::Ostream& Foam::UOPstreamBase::write(const char c)
{
    if (!isspace(c))
    {
        putChar(c);
    }

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::writeQuoted
(
    const char* str,
    std::streamsize len,
    const bool quoted
)
{
    if (quoted)
    {
        putChar(token::tokenType::STRING);
    }
    else
    {
        putChar(token::tokenType::WORD);
    }
    putString(str, len);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const char* str)
{
    const word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        return write(nonWhiteChars[0]);
    }
    else if (nonWhiteChars.size())
    {
        return write(nonWhiteChars);
    }

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const word& str)
{
    putChar(token::tokenType::WORD);
    putString(str);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const std::string& str)
{
    putChar(token::tokenType::STRING);
    putString(str);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const int32_t val)
{
    putChar(token::tokenType::INTEGER_32);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const int64_t val)
{
    putChar(token::tokenType::INTEGER_64);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const uint32_t val)
{
    putChar(token::tokenType::UNSIGNED_INTEGER_32);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const uint64_t val)
{
    putChar(token::tokenType::UNSIGNED_INTEGER_64);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const float val)
{
    putChar(token::tokenType::FLOAT);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write(const double val)
{
    putChar(token::tokenType::DOUBLE);
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::write
(
    const char* data,
    std::streamsize count
)
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    // For count=0 (no-op, also no alignment) and data=nullptr are handled
    // properly in writeToBuffer

    // Align on word boundary (64-bit)
    writeToBuffer(data, count, 8);

    return *this;
}


Foam::Ostream& Foam::UOPstreamBase::writeRaw
(
    const char* data,
    std::streamsize count
)
{
    // No check for IOstreamOption::BINARY since this is either done in the
    // beginRawWrite() method, or the caller knows what they are doing.

    // For count=0 (no-op, also no alignment) and data=nullptr are handled
    // properly in writeToBuffer

    // Previously aligned and sizes reserved via beginRawWrite() etc.
    writeToBuffer(data, count, 1);

    return *this;
}


bool Foam::UOPstreamBase::beginRawWrite(std::streamsize count)
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    // Align on word boundary (64-bit)
    // - as per write(const char*, streamsize)
    // - for count=0 this is a no-op
    prepareBuffer(count, 8);

    return true;
}


void Foam::UOPstreamBase::rewind()
{
    sendBuf_.clear();  // Overwrite into buffer
    setOpened();
    setGood();
}


void Foam::UOPstreamBase::print(Ostream& os) const
{
    os  << "Writing to processor " << toProcNo_
        << " from processor " << myProcNo() << " in communicator " << comm_
        << " and tag " << tag_ << Foam::endl;
}


// ************************************************************************* //
