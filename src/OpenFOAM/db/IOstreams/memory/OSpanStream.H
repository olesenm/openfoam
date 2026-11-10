/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

Class
    Foam::OSpanStream

Description
    Similar to OCharStream but using an externally managed buffer for
    its output.

    This allows the output buffer to be reused and can make it easier when
    writing out data.  It is the user's responsibility to ensure proper
    synchronization in the sizes. Provided that the external buffer is large
    enough that overflow does not occur, the following usage pattern
    works.

    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        {
            OSpanStream os(buffer);
            os << "content1" << " and more content";
            buffer.resize(os.size());       // synchronize sizes
        }

        something.write(buffer, buffer.size());
    \endcode

    Although the OSpanStream is quite lightweight, there may be cases
    where it is preferable to reuse the stream as well.
    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        OSpanStream os(buffer);
        os << "content1" << " and more content";
        buffer.resize(os.size());           // synchronize sizes

        something.write(buffer, buffer.size());

        os.rewind();
        os << "content2";
        buffer.resize(os.size());      // synchronize sizes

        something.write(buffer, buffer.size());

        // or simply using the output size directly (without sync)
        os.rewind();
        os << "content3";

        something.write(buffer, os.size());
    \endcode

See Also
    Foam::ICharStream
    Foam::ISpanStream
    Foam::OCharStream

\*---------------------------------------------------------------------------*/

#ifndef Foam_OSpanStream_H
#define Foam_OSpanStream_H

#include "memoryStreamBuffer.H"
#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ospanstream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::ostringstream, but with an externally managed output buffer
//- which makes it most similar to std::ospanstream (C++23)
class ospanstream
:
    virtual public std::ios,
    protected Foam::memorybuf::out_base,
    public std::ostream
{
    typedef Foam::memorybuf::out_base buffer_type;
    typedef std::ostream stream_type;

public:

    // Constructors

        //- Default construct - empty
        ospanstream()
        :
            buffer_type(),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Construct for character array and number of bytes
        ospanstream(char* buffer, size_t nbytes)
        :
            buffer_type(buffer, nbytes),
            stream_type(static_cast<buffer_type*>(this))
        {}


    // Member Functions

        //- The current output position within the buffer (tellp)
        std::streampos output_pos() const
        {
            return buffer_type::span_tellp();
        }

        //- The number of bytes outputted
        std::streamsize count() const
        {
            return buffer_type::size_bytes();
        }

        //- The put buffer capacity
        std::streamsize capacity() const
        {
            return buffer_type::span_capacity();
        }

        //- Rewind the stream, clearing any old errors
        void rewind()
        {
            buffer_type::pubseekpos(0, std::ios_base::out);
            stream_type::clear();  // Clear old errors
        }

        //- Reposition the stream from the start
        void seek(std::streampos pos)
        {
            if (buffer_type::in_range(pos))
            {
                buffer_type::pubseekpos(pos, std::ios_base::out);
                stream_type::clear();  // Clear old errors
            }
        }

        //- A string_view of buffer contents
        std::string_view view() const { return buffer_type::view(); }

        //- A sub-slice string view of the buffer contents
        std::string_view view(size_t pos, size_t len = std::string::npos) const
        {
            return buffer_type::view(pos, len);
        }

        //- For ostringstream compatibility, return the buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        std::string str() const
        {
            return std::string
            (
                buffer_type::data_bytes(),
                buffer_type::size_bytes()
            );
        }

        //- A modifiable span (view) of the buffer contents.
        //- The caller really knows what they are doing.
        stdFoam::span<char> span()
        {
            return stdFoam::span<char>
            (
                buffer_type::data_bytes(),
                buffer_type::size_bytes()
            );
        }


        //- Reset the put buffer area
        void reset(char* buffer, size_t nbytes)
        {
            buffer_type::resetp(buffer, nbytes);
            stream_type::clear();  // Clear old errors
        }

        //- Reset the put buffer area to use the data area from a string
        void reset(std::string& s)
        {
            s.resize(s.capacity());
            buffer_type::resetp(s.data(), s.size());
            stream_type::clear();  // Clear old errors
        }

        //- Some information about the output buffer position/capacity
        void debug_info(Ostream& os) const
        {
            os << "put=" << output_pos() << '/' << capacity();
        }

        //- Information about stream
        void print(Ostream& os) const { debug_info(os); os << '\n'; }


    // Extra/Convenience Methods

        //- Overwrite a single character
        void overwrite(std::streampos pos, char c)
        {
            buffer_type::overwrite(pos, c);
        }

        //- Overwrite a sub-slice with character content
        void overwrite
        (
            std::streampos pos,
            const char* data,
            std::streamsize count
        )
        {
            buffer_type::overwrite(pos, data, count);
        }

        //- The output data (start of output characters)
        const char* cdata_bytes() const { return buffer_type::data_bytes(); }

        //- The output data (start of output characters)
        char* data_bytes() { return buffer_type::data_bytes(); }

        //- The current number of output characters
        std::streamsize size_bytes() const { return buffer_type::size_bytes(); }
};


/*---------------------------------------------------------------------------*\
                         Class OSpanStream Declaration
\*---------------------------------------------------------------------------*/

class OSpanStream
:
    public Foam::Detail::StreamAllocator<Foam::ospanstream>,
    public Foam::OSstream
{
    typedef
        Foam::Detail::StreamAllocator<Foam::ospanstream>
        allocator_type;

public:

    // Constructors

        //- Default construct (empty output)
        explicit OSpanStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "output", streamOpt.format(), streamOpt.version())
        {}

        //- Use data area from string content
        explicit OSpanStream
        (
            std::string& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(streamOpt)
        {
            stream_.reset(buffer);
        }

        //- Construct using specified buffer and number of bytes
        OSpanStream
        (
            char* buffer,
            size_t nbytes,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(streamOpt)
        {
            stream_.reset(buffer, nbytes);
        }

        //- Construct using data area from a List and its inherent storage size
        explicit OSpanStream
        (
            ::Foam::UList<char>& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(buffer.data(), buffer.size(), streamOpt)
        {}

        //- Construct using full data area from DynamicList
        template<int SizeMin>
        explicit OSpanStream
        (
            ::Foam::DynamicList<char,SizeMin>& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSpanStream(buffer.data(), buffer.capacity(), streamOpt)
        {
            buffer.resize(buffer.capacity());  // Uses entire space
        }


    // Member Functions

        //- Position of the put buffer
        std::streampos tellp() const { return stream_.output_pos(); }

        //- The current output position within the buffer (tellp)
        std::streampos output_pos() const { return stream_.output_pos(); }

        //- The number of bytes outputted
        std::streamsize count() const { return stream_.count(); }

        //- The current output size. Same as count(), output_pos(), tellp().
        std::streamsize size() const { return stream_.count(); }

        //- The put buffer capacity
        std::streamsize capacity() const { return stream_.capacity(); }

        //- A string_view of buffer contents
        std::string_view view() const { return stream_.view(); }

        //- A sub-slice string view of the buffer contents
        std::string_view view(size_t pos, size_t len = std::string::npos) const
        {
            return stream_.view(pos, len);
        }

        //- For OStringStream compatibility, return buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        auto str() const { return stream_.str(); }

        //- A modifiable span (view) of the buffer contents.
        //- The caller really knows what they are doing.
        stdFoam::span<char> span() { return stream_.span(); }


        //- Reset the put area
        void reset(char* buffer, size_t nbytes)
        {
            stream_.reset(buffer, nbytes);
            syncState();
        }

        //- Reset the put buffer area to use the data area from a string
        void reset(std::string& s)
        {
            stream_.reset(s);
            syncState();
        }

        //- Rewind the stream, clearing any old errors
        virtual void rewind()
        {
            stream_.rewind();
            syncState();
        }

        //- Reposition the stream from the start
        void seek(std::streampos pos)
        {
            stream_.seek(pos);
            syncState();
        }

        //- Print stream description
        virtual void print(Ostream& os) const override
        {
            os << "ospanstream: ";
            stream_.debug_info(os);
            os << '\n';
        }


    // Extra/Convenience Methods

        //- Append character content - same as writeRaw()
        void append(const char* data, std::streamsize count)
        {
            if (data && count > 0)
            {
                writeRaw(data, count);
            }
        }

        //- Append character content - same as writeRaw()
        void append(std::string_view sv)
        {
            append(sv.data(), sv.size());
        }

        //- Overwrite a single character
        void overwrite(std::streampos pos, char c)
        {
            stream_.overwrite(pos, c);
        }

        //- Overwrite a sub-slice with character content
        void overwrite
        (
            std::streampos pos,
            const char* data,
            std::streamsize count
        )
        {
            stream_.overwrite(pos, data, count);
        }

        //- Overwrite a sub-slice with character content
        void overwrite(std::streampos pos, std::string_view sv)
        {
            stream_.overwrite(pos, sv.data(), sv.size());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
