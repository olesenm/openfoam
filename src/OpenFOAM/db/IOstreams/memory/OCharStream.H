/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2025 OpenCFD Ltd.
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
    Foam::OCharStream

Description
    An output stream that writes to a List and manages the List storage.
    Similar to OStringStream but with a List for its storage instead of
    as string to allow reuse of List contents without copying.

    Internally imposes a 512 byte min-size and uses capacity doubling.

See Also
    Foam::ICharStream
    Foam::OSpanStream
    Foam::ISpanStream

\*---------------------------------------------------------------------------*/

#ifndef Foam_OCharStream_H
#define Foam_OCharStream_H

#include "OSpanStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ocharstream;
class OCharStream;

// Older names (prior to 2023-08)
typedef OCharStream OListStream;


/*---------------------------------------------------------------------------*\
                        Class ocharstream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::ostringstream, but with the ability to swap
//- character content.
//- Has some similarity to std::ospanstream (C++23)
class ocharstream
:
    virtual public std::ios,
    protected Foam::memorybuf::out_dynamic,
    public std::ostream
{
    typedef Foam::memorybuf::out_dynamic buffer_type;
    typedef std::ostream stream_type;

public:

    // Constructors

        //- Default construct - empty
        ocharstream()
        :
            buffer_type(),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Move construct from List
        ocharstream(List<char>&& buffer)
        :
            ocharstream()
        {
            swap(buffer);
        }

        //- Move construct from DynamicList
        template<int SizeMin>
        ocharstream(DynamicList<char,SizeMin>&& buffer)
        :
            ocharstream()
        {
            swap(buffer);
        }


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

        //- Reserve output space for at least this amount.
        //- Applies a min-size and capacity doubling.
        void reserve(std::streamsize n)
        {
            buffer_type::reserve(n);
        }

        //- Reserve output space for at least this amount.
        //- Does not apply min-size or capacity doubling etc.
        void reserve_exact(std::streamsize n)
        {
            buffer_type::reserve_exact(n);
        }

        //- Increase (reserve) space for another \c count entries
        void extend(std::streamsize count)
        {
            buffer_type::extend(count);
        }

        //- Increase (reserve) space for another \c count entries
        void extend_exact(std::streamsize count)
        {
            buffer_type::extend_exact(count);
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

        //- A list \em span of current output contents (is modifiable!!)
        UList<char> list() const
        {
            return UList<char>
            (
                buffer_type::data_bytes(),
                label(buffer_type::size_bytes())
            );
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

        //- Exchange stream content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            buffer_type::swap(other);
            stream_type::clear();  // Clear old errors
        }

        //- Exchange stream content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            buffer_type::swap(other);
            stream_type::clear();  // Clear old errors
        }

        //- Reset buffer and return contents
        DynamicList<char> release()
        {
            DynamicList<char> chars(buffer_type::release());
            stream_type::clear();  // Clear old errors
            return chars;
        }

        //- Some information about the output buffer position/capacity
        void debug_info(Ostream& os) const
        {
            os << "put=" << output_pos() << '/' << capacity();
        }

        //- Information about stream
        void print(Ostream& os) const { debug_info(os); os << '\n'; }


    // Extra/Convenience Methods

        //- Append a single character to the end
        void push_back(char c) { stream_type::put(c); }

        //- Rewind the end by 1 or more elements
        void pop_back(int n = 1) { buffer_type::pop_back(n); }

        //- Append repeated character content
        void append(std::streamsize count, char c)
        {
            if (count > 0)
            {
                buffer_type::extend(count);
                while (count-- > 0)
                {
                    stream_type::put(c);
                }
            }
        }

        //- Append character content - like a plain write()
        void append(const char* data, std::streamsize count)
        {
            if (data && count > 0)
            {
                buffer_type::extend(count);
                write(data, count);
            }
        }

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
                         Class OCharStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream with internal List storage
class OCharStream
:
    public Foam::Detail::StreamAllocator<Foam::ocharstream>,
    public Foam::OSstream
{
    typedef
        Foam::Detail::StreamAllocator<Foam::ocharstream>
        allocator_type;

public:

    // Constructors

        //- Default construct (empty output)
        explicit OCharStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "output", streamOpt.format(), streamOpt.version())
        {}

        //- Move construct from a List of initial storage
        explicit OCharStream
        (
            ::Foam::List<char>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OCharStream(streamOpt)
        {
            stream_.swap(buffer);
        }

        //- Move construct from a DynamicList of initial storage
        //- (uses entire capacity)
        template<int SizeMin>
        explicit OCharStream
        (
            ::Foam::DynamicList<char,SizeMin>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OCharStream(streamOpt)
        {
            stream_.swap(buffer);
        }


    // Member Functions

        //- Position of the put buffer
        std::streampos tellp() const { return stream_.output_pos(); }

        //- The current output position within the buffer (tellp)
        std::streampos output_pos() const { return stream_.output_pos(); }

        //- The number of bytes outputted
        std::streamsize count() const { return stream_.count(); }

        //- The current output size. Same as count(), output_pos(), tellp().
        label size() const { return label(stream_.count()); }

        //- The put buffer capacity
        std::streamsize capacity() const { return stream_.capacity(); }

        //- Reserve output space for at least this amount
        void reserve(std::streamsize n) { stream_.reserve(n); }

        //- Reserve output space for at least this amount.
        //- Does not apply min-size or capacity doubling etc.
        void reserve_exact(std::streamsize n) { stream_.reserve_exact(n); }

        //- Increase (reserve) space for another \c n entries
        void extend(std::streamsize n) { stream_.extend(n); }

        //- Increase (reserve) space for another \c n entries
        void extend_exact(std::streamsize n) { stream_.extend_exact(n); }

        //- A string_view of buffer contents
        std::string_view view() const { return stream_.view(); }

        //- A sub-slice string view of the buffer contents
        std::string_view view(size_t pos, size_t len = std::string::npos) const
        {
            return stream_.view(pos, len);
        }

        //- A list \em span of the current output characters (is modifiable!)
        UList<char> list() const { return stream_.list(); }

        //- For OStringStream compatibility, return the buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        auto str() const { return stream_.str(); }

        //- Exchange stream content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            stream_.swap(other);
            syncState();
        }

        //- Exchange stream content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            stream_.swap(other);
            syncState();
        }

        //- Reset buffer and return contents
        DynamicList<char> release()
        {
            DynamicList<char> chars(stream_.release());
            syncState();
            return chars;
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
            os << "ocharstream: ";
            stream_.debug_info(os);
            os << '\n';
        }


    // Extra/Convenience Methods

        //- Append a single character to the end
        void push_back(char c) { stream_.push_back(c); }

        //- Rewind the end by 1 or more elements
        void pop_back(int n = 1) { stream_.pop_back(n); }

        //- Append repeated character content
        void append(std::streamsize count, char c)
        {
            stream_.append(count, c);
        }

        //- Append character content
        void append(const char* data, std::streamsize count)
        {
            stream_.append(data, count);
        }

        //- Append character content
        void append(std::string_view sv)
        {
            stream_.append(sv.data(), sv.size());
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


    // Housekeeping

        //- Block size was used in OpenFOAM-v2306 and earlier
        void setBlockSize(int n) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
