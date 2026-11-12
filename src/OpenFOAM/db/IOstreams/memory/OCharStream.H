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
    An output stream with internal array storage.
    Similar to OStringStream but with an array storage
    to enable reuse of its contents without copying.

See Also
    Foam::ICharStream
    Foam::ISpanStream
    Foam::OSpanStream

\*---------------------------------------------------------------------------*/

#ifndef Foam_OCharStream_H
#define Foam_OCharStream_H

#include "OSpanStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
                     Class Detail::ocharstream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::ostringstream, but with the ability to swap
//- character content.
//- Has some similarity to std::ospanstream (C++23)
template<class BufferType>
class ocharstream
:
    virtual public std::ios,
    protected BufferType,
    public std::ostream
{
    typedef BufferType buffer_type;
    typedef std::ostream stream_type;

public:

    // Data Types

        //- The dynamic container type of the buffer character content
        typedef typename buffer_type::dynamic_type dynamic_type;

        //- True if the buffer uses std::vector<char> storage
        static constexpr bool is_vector_buffer =
        std::is_base_of_v<std::vector<char>, dynamic_type>;


    // Constructors

        //- Default construct - empty
        ocharstream()
        :
            buffer_type(),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Move construct from std::vector (for vector buffer only)
        template
        <
            bool IsCharVector = is_vector_buffer,
            std::enable_if_t<IsCharVector, bool> = true
        >
        ocharstream(std::vector<char>&& buffer) : ocharstream()
        {
            buffer_type::swap(buffer);
            stream_type::clear();  // Clear old errors
        }

        //- Move construct from List (for List buffer only)
        template
        <
            bool IsCharList = !is_vector_buffer,
            std::enable_if_t<IsCharList, bool> = true
        >
        ocharstream(List<char>&& buffer) : ocharstream()
        {
            buffer_type::swap(buffer);
            stream_type::clear();  // Clear old errors
        }

        //- Move construct from DynamicList (for List buffer only)
        template
        <
            int SizeMin,
            bool IsCharList = !is_vector_buffer,
            std::enable_if_t<IsCharList, bool> = true
        >
        ocharstream(DynamicList<char,SizeMin>&& buffer) : ocharstream()
        {
            buffer_type::swap(buffer);
            stream_type::clear();  // Clear old errors
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

        //- Exchange stream content and parameter contents, reset positions
        //- (for vector buffer only)
        template
        <
            bool IsCharVector = is_vector_buffer,
            std::enable_if_t<IsCharVector, bool> = true
        >
        void swap(std::vector<char>& other)
        {
            buffer_type::swap(other);
            stream_type::clear();  // Clear old errors
        }

        //- Exchange stream content and parameter contents, reset positions
        //- (for List buffer only)
        template
        <
            bool IsCharList = !is_vector_buffer,
            std::enable_if_t<IsCharList, bool> = true
        >
        void swap(::Foam::List<char>& other)
        {
            buffer_type::swap(other);
            stream_type::clear();  // Clear old errors
        }

        //- Exchange stream content and parameter contents, reset positions
        //- (for List buffer only)
        template
        <
            int SizeMin,
            bool IsCharList = !is_vector_buffer,
            std::enable_if_t<IsCharList, bool> = true
        >
        void swap(::Foam::DynamicList<char,SizeMin>& other)
        {
            buffer_type::swap(other);
            stream_type::clear();  // Clear old errors
        }

        //- Reset stream and return contents as \c dynamic_type
        auto release()
        {
            auto chars(buffer_type::release());
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

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                         Class OCharStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream with internal array storage.
class OCharStream
:
    public Foam::Detail::StreamAllocator
    <Foam::Detail::ocharstream<Foam::memorybuf::out_dynamic>>,
    public Foam::OSstream
{
    //- The wrapped stream type (eg, ocharstream)
    using wrapped_type = stdFoam::remove_cvref_t<decltype(stream_)>;

public:

    // Data Types

        //- The dynamic container type of the buffer character content
        typedef typename wrapped_type::dynamic_type dynamic_content_type;

        //- True if the buffer uses std::vector<char> storage
        static constexpr bool is_vector_buffer = wrapped_type::is_vector_buffer;


    // Constructors

        //- Default construct (empty output)
        explicit OCharStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OSstream(stream_, "output", streamOpt.format(), streamOpt.version())
        {}

        //- Move construct from std::vector of initial storage
        //- (for vector buffer only)
        template
        <
            class StreamType = wrapped_type,
            std::enable_if_t<StreamType::is_vector_buffer, bool> = true
        >
        explicit OCharStream
        (
            std::vector<char>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OCharStream(streamOpt)
        {
            static_cast<StreamType&>(stream_).swap(buffer);
        }

        //- Move construct from a List of initial storage
        //- (for List buffer only)
        template
        <
            class StreamType = wrapped_type,
            std::enable_if_t<!StreamType::is_vector_buffer, bool> = true
        >
        explicit OCharStream
        (
            ::Foam::List<char>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OCharStream(streamOpt)
        {
            static_cast<StreamType&>(stream_).swap(buffer);
        }

        //- Move construct from DynamicList of initial storage
        //- (uses entire capacity).
        //- (for List buffer only).
        template
        <
            int SizeMin,
            class StreamType = wrapped_type,
            std::enable_if_t<!StreamType::is_vector_buffer, bool> = true
        >
        explicit OCharStream
        (
            ::Foam::DynamicList<char,SizeMin>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            OCharStream(streamOpt)
        {
            static_cast<StreamType&>(stream_).swap(buffer);
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

        //- For OStringStream compatibility, return the buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        auto str() const { return stream_.str(); }

        //- A modifiable span (view) of the buffer contents.
        //- The caller really knows what they are doing.
        stdFoam::span<char> span() { return stream_.span(); }


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


    // Additional content management

        //- Exchange stream content and parameter contents, reset positions
        template
        <
            class StreamType = wrapped_type,
            std::enable_if_t<StreamType::is_vector_buffer, bool> = true
        >
        void swap(std::vector<char>& other)
        {
            static_cast<StreamType&>(stream_).swap(other);
            syncState();
        }

        //- Exchange stream content and parameter contents, reset positions
        template
        <
            class StreamType = wrapped_type,
            std::enable_if_t<!StreamType::is_vector_buffer, bool> = true
        >
        void swap(Foam::List<char>& other)
        {
            static_cast<StreamType&>(stream_).swap(other);
            syncState();
        }

        //- Exchange stream content and parameter contents, reset positions
        //- (for List buffer only)
        template
        <
            int SizeMin,
            class StreamType = wrapped_type,
            std::enable_if_t<!StreamType::is_vector_buffer, bool> = true
        >
        void swap(Foam::DynamicList<char,SizeMin>& other)
        {
            static_cast<StreamType&>(stream_).swap(other);
            syncState();
        }

        //- Reset stream and return character content
        //- as \c dynamic_content_type
        auto release()
        {
            auto chars(stream_.release());
            syncState();
            return chars;
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Typedefs

//- Simple ocharstream with List storage
typedef Detail::ocharstream<Foam::memorybuf::out_dynamic> ocharstream;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
