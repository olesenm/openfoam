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
    Foam::ICharStream

Description
    An input stream that reads from a List and manages the List storage.
    Similar to IStringStream but with a List for its storage instead of
    as string to allow reuse of List contents without copying.

See Also
    Foam::ISpanStream
    Foam::OCharStream
    Foam::OSpanStream

\*---------------------------------------------------------------------------*/

#ifndef Foam_ICharStream_H
#define Foam_ICharStream_H

#include "ISpanStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class icharstream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::istringstream, but with the ability to swap
//- character content.
//- Has some similarity to std::ispanstream (C++23)
class icharstream
:
    virtual public std::ios,
    protected Foam::memorybuf::in_dynamic,
    public std::istream
{
    typedef Foam::memorybuf::in_dynamic buffer_type;
    typedef std::istream stream_type;

public:

    // Data Types

        //- The dynamic container type of the buffer character content
        typedef typename buffer_type::dynamic_type dynamic_type;

        //- True if the buffer uses std::vector<char> storage
        static constexpr bool is_vector_buffer =
        std::is_base_of_v<std::vector<char>, dynamic_type>;


    // Constructors

        //- Default construct - empty
        icharstream()
        :
            buffer_type(),
            stream_type(static_cast<buffer_type*>(this))
        {}

        //- Copy construct from content
        icharstream(const char* buffer, size_t nbytes) : icharstream()
        {
            reset(buffer, nbytes);
        }

        //- Move construct from List
        icharstream(List<char>&& buffer) : icharstream()
        {
            swap(buffer);
        }

        //- Move construct from DynamicList
        template<int SizeMin>
        icharstream(DynamicList<char,SizeMin>&& buffer) : icharstream()
        {
            swap(buffer);
        }


    // Member Functions

        //- The current get position within the buffer (tellg)
        std::streampos input_pos() const
        {
            return buffer_type::span_tellg();
        }

        //- The get buffer capacity
        std::streamsize capacity() const
        {
            return buffer_type::span_capacity();
        }

        //- The number of characters remaining in the get area.
        //- Same as (capacity() - input_pos())
        std::streamsize remaining() const
        {
            return buffer_type::span_remaining();
        }

        //- A string_view of buffer contents
        std::string_view view() const { return buffer_type::view(); }

        //- A sub-slice string view of the buffer contents
        std::string_view view(size_t pos, size_t len = std::string::npos) const
        {
            return buffer_type::view(pos, len);
        }

        //- For istringstream compatibility, return the buffer as string copy.
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

        //- Rewind the stream, clearing any old errors
        void rewind()
        {
            buffer_type::pubseekpos(0, std::ios_base::in);
            stream_type::clear();  // Clear old errors
        }

        //- Reposition the stream from the start
        void seek(std::streampos pos)
        {
            if (buffer_type::in_range(pos))
            {
                buffer_type::pubseekpos(pos, std::ios_base::in);
                stream_type::clear();  // Clear old errors
            }
        }

        //- Reset stream content (copy), reset positions
        void reset(const char* buffer, size_t nbytes)
        {
            buffer_type::reset(buffer, nbytes);
            stream_type::clear();  // Clear old errors
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

        //- Reset buffer and return contents as \c dynamic_type
        auto release()
        {
            auto chars(buffer_type::release());
            stream_type::clear();  // Clear old errors
            return chars;
        }

        //- Some information about the input buffer position/capacity
        void debug_info(Ostream& os) const
        {
            os << "get=" << input_pos() << '/' << capacity();
        }

        //- Information about stream
        void print(Ostream& os) const { debug_info(os); os << '\n'; }
};


/*---------------------------------------------------------------------------*\
                         Class ICharStream Declaration
\*---------------------------------------------------------------------------*/

//- An ISstream with internal List storage. Always UNCOMPRESSED.
class ICharStream
:
    public Foam::Detail::StreamAllocator<Foam::icharstream>,
    public Foam::ISstream
{
    //- The wrapped stream type (eg, icharstream)
    using wrapped_type = stdFoam::remove_cvref_t<decltype(stream_)>;

public:

    // Data Types

        //- The dynamic container type of the buffer character content
        typedef typename wrapped_type::dynamic_type dynamic_content_type;

        //- True if the buffer uses std::vector<char> storage
        static constexpr bool is_vector_buffer = wrapped_type::is_vector_buffer;


    // Constructors

        //- Default construct (empty), optionally with specified stream option
        explicit ICharStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Copy construct from string content
        explicit ICharStream
        (
            const std::string& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ICharStream(streamOpt)
        {
            stream_.reset(buffer.data(), buffer.size());
        }

        //- Move construct from List
        explicit ICharStream
        (
            List<char>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ICharStream(streamOpt)
        {
            stream_.swap(buffer);
        }

        //- Move construct from DynamicList (uses current size)
        template<int SizeMin>
        explicit ICharStream
        (
            DynamicList<char,SizeMin>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ICharStream(streamOpt)
        {
            stream_.swap(buffer);
        }


    // Member Functions

    // Access/Query

        //- Position of the get buffer
        std::streampos tellg() const { return stream_.input_pos(); }

        //- The current get position within the buffer (tellg)
        std::streampos input_pos() const { return stream_.input_pos(); }

        //- The input list size. Same as capacity()
        std::streamsize size() const { return stream_.capacity(); }

        //- The get buffer capacity
        std::streamsize capacity() const { return stream_.capacity(); }

        //- The number of characters remaining in the get area.
        //- Same as (capacity() - input_pos())
        std::streamsize remaining() const { return stream_.remaining(); }

        //- A string_view of buffer contents
        std::string_view view() const { return stream_.view(); }

        //- A sub-slice string view of the buffer contents
        std::string_view view(size_t pos, size_t len = std::string::npos) const
        {
            return stream_.view(pos, len);
        }

        //- For IStringStream compatibility, return the buffer as string copy.
        //  Use sparingly - it creates a full copy!!
        auto str() const { return stream_.str(); }

        //- A modifiable span (view) of the buffer contents.
        //- The caller really knows what they are doing.
        stdFoam::span<char> span() { return stream_.span(); }


    // Edit

        //- Reset content (copy)
        void reset(const char* buffer, size_t nbytes)
        {
            stream_.reset(buffer, nbytes);
            syncState();
        }

        //- Rewind the stream, clearing any old errors
        virtual void rewind() override
        {
            Istream::rewind();
            stream_.rewind();
            syncState();
        }

        //- Reposition the stream from the start
        void seek(std::streampos pos)
        {
            stream_.seek(pos);
            syncState();
        }


    // Additional content management

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

        //- Reset stream and return character content
        //- as \c dynamic_content_type
        auto release()
        {
            auto chars(stream_.release());
            syncState();
            return chars;
        }


    // Other

        //- Print stream description
        virtual void print(Ostream& os) const override
        {
            os << "icharstream: ";
            stream_.debug_info(os);
            os << '\n';
        }


    // Member Operators

        //- A non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary
        Istream& operator()() const
        {
            // Could also rewind
            return const_cast<ICharStream&>(*this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
