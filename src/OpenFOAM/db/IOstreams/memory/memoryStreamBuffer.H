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
    Foam::memorybuf

Description
    A std::streambuf used for memory buffer streams such as
    ispanstream, ocharstream, etc.

\*---------------------------------------------------------------------------*/

#ifndef Foam_memoryStreamBuffer_H
#define Foam_memoryStreamBuffer_H

#include "DynamicList.H"

#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class memorybuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf for memory similar to std::spanbuf (C++23)
class memorybuf
:
    public std::streambuf
{
protected:

    //- Set position pointer to relative position
    virtual std::streampos seekoff
    (
        std::streamoff off,
        std::ios_base::seekdir way,
        std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
    )
    {
        const bool testin  = which & std::ios_base::in;
        const bool testout = which & std::ios_base::out;

        if (way == std::ios_base::beg)
        {
            if (testin)
            {
                setg(eback(), eback(), egptr());
                gbump(off);
            }
            if (testout)
            {
                setp(pbase(), epptr());
                pbump(off);
            }
        }
        else if (way == std::ios_base::cur)
        {
            if (testin)
            {
                gbump(off);
            }
            if (testout)
            {
                pbump(off);
            }
        }
        else if (way == std::ios_base::end)
        {
            if (testin)
            {
                setg(eback(), eback(), egptr());
                gbump(egptr() - eback() - off);
            }
            if (testout)
            {
                setp(pbase(), epptr());
                pbump(epptr() - pbase() - off);
            }
        }

        if (testin)
        {
            return (gptr() - eback()); // span_tellg()
        }
        if (testout)
        {
            return (pptr() - pbase()); // span_tellp()
        }

        return -1;
    }


    //- Set position pointer to absolute position
    virtual std::streampos seekpos
    (
        std::streampos pos,
        std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
    )
    {
        return seekoff(pos, std::ios_base::beg, which);
    }

public:

    // Forward Declarations
    class in_base;
    class in_dynamic;
    class out_base;
    class out_dynamic;
};


/*---------------------------------------------------------------------------*\
                     Class memorybuf::in_base Declaration
\*---------------------------------------------------------------------------*/

//- The base input streambuf with memory access
class memorybuf::in_base
:
    public memorybuf
{
protected:

    //- Get sequence of characters from a fixed region
    virtual std::streamsize xsgetn(char* s, std::streamsize n)
    {
        std::streamsize count = 0;
        while (count < n && gptr() < egptr())
        {
            *(s + count++) = *(gptr());
            gbump(1);
        }
        return count;
    }

public:

    // Constructors

        //- Default construct
        in_base() = default;

        //- Construct for character array (can be nullptr) and number of bytes
        in_base(char* s, std::streamsize n)
        {
            resetg(s, n);
        }


    // Member Functions

        // //- Reset get buffer pointer to the beginning of existing span
        // void rewind() { setg(eback(), eback(), egptr()); }

        //- Reset get buffer with character data (can be nullptr) and count
        //  Sets get pointer to the begin.
        void resetg(char* s, std::streamsize n)
        {
            if (s)
            {
                setg(s, s, s + n);
            }
            else
            {
                setg(nullptr, nullptr, nullptr);
            }
        }

        //- The current buffer get position
        std::streamsize span_tellg() const { return (gptr() - eback()); }

        //- The get buffer capacity
        std::streamsize span_capacity() const { return (egptr() - eback()); }

        //- The number of characters remaining in the get area
        std::streamsize span_remaining() const
        {
            return (gptr() < egptr()) ? (egptr() - gptr()) : 0;
        }

        //- The span data (start of input characters)
        char* data_bytes() const { return eback(); }

        //- The span size (number of input characters)
        std::streamsize size_bytes() const { return (egptr() - eback()); }

        //- True if position is within the current input range
        bool in_range(std::streampos pos) const
        {
            return (pos >= 0 && pos < span_capacity());
        }

        //- A string view of the current input region
        std::string_view view() const
        {
            return std::string_view(data_bytes(), size_bytes());
        }

        //- A sub-slice string view of the current input region
        std::string_view view(size_t pos, size_t len) const
        {
            // Restrict intersection count to the current content range
            if (len && in_range(pos))
            {
                if ((len == std::string::npos) || (len > (size_bytes()-pos)))
                {
                    len = (size_bytes()-pos);
                }
            }
            else
            {
                len = 0;  // Ignore out-of-range
            }

            if (len)
            {
                return std::string_view(data_bytes()+pos, len);
            }
            else
            {
                return std::string_view();
            }
        }

        //- Some information about the input buffer position/capacity
        void info(Ostream& os) const
        {
            os  << "get=" << span_tellg() << '/' << span_capacity();
        }
};


/*---------------------------------------------------------------------------*\
                    Class memorybuf::in_dynamic Declaration
\*---------------------------------------------------------------------------*/

//- An input streambuf for memory access, backed by Foam::List<char>
class memorybuf::in_dynamic
:
    public memorybuf::in_base
{
private:

    //- Character storage
    List<char> storage_;

public:

    // Data Types

        //- The dynamic container type when swapping out character content
        typedef ::Foam::DynamicList<char> dynamic_type;


    // Constructors

        //- Default construct - empty
        in_dynamic() = default;

        //- Copy construct from content
        in_dynamic(const char* s, std::streamsize n)
        {
            if (s && n)
            {
                storage_.resize_nocopy(n);
                std::copy(s, (s + n), storage_.data());
            }
            sync_gbuffer();
        }

        //- Move construct from List
        in_dynamic(::Foam::List<char>&& buffer)
        :
            storage_(std::move(buffer))
        {
            sync_gbuffer();
        }

        //- Move construct from DynamicList (added length only)
        template<int SizeMin>
        in_dynamic(::Foam::DynamicList<char,SizeMin>&& buffer)
        {
            storage_.transfer(buffer);  // Implies shrink_to_fit
            sync_gbuffer();
        }


    // Member Functions

        //- Sync get buffer pointers to agree with list dimensions
        //  Sets get pointer to the begin (rewind).
        void sync_gbuffer()
        {
            resetg(storage_.data(), storage_.size());
        }

        //- Reset content (copy)
        void reset(const char* s, std::streamsize n)
        {
            if (s && n)
            {
                storage_.resize_nocopy(n);
                std::copy(s, (s + n), storage_.data());
            }
            else
            {
                storage_.clear();
            }
            sync_gbuffer();
        }

        //- Exchange buffer content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            other.swap(storage_);  // Swap contents
            sync_gbuffer();
        }

        //- Exchange buffer content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            // NB: not storage_.swap(other)! - incorrect slicing
            other.swap(storage_);  // Swap contents: implies shrink_to_fit
            sync_gbuffer();
        }

        //- Reset buffer and return contents.
        //- The list size and capacity are identical
        DynamicList<char> release()
        {
            DynamicList<char> chars(std::move(storage_));
            sync_gbuffer();
            return chars;
        }
};


/*---------------------------------------------------------------------------*\
                     Class memorybuf::out_base Declaration
\*---------------------------------------------------------------------------*/

//- An output streambuf for memory access
class memorybuf::out_base
:
    public memorybuf
{
protected:

    //- Put sequence of characters to a fixed region
    virtual std::streamsize xsputn(const char* s, std::streamsize n)
    {
        std::streamsize count = 0;
        while (count < n && pptr() < epptr())
        {
            *(pptr()) = *(s + count++);
            pbump(1);
        }
        return count;
    }

public:

    // Constructors

        //- Default construct
        out_base() = default;

        //- Construct for character array (can be nullptr) and number of bytes
        out_base(char* s, std::streamsize n)
        {
            resetp(s, n);
        }


    // Member Functions

        // //- Reset put buffer pointer to the beginning of existing span
        // void rewind() { setp(pbase(), epptr()); }

        //- Reset put buffer with character data (can be nullptr) and count
        //  Sets put pointer to the begin.
        inline void resetp(char* s, std::streamsize n)
        {
            if (s)
            {
                // As per (std::ios_base::out && !std::ios_base::ate)
                setp(s, s + n);
                // No treatment for (std::ios_base::out && std::ios_base::ate)
            }
            else
            {
                setp(nullptr, nullptr);
            }
        }

        //- The current buffer put position
        std::streamsize span_tellp() const { return (pptr() - pbase()); }

        //- The put buffer capacity
        std::streamsize span_capacity() const { return (epptr() - pbase()); }

        //- The span data (start of output characters)
        char* data_bytes() const { return pbase(); }

        //- The span size (size of output buffer)
        std::streamsize size_bytes() const { return (pptr() - pbase()); }

        //- True if position is within the current output range
        bool in_range(std::streampos pos) const
        {
            return (pos >= 0 && pos < span_tellp());
        }

        //- A string view of the current output region
        std::string_view view() const
        {
            return std::string_view(data_bytes(), size_bytes());
        }

        //- A sub-slice string view of the current output region
        std::string_view view(size_t pos, size_t len) const
        {
            // Restrict intersection count to the current content range
            if (len && in_range(pos))
            {
                if ((len == std::string::npos) || (len > (size_bytes()-pos)))
                {
                    len = (size_bytes()-pos);
                }
            }
            else
            {
                len = 0;  // Ignore out-of-range
            }

            if (len)
            {
                return std::string_view(data_bytes()+pos, len);
            }
            else
            {
                return std::string_view();
            }
        }

        //- Decrease the put area by 1 or more elements
        void pop_back(int n = 1)
        {
            for (; (n > 0 && (pbase() < pptr())); --n)
            {
                pbump(-1);
            }
        }

        //- Overwrite a sub-slice with character content
        void overwrite
        (
            std::streampos pos,
            const char* data,
            std::streamsize count
        )
        {
            if (data && (count > 0) && in_range(pos))
            {
                // Restrict to intersection with current content
                if (span_tellp() <= pos+std::streampos(count))
                {
                    count = (span_tellp() - pos);
                }
                std::copy(data, (data+count), (data_bytes()+pos));
            }
        }

        //- Overwrite a single character
        void overwrite(std::streampos pos, char c)
        {
            if (c && in_range(pos))
            {
                *(data_bytes()+pos) = c;
            }
        }

        //- Some information about the output buffer position/capacity
        void info(Ostream& os) const
        {
            os  << "put=" << span_tellp() << '/' << span_capacity();
        }
};


/*---------------------------------------------------------------------------*\
                   Class memorybuf::out_dynamic Declaration
\*---------------------------------------------------------------------------*/

//- An output streambuf for memory access, backed by Foam::List<char>
class memorybuf::out_dynamic
:
    public memorybuf::out_base
{
private:

    //- Character storage.
    //  Internally manage like a DynamicList, with its capacity known
    //  from the list size and the addressable size known through the
    //  stream pointers.
    List<char> storage_;

protected:

    //- Handle overflow
    virtual int overflow(int_type c = traits_type::eof())
    {
        if (c != traits_type::eof())
        {
            // Needed more space?
            extend(1);

            *(pptr()) = c;
            pbump(1);
        }
        return c;
    }

    //- Put sequence of characters
    virtual std::streamsize xsputn(const char* s, std::streamsize n)
    {
        // Enough space so that appends work without problem
        extend(n);

        std::streamsize count = 0;
        while (count < n && pptr() < epptr())
        {
            *(pptr()) = *(s + count++);
            pbump(1);
        }

        return count;
    }


public:

    // Data Types

        //- The dynamic container type when swapping out character content
        typedef ::Foam::DynamicList<char> dynamic_type;


    // Constructors

        //- Default construct - no initial reserved number of bytes.
        out_dynamic()
        {
            sync_pbuffer();
        }

        //- Default construct with initial reserved number of bytes.
        out_dynamic(size_t nbytes)
        :
            storage_(label(nbytes))
        {
            sync_pbuffer();
        }

        //- Move construct from List
        out_dynamic(::Foam::List<char>&& buffer)
        :
            storage_(std::move(buffer))
        {
            sync_pbuffer();
        }

        //- Move construct from DynamicList (uses entire capacity)
        template<int SizeMin>
        out_dynamic(::Foam::DynamicList<char,SizeMin>&& buffer)
        {
            buffer.resize(buffer.capacity());  // Use entire space
            storage_.transfer(buffer);
            sync_pbuffer();
        }


    // Member Functions

        //- Normal lower capacity limit.
        //  512 bytes is a bit arbitrary but consistent with std::stringstream
        static constexpr label min_size() noexcept { return 512; }

        //- The largest storage size
        static constexpr label max_size() noexcept
        {
            return UList<char>::max_size();
        }

        //- The 1/2 of max_size() - rounded to power-of-two
        static constexpr label max_size_2() noexcept
        {
            return (1+max_size()/2);
        }

        //- The 1/4 of max_size() - rounded to power-of-two
        static constexpr label max_size_4() noexcept
        {
            return (max_size_2()/2);
        }

        //- Sync put buffer pointers to agree with list dimensions.
        //  Sets put pointer to the begin (rewind).
        void sync_pbuffer()
        {
            resetp(storage_.data(), storage_.size());
        }

        //- Increase capacity (if needed) and adjust buffer pointers.
        //- Applies a min-size and capacity doubling.
        void reserve(const std::streamsize len)
        {
            if (storage_.size() < len)
            {
                const auto cur = span_tellp();  // Output position

                label size = min_size();

                if (size < len)
                {
                    // Increase capacity, but grow more slowly at the largest
                    // sizes. Expect a buffer of char to approach max_size()
                    // more commonly than buffers of other data types.

                    if (storage_.size() <= max_size_4())
                    {
                        // (0 < capacity <= 0.25) : fast growth (2)
                        size = label(2*storage_.size());
                    }
                    else if (storage_.size() <= max_size_2())
                    {
                        // (0.25 < capacity <= 0.5) : slower growth (1.5)
                        size = label((storage_.size()/2)*3);
                    }
                    else if (storage_.size() <= (max_size_2() + max_size_4()))
                    {
                        // (0.5 < capacity <= 0.75) : very slow growth (1.25)
                        size = label((storage_.size()/4)*5);
                    }
                    else
                    {
                        // Already very large - use max (or just len?)
                        size = max_size();
                    }

                    // max(size, len)
                    if (size < len)
                    {
                        size = len;
                    }
                }

                // std::cerr
                //     <<"request:" << len
                //     << " old cap:" << storage_.size()
                //     << " new cap:" << size
                //     << " pos:" << cur << endl;

                storage_.resize_copy(cur, size);
                sync_pbuffer();
                pbump(cur);
            }
        }

        //- Increase capacity for at least this size.
        //- Does not apply min-size or capacity doubling etc.
        void reserve_exact(const std::streamsize len)
        {
            if (storage_.size() < len)
            {
                const auto cur = span_tellp();  // Output position

                storage_.resize_copy(cur, len);
                sync_pbuffer();
                pbump(cur);
            }
        }

        //- Increase (reserve) space for another \c count entries
        void extend(std::streamsize count)
        {
            reserve(span_tellp() + count);
        }

        //- Increase (reserve) space for another \c count entries
        void extend_exact(std::streamsize count)
        {
            reserve_exact(span_tellp() + count);
        }

        //- Clear storage
        void clearStorage()
        {
            storage_.clear();
            sync_pbuffer();
        }

        //- Shrink storage to addressed storage
        void shrink_to_fit()
        {
            const auto cur = span_tellp();  // Addressed length

            storage_.resize(cur);
            sync_pbuffer();
            pbump(cur);
        }

        //- Exchange buffer content and parameter contents, reset positions
        void swap(List<char>& other)
        {
            const auto cur = span_tellp();  // Addressed length
            other.swap(storage_);
            other.resize(cur);              // Truncate to output length
            sync_pbuffer();
        }

        //- Exchange buffer content and parameter contents, reset positions
        template<int SizeMin>
        void swap(DynamicList<char,SizeMin>& other)
        {
            const auto cur = span_tellp();  // Addressed length

            other.resize(other.capacity()); // Use entire space
            other.swap(storage_);  // NB: not storage_.swap(other)
            other.resize(cur);              // Restrict to output length
            sync_pbuffer();
        }

        //- Reset buffer and return contents as a DynamicList.
        //- The list size corresponds to the region of output.
        DynamicList<char> release()
        {
            const auto cur = span_tellp();  // Addressed length
            DynamicList<char> chars(std::move(storage_));
            chars.resize(cur);              // Restrict to output length

            if (chars.empty()) chars.clearStorage();  // Can destroy now
            sync_pbuffer();
            return chars;
        }


    // Housekeeping

        //- Same as shrink_to_fit()
        void shrink() { shrink_to_fit(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
