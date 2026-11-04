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

InClass
    Foam::StringStream

Description
    Input/output from string buffers.

\*---------------------------------------------------------------------------*/

#ifndef Foam_StringStream_H
#define Foam_StringStream_H

#include "ISstream.H"
#include "OSstream.H"
#include <sstream>
#include <string>
#include <string_view>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IStringStream;
class OStringStream;

/*---------------------------------------------------------------------------*\
                        Class IStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Input from string buffer, using a ISstream. Always UNCOMPRESSED.
class IStringStream
:
    public Foam::Detail::StreamAllocator<std::istringstream>,
    public Foam::ISstream
{
    typedef
        Foam::Detail::StreamAllocator<std::istringstream>
        allocator_type;

public:

    // Constructors

        //- Default construct or with specified stream option
        explicit IStringStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Construct from std::string
        explicit IStringStream
        (
            const std::string& s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {
            stream_.str(s);
        }

        //- Construct from char*
        explicit IStringStream
        (
            const char* s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {
            stream_.str(s);
        }

        //- Copy construct, copies content and format
        IStringStream(const IStringStream& ss)
        :
            allocator_type(),
            ISstream(stream_, ss.name(), static_cast<IOstreamOption>(ss))
        {
            stream_.str(ss.stream_.str());
        }


    // Member Functions

        //- Get the current buffer contents as a \c std::string.
        //  \note OpenFOAM-v2606 and earlier: returned a \c Foam::string but
        //   marked with the commment "may change in future" since FEB-2024
        std::string str() const { return stream_.str(); }

        //- Set the string
        void str(const std::string& s) { stream_.str(s); }

        #if defined(__cpp_lib_sstream_from_string_view)
        //- A string_view of the current buffer contents.
        std::string_view view() const { return stream_.view(); }
        #endif

        //- Reset the input buffer and rewind the stream
        virtual void reset(const std::string& s)
        {
            this->str(s);
            this->rewind();
        }

        //- Print stream description
        virtual void print(Ostream& os) const override
        {
            os << "istringstream : buffer =\n" << stream_.str() << '\n';
            IOstream::print(os);
            IOstream::print(os, stream_.rdstate());
        }


    // Member Operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary.
        Istream& operator()() const
        {
            // Could also rewind
            return const_cast<IStringStream&>(*this);
        }


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct empty with given format
        explicit IStringStream(IOstreamOption::streamFormat fmt)
        :
            IStringStream(IOstreamOption(fmt))
        {}

        //- Construct from std::string with given format
        IStringStream
        (
            const std::string& s,
            IOstreamOption::streamFormat fmt
        )
        :
            IStringStream(s, IOstreamOption(fmt))
        {}

        //- Construct from char* with given format
        IStringStream
        (
            const char* s,
            IOstreamOption::streamFormat fmt
        )
        :
            IStringStream(s, IOstreamOption(fmt))
        {}

    #endif /* Foam_IOstream_extras */
};


/*---------------------------------------------------------------------------*\
                        Class OStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Output to string buffer, using a OSstream. Always UNCOMPRESSED.
class OStringStream
:
    public Foam::Detail::StreamAllocator<std::ostringstream>,
    public Foam::OSstream
{
    typedef
        Foam::Detail::StreamAllocator<std::ostringstream>
        allocator_type;

public:

    // Constructors

        //- Default construct or with specified stream option
        explicit OStringStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "output", streamOpt.format(), streamOpt.version())
        {}

        //- Copy construct, copies content and format
        OStringStream(const OStringStream& ss)
        :
            allocator_type(),
            OSstream(stream_, ss.name(), static_cast<IOstreamOption>(ss))
        {
            stream_.str(ss.stream_.str());
        }


    // Member Functions

        //- Get the current buffer contents as a \c std::string.
        //  \note OpenFOAM-v2606 and earlier: returned a \c Foam::string but
        //   marked with the commment "may change in future" since FEB-2024
        std::string str() const { return stream_.str(); }

        //- Set the string
        void str(const std::string& s) { stream_.str(s); }

        #if defined(__cpp_lib_sstream_from_string_view)
        //- A string_view of the current buffer contents.
        std::string_view view() const { return stream_.view(); }
        #endif

        //- Reset the output buffer and rewind the stream
        void reset()
        {
            stream_.str("");  // No other way to reset the end
            this->rewind();
        }

        //- Rewind the output stream
        virtual void rewind()
        {
            stream_.rdbuf()->pubseekpos(0, std::ios_base::out);
        }

        //- Print stream description
        virtual void print(Ostream& os) const override
        {
            os << "ostringstream : buffer =\n" << stream_.str() << '\n';
            IOstream::print(os);
            IOstream::print(os, stream_.rdstate());
        }


    // Older style, without stream option (including 2012 release)
    #ifdef Foam_IOstream_extras

        //- Construct empty with given format
        explicit OStringStream(IOstreamOption::streamFormat fmt)
        :
            OStringStream(IOstreamOption(fmt))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
