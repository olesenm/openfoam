/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "SpanStream.H"
#include "wordList.H"
#include "IOstreams.H"
#include "argList.H"

#include <charconv>
#include <cctype>
#include <cstdio>
#include <limits>
#include <iomanip>

using namespace Foam;

Ostream& printString(Ostream& os, const char* first, const char* last)
{
    os << '"';
    for (; first != last; (void)++first)
    {
        os << *first;
    }
    os << '"';

    return os;
}


Ostream& printView(Ostream& os, const char* first, const char* last)
{
    char buf[4];
    os << label(last-first) << '(';

    for (; first != last; (void)++first)
    {
        const char c = *first;

        if (isprint(c))
        {
            os << c;
        }
        else if (c == '\t')
        {
            os << "\\t";
        }
        else if (c == '\n')
        {
            os << "\\n";
        }
        else
        {
            ::snprintf(buf, 4, "%02X", c);
            os << "\\x" << buf;
        }
    }
    os << ')';

    return os;
}


Ostream& printView(Ostream& os, std::string_view s)
{
    return printView(os, s.begin(), s.end());
}


Ostream& printView(Ostream& os, const UList<char>& list)
{
    return printView(os, list.begin(), list.end());
}


Ostream& writeList(Ostream& os, const UList<char>& list)
{
    return printView(os, list);
}


Ostream& toString(Ostream& os, const UList<char>& list)
{
    return printString(os, list.begin(), list.end());
}


Ostream& toString(Ostream& os, std::string_view s)
{
    return printString(os, s.begin(), s.end());
}


template<class BufType>
void printInfo(const BufType& buf)
{
    Info<< nl << "=========================" << endl;
    buf.print(Info);
    Info<< "addr: " << Foam::name(buf.view().data()) << nl;
    toString(Info, buf.view());
    Info<< nl << "=========================" << endl;
}


// Return a left-padded integer as "word"
template<class IntType>
std::string leftpadded(IntType val, char fillch = ' ')
{
    std::string buf;
    buf.resize((std::numeric_limits<IntType>::digits10+1), fillch);

    auto first = (buf.data());
    auto last  = (buf.data() + buf.size());

    auto result = std::to_chars(first, last, val);

    if (result.ec == std::errc{})
    {
        auto* iter = result.ptr;
        int count = std::distance(iter, last);

        std::cout << "did not fill: " << count << " chars\n";

        // With two spaces before comments
        if (count > 0) { *iter++ = ' '; --count; }
        if (count > 0) { *iter++ = ' '; --count; }
        for (char c = (count >= 2 ? '/' : ' '); count > 0; --count)
        {
            *iter++ = c;
        }
    }

    return buf;
}


template<class IntType>
void leftpad(std::ostream& os, IntType val, char fillch = ' ')
{
    // set fill char and width
    os.setf(std::ios_base::left, std::ios_base::adjustfield);
    fillch = os.fill(fillch);
    os.width(std::numeric_limits<IntType>::digits10+1);
    os  << val;
    // restore fill char
    os.fill(fillch);
}


template<class IntType>
void rightpad(std::ostream& os, IntType val, char fillch = ' ')
{
    // set fill char and width
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    fillch = os.fill(fillch);
    os.width(std::numeric_limits<IntType>::digits10+1);
    os  << val;
    // restore fill char
    os.fill(fillch);
}


// Left-padded value with trailing comment slashes
template<class IntType>
void leftpad(ocharstream& os, IntType val)
{
    const auto beg = os.tellp();
    os << val;
    int count = (std::numeric_limits<IntType>::digits10+1) - (os.tellp() - beg);

    // With two spaces before comments
    if (count > 0) { os << ' '; --count; }
    if (count > 0) { os << ' '; --count; }
    for (const char c = (count >= 2 ? '/' : ' '); count > 0; --count)
    {
        os << c;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addBoolOption("fake-zerosize", "Fake overwriting with zero data");
    argList::addBoolOption("dict-format", "Format as dictionary entry");

    #include "setRootCase.H"

    const bool optFakeZerosize = args.found("fake-zerosize");
    const bool isDictFormat = args.found("dict-format");

    // const constexpr int width = (std::numeric_limits<label>::digits10+1);

    // experiment with to_chars instead of streaming
    {
        // Some value
        label val(1234);

        auto fixed = leftpadded(val);
        Info<< "leftpadded " << val << " : " << fixed << nl;
    }

    ocharstream labelbuf;
    labelbuf.reserve_exact(32);

    // Some value
    labelbuf.rewind();
    rightpad(labelbuf, label(10));

    printInfo(labelbuf);

    OCharStream obuf;
    obuf.reserve_exact(48);

    printInfo(obuf);

    obuf.push_back('>');
    obuf.append(" string_view ");
    obuf.push_back('<');
    printInfo(obuf);

    obuf.pop_back(8);
    printInfo(obuf);

    obuf.pop_back(100);
    printInfo(obuf);


    // Fill with some content
    for (int i = 0; i < 26; ++i)
    {
        obuf<< char('A' + i);
    }

    // Change letter 'O' to '_'
    if (auto i = obuf.view().find('O'); i != std::string::npos)
    {
        obuf.overwrite(i, '_');
    }

    // append and push_back some content
    obuf.append(5, '<');
    obuf.push_back('#');
    obuf.append(5, '>');

    printInfo(obuf);

    obuf.pop_back(8);
    printInfo(obuf);

    obuf.overwrite(4, labelbuf.view());
    printInfo(obuf);

    obuf.overwrite(20, "####");
    printInfo(obuf);

    Info<< "operation Ignored..." << nl;
    obuf.overwrite(45, "????");
    printInfo(obuf);

    // Update with new value
    {
        labelbuf.rewind();
        rightpad(labelbuf, label(200), '.');
        obuf.overwrite(4, labelbuf.view());

        printInfo(obuf);
    }

    // With yet another value (non-fixed width)
    {
        labelbuf.rewind();
        labelbuf << label(15);
        obuf.overwrite(4, labelbuf.view());

        printInfo(obuf);
    }

    // Slightly harder test
    {
        std::string chars(26, '?');

        for (int i = 0; i < 26; ++i)
        {
            chars[i] = char('A' + i);
        }

        auto& os = obuf;
        os.rewind();

        const word procName("processor0");

        // Write as primitiveEntry or commented content
        // // constexpr bool isDictFormat = true;

        // if constexpr (isDictFormat)
        if (isDictFormat)
        {
            // Like writeKeyword() with compoundToken
            os << nl << procName << ' ' << word("List<char>") << nl;
        }
        else
        {
            // Human-readable comments
            os << nl << "// " << procName << nl;
        }

        // This is the code we want to have, but assume we don't know
        // the size or data beforehand.
        //
        // if (str && len > 0)
        // {
        //     // Special treatment for char data (binary I/O only)
        //     const auto oldFmt = os.format(IOstreamOption::BINARY);
        //
        //     os << label(len) << nl;
        //     os.write(str, len);
        //     os << nl;
        //
        //     os.format(oldFmt);
        // }
        // else
        // {
        //     os << label(0) << nl;
        // }

        // Position before writing the label
        const auto labelBegin = os.tellp();

        // Replace: os << label(len) << nl;
        // with a fixed-length version
        {
            labelbuf.rewind();
            rightpad(labelbuf, 0);

            os.append(labelbuf.view());
            os << nl;
        }

        constexpr bool testUnknown = true;

        label dataCount = 0;

        if constexpr (testUnknown)
        {
            // Pretend we don't know the number of characters a priori

            const auto oldFmt = os.format(IOstreamOption::BINARY);

            const auto lineNumber = os.lineNumber();

            // count is unknown but irrelevant for serial
            os.beginRawWrite(0);

            // Position before raw binary data
            const auto dataBegin = os.tellp();

            // Some type of output, streaming etc
            os.writeRaw(chars.data(), chars.size());

            // How many chars of binary data written?
            dataCount = (os.tellp() - dataBegin);

            os.endRawWrite();
            os.lineNumber() = lineNumber;
            os << nl;

            os.format(oldFmt);
        }
        else
        {
            // If we had all data collected a priori

            dataCount = chars.size();

            const auto oldFmt = os.format(IOstreamOption::BINARY);

            if (dataCount > 0)
            {
                os.write(chars.data(), chars.size());
                os << nl;
            }
            os.format(oldFmt);
        }

        if (optFakeZerosize)
        {
            dataCount = 0;  // fake zero-size
        }

        printInfo(os);

        // Update the data count with the correct value

        if (dataCount > 0)
        {
            labelbuf.rewind();
            leftpad(labelbuf, label(dataCount));

            os.overwrite(labelBegin, labelbuf.view());
        }
        else
        {
            os.seek(int64_t(labelBegin)-1);

            // if constexpr (isDictFormat)
            if (isDictFormat)
            {
                os << ' ' << label(0);
            }
            else
            {
                os << nl << label(0) << nl;
            }
        }

        // if constexpr (isDictFormat)
        if (isDictFormat)
        {
            os.endEntry();
        }

        printInfo(os);

        Info<< "view: " << os.view(4, 8) << nl;
        Info<< "view: " << os.view(32) << nl;
        // Ignores out-of-range
        Info<< "view: " << os.view(1000) << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
