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

\*---------------------------------------------------------------------------*/

#include "collatedFileOperation.H"
#include "decomposedBlockData.H"
#include "Time.H"
#include "OCharStream.H"
#include "UPstreamFile.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// A fixed-width 0-padded integer
template<class IntType>
void zeropadded(std::ostream& os, IntType val)
{
    // set fill char and width
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    char fillch = os.fill('0');
    os.width(std::numeric_limits<IntType>::digits10+1);
    os  << val;
    // restore fill char
    os.fill(fillch);
}

//
// Some fancy rewriting of the header content to include
// <block.start>, <block.count>, <block.sizes> entries
//
// Uses a fixed-width integer format for the <block.start> entry,
// which allows in-place rewriting.
//
template<class IntType>
void rewriteBlockHeaderInfo
(
    Foam::OCharStream& header,
    const Foam::UList<IntType>& blockSizes,
    int32_t blockCount = 0  // Fallback (when blockSizes is empty)
)
{
    using namespace Foam;

    // An int32 (or even smaller) is large enough for the
    // size of the header content, which is the offset to the first block
    typedef int32_t headerOffsetType;

    // Get the block count from blockSizes
    if (!blockSizes.empty())
    {
        blockCount = static_cast<int32_t>(blockSizes.size());
    }

    if
    (
        const auto paste = header.view().rfind('}');
        paste != std::string::npos
    )
    {
        // Keep everything in ASCII
        const auto oldFmt = header.format(IOstreamOption::ASCII);

        // Fixed-width label entry
        Foam::ocharstream labelbuf;
        labelbuf.reserve_exact(32);

        // Everything trailing after the last '}' from 'FoamFile {}'
        std::string trailing(header.view().substr(paste));
        header.seek(paste);

        // <block.start> entry
        header.append("    block.start ");

        // Position before writing the label for <block.start>
        const auto labelBegin = header.tellp();

        // Fixed-length integer 00000000
        {
            labelbuf.rewind();
            zeropadded(labelbuf, headerOffsetType(0));

            header.append(labelbuf.view());
            header.endEntry();
        }

        // <block.count> entry
        if (blockCount > 0)
        {
            header.append("    block.count ");
            header << blockCount;
            header.endEntry();
        }

        // <block.sizes> entries - flat list on a single line
        if (!blockSizes.empty())
        {
            header.append("    block.sizes\n");
            blockSizes.writeList(header);
            header.endEntry();
        }

        // Reattach old content
        header.append(trailing);

        // Update (overwrite) the <block.start> value
        {
            labelbuf.rewind();
            zeropadded(labelbuf, headerOffsetType(header.view().size()));

            header.overwrite(labelBegin, labelbuf.view());
        }

        // Restore format
        header.format(oldFmt);
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileOperations::collatedFileOperation::writeObject_mpiio
(
    const fileName& pathName,
    const regIOobject& io,
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    const Time& tm = io.time();
    const fileName& inst = io.instance();

    // ----------------------------------------------------------------
    // Compile-time selection of header styles:
    // - plain : no extra information
    // - start : include the starting offset to the first block
    // - sizes : also include the block sizes

    enum class headerTypes { PLAIN, START, SIZES };
    constexpr auto headerType = headerTypes::START;

    // ----------------------------------------------------------------
    // Compile-time selection of footer styles:
    // - none :
    // - agglom : include the footer in the [N-1] block

    enum class footerTypes { NONE, AGGLOM };
    constexpr auto footerType = footerTypes::AGGLOM;

    // ----------------------------------------------------------------

    if
    (
        (inst.isAbsolute() || !tm.processorCase())
     || (io.global() || io.globalObject())
     || (!UPstream::is_parallel(comm_))
    )
    {
        FatalErrorInFunction
            << "Should not have been called for any of these conditions:"
            << " - isAbsolute" << nl
            << " - not processorCase" << nl
            << " - global or globalObject" << nl
            << " - not parallel" << nl
            << abort(FatalError);

        return false;
    }
    else if (!UPstream::File::supported())
    {
        FatalErrorInFunction
            << "Should not have been called without MPI-IO support" << nl
            << abort(FatalError);

        return false;
    }
    else
    {
        // Stream to memory and then write with MPI/IO

        const auto blocki = UPstream::myProcNo(comm_);
        const auto nblock = UPstream::nProcs(comm_);

        // Fixed-width label entry
        ocharstream labelbuf;
        labelbuf.reserve_exact(32);

        // The content buffer
        OCharStream os(streamOpt);

        // The file-level header (when needed separately)
        OCharStream header(streamOpt);

        // The file-level footer.
        // Format everywhere so that its size is known and defined without
        // additional communication, regardless where/how it is written.

        OCharStream footer(streamOpt);
        if (footerType != footerTypes::NONE)
        {
            IOobject::writeEndDivider(footer);
        }
        const int64_t footerSize(footer.view().size());


        // Header for the outer container
        // - separately or combine into the first block
        if (UPstream::master(comm_))
        {
            OCharStream& dest =
            (
                // Only the SIZES style needs a separate header
                (headerType == headerTypes::SIZES)
              ? header
              : os
            );

            // Need binary for the outer container
            const auto oldFmt = dest.format(IOstreamOption::BINARY);

            decomposedBlockData::writeHeader
            (
                dest,
                streamOpt,
                io
            );

            dest.format(oldFmt);

            // Some header styles can already be finalized
            if (headerType == headerTypes::START)
            {
                rewriteBlockHeaderInfo
                (
                    dest,
                    UList<int>(),   // No size information
                    nblock          // <block.count> information
                );
            }
        }

        // Write as commented content, or as primitiveEntry
        constexpr bool isDictFormat = false;

        // The introducers for the content
        {
            const word procName("processor" + Foam::name(blocki));
            if constexpr (isDictFormat)
            {
                // Like writeKeyword()
                os << nl << procName << nl;
            }
            else
            {
                // Human-readable comments
                os << nl << "// " << procName << nl;
            }
        }

        // Begin of block content  LABEL(...)

        // Position before writing the label
        const auto labelBegin = os.tellp();

        // Replace: os << label(len) << nl;
        // with a fixed-length version
        {
            labelbuf.rewind();
            zeropadded(labelbuf, label(0));

            os.append(labelbuf.view());
            os << nl;
        }

        const auto lineNumber = os.lineNumber();

        // Begin binary blob
        {
            const auto oldFmt = os.format(IOstreamOption::BINARY);

            // count is unknown but irrelevant for serial stream
            os.beginRawWrite(0);

            os.format(oldFmt);
        }

        // Position of binary blob - after the '(' begin
        const auto blobBegin = os.tellp();

        // Generate content
        bool ok = true;

        // Block 0 gets a FoamFile header, without comment banner
        if (UPstream::master(comm_))
        {
            const bool oldBanner = IOobject::bannerEnabled(false);

            ok = ok && io.writeHeader(os);

            IOobject::bannerEnabled(oldBanner);
        }

        if (writeOnProc)
        {
            ok = ok && io.writeData(os);
        }

        // How many chars of binary data were written?
        const int64_t blobCount(os.tellp() - blobBegin);

        // Finalize the binary blob - closing ')'
        os.endRawWrite();
        os.lineNumber(lineNumber);
        os << nl;

        // Update the size information for the binary blob
        if (blobCount > 0)
        {
            labelbuf.rewind();
            zeropadded(labelbuf, label(blobCount));

            os.overwrite(labelBegin, labelbuf.view());
        }
        else
        {
            // Seek with begin-1 to also overwrite newline with space
            os.seek(int64_t(labelBegin)-1);

            if constexpr (isDictFormat)
            {
                os << ' ' << label(0);
            }
            else
            {
                os << nl << label(0) << nl;
            }
        }

        if constexpr (isDictFormat)
        {
            os.endEntry();
        }


        // ========================================
        // All content now exists, still need the following:
        // - the per-rank file begin offset
        // - the total output size

        // The output start for the rank
        int64_t blockOffset = 0;

        // The total size (including any header and footer)
        int64_t totalSize = 0;

        if (headerType != headerTypes::SIZES)
        {
            // Header without a master overview of information.
            // This is the simplest to determine
            // => Exscan and broadcast combination

            const int64_t headerSize(header.view().size());

            int64_t localSize(os.view().size());

            // Include any separate header into the size of block [0]
            if (UPstream::master(comm_))
            {
                localSize += headerSize;
            }

            // Like globalObject::calculate(), except that we
            // would discard the size anyhow.
            blockOffset = UPstream::find_offset(localSize, comm_);

            // Adjust begin offset of block [0] for the presence of any
            // separate header. No need to adjust any of the other blocks
            // since the her header size was added into the local size of
            // block [0] prior to the Exscan_sum

            if (UPstream::master(comm_))
            {
                blockOffset = headerSize;

                // Is this a possible corner case?
                // Collated with one rank and separated header??
                // Unlikely: if (nblock == 1) { localSize -= headerSize; }
            }

            // The last block provides the total size
            if (blocki == (nblock-1))
            {
                totalSize = (blockOffset + localSize);

                if (footerType == footerTypes::AGGLOM)
                {
                    totalSize += footerSize;
                }
            }

            // Broadcast totalSize with [N-1] as the root
            UPstream::broadcast(&totalSize, 1, comm_, nblock-1);
        }
        else
        {
            // With <block.sizes> in the header, master needs an overview
            // but nobody else does.
            //
            // == gather, calculate, scatter + broadcast.
            // This will still be better than an Allgather.

            // The sizes (without header or footer!) gathered onto the master

            List<int64_t> sizes
            (
                UPstream::listGatherValues<int64_t>(os.view().size(), comm_)
            );

            if (UPstream::master(comm_))
            {
                // Update header with <block.sizes>, ... information.
                rewriteBlockHeaderInfo(header, sizes);

                // Now rewrite sizes to be the begin offsets.

                // Collective output begins after any separate header
                totalSize = header.view().size();

                for (auto& val : sizes)
                {
                    const auto count = val;
                    val = totalSize;  // The begin offset
                    totalSize += count;
                }

                // If the footer will be included in the last block:
                if (footerType == footerTypes::AGGLOM)
                {
                    totalSize += footerSize;
                }
            }

            // The file begin offset for the local block
            blockOffset = UPstream::listScatterValues(sizes, comm_);

            // Broadcast totalSize with master as the root
            UPstream::broadcast(&totalSize, 1, comm_, UPstream::masterNo());
        }

        // Agglomerate the footer into the last block
        if (footerType == footerTypes::AGGLOM && blocki == (nblock-1))
        {
            os.extend_exact(footer.view().size());
            os.append(footer.view());
        }


        // ========================================
        // All preparation is done and we have the following:
        // - header : the separate header (only for SIZES format)
        // - os : the content
        // - blockOffset : the output start for the rank
        // - totalSize : size of content and any header/footer
        //
        // The local blockSize to be written is known directly from the
        // stream view.
        // ========================================

        // Make directory, open file
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        const fileName fName
        (
            objectPath(io, decomposedBlockData::typeName)
        );


        // Using mkDir not Foam::mkDir
        mkDir(fName.path());

        // Write file contents
        {
            UPstream::File file;

            file.open_write(comm_, fName, IOstreamOption::ATOMIC);

            // Any separate header to be written?
            if
            (
                const auto view = header.view();
                (!view.empty() && UPstream::master(comm_))
            )
            {
                // Master only write
                file.write_at(0, view);
            }

            // Collective write
            file.write_at_all(blockOffset, os.view());

            // The file size (already includes header and footer sizes)
            file.set_size(totalSize);
            file.close();
        }

        return ok;
    }
}


// ************************************************************************* //
