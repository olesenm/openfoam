/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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

#include "decomposedBlockData.H"
#include "Fstream.H"
#include "IPstream.H"
#include "OPstream.H"
#include "SpanStream.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::decomposedBlockData::isCollatedType(const word& objectType)
{
    return (objectType == decomposedBlockData::typeName);
}


bool Foam::decomposedBlockData::isCollatedType(const IOobject& io)
{
    // same: return io.isHeaderClass(decomposedBlockData::typeName);
    // same: return isCollatedType(io.headerClassName());
    return io.isHeaderClass<decomposedBlockData>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decomposedBlockData::decomposedBlockData
(
    const label comm,
    const IOobject& io,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(comm)
{
    // Temporary warning
    if (readOpt() == IOobjectOption::READ_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with READ_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        read();
    }
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    const label blocki,
    const char* str,
    const size_t len
)
{
    // Offset to the beginning of this output
    // This should generally be OK for non-compressed streams
    // (eg, std::ofstream)

    std::streamoff blockOffset = os.stdStream().tellp();

    const word procName("processor" + Foam::name(blocki));

    // Write as primitiveEntry or commented content
    constexpr bool isDictFormat = false;

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

    // Data
    if (str && len > 0)
    {
        // Special treatment for char data (binary I/O only)
        const auto oldFmt = os.format(IOstreamOption::BINARY);

        os << label(len) << nl;
        os.write(str, len);
        os << nl;

        os.format(oldFmt);
    }
    else
    {
        os << label(0) << nl;
    }

    if constexpr (isDictFormat)
    {
        os.endEntry();
    }

    return blockOffset;
}


bool Foam::decomposedBlockData::readBlockEntry
(
    Istream& is,
    List<char>& charData
)
{
    // Handle any of these:

    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;
    // 3.  processorN  NCHARS (...) ;

    is.fatalCheck(FUNCTION_NAME);
    token tok(is);
    is.fatalCheck(FUNCTION_NAME);

    // Dictionary format?
    // - has a keyword, which is not a compound token:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (!isDictFormat && tok.good())
    {
        is.putBack(tok);
    }
    charData.readList(is);

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return true;
}


bool Foam::decomposedBlockData::skipBlockEntry(Istream& is)
{
    // As per readBlockEntry but seeks instead of reading.
    // Internals like charList::readList - ie, always binary

    // Handle any of these:
    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;
    // 3.  processorN  NCHARS (...) ;

    if (!is.good()) return false;
    token tok(is);
    if (!is.good()) return false;

    // Dictionary format?
    // - has a keyword, which is not a compound token:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (isDictFormat)
    {
        is >> tok;
        if (!is.good()) return false;
    }


    bool handled = false;

    // Like charList::readList
    if (tok.isCompound())
    {
        handled = true;
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..) or just a plain '0'

        const label len = tok.labelToken();

        // Special treatment for char data (binary I/O only)
        const auto oldFmt = is.format(IOstreamOption::BINARY);

        if (len)
        {
            // read(...) includes surrounding start/end delimiters.

            // Note: nullptr to ignore instead of reading
            is.read(nullptr, std::streamsize(len));
        }
        is.format(oldFmt);

        handled = true;
    }
    else
    {
        // Incorrect token
        return false;
    }

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return handled;
}


Foam::label Foam::decomposedBlockData::getNumBlocks
(
    Istream& is,
    const label maxNumBlocks
)
{
    label nBlocks = 0;

    // Handle OpenFOAM header if it is the first entry
    if (is.good())
    {
        token tok(is);

        if (is.good() && tok.isWord("FoamFile"))
        {
            dictionary headerDict(is);  // Read sub-dictionary content

            if (headerDict.readIfPresent("version", tok))
            {
                is.version(tok);
            }

            word formatName;
            if (headerDict.readIfPresent("format", formatName))
            {
                is.format(formatName);
            }

            //// Obtain number of blocks directly
            ///  This may not be reliable...
            //if (headerDict.readIfPresent("blocks", nBlocks))
            //{
            //    return nBlocks;
            //}
        }
        else if (tok.good())
        {
            is.putBack(tok);
        }
    }

    while (is.good() && skipBlockEntry(is))
    {
        ++nBlocks;

        if (maxNumBlocks == nBlocks)
        {
            break;
        }
    }

    return nBlocks;
}


bool Foam::decomposedBlockData::hasBlock(Istream& is, const label blockNumber)
{
    return
    (
        blockNumber >= 0
     && (blockNumber < getNumBlocks(is, blockNumber+1))
    );
}


std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    IOstreamOption streamOptData,
    const regIOobject& io,
    const label blocki,
    const bool withLocalHeader
)
{
    // Serialize content to write
    DynamicList<char> serialized;
    {
        OCharStream buf(streamOptData);
        buf.reserve(4*1024);  // Start with a slightly larger buffer

        bool ok = true;

        // Generate FoamFile header on master, without comment banner
        if (withLocalHeader)
        {
            const bool old = IOobject::bannerEnabled(false);

            ok = io.writeHeader(buf);

            IOobject::bannerEnabled(old);
        }

        // Serialize the output
        ok = ok && io.writeData(buf);

        if (!ok)
        {
            return std::streamoff(-1);
        }

        // Take ownership of serialized content
        serialized = buf.release();
    }

    return decomposedBlockData::writeBlockEntry(os, blocki, serialized);
}


Foam::autoPtr<Foam::ISstream>
Foam::decomposedBlockData::readBlock
(
    const label blocki,
    ISstream& is,
    IOobject& headerIO
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlock:"
            << " stream:" << is.name() << " attempt to read block " << blocki
            << endl;
    }

    // The character input stream for the specified block
    autoPtr<ISstream> blockIsPtr;

    // Extracted header information
    IOstreamOption streamOptData;
    unsigned labelWidth = is.labelByteSize();
    unsigned scalarWidth = is.scalarByteSize();

    // Read master for header
    List<char> data;
    decomposedBlockData::readBlockEntry(is, data);

    if (blocki == 0)  // ie, UPstream::masterNo()
    {
        blockIsPtr.reset(new ICharStream(std::move(data)));
        blockIsPtr->name() = is.name();

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*blockIsPtr))
            {
                FatalIOErrorInFunction(*blockIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }
    }
    else
    {
        {
            // Read header from first block,
            // without advancing the stream position
            ISpanStream headerStream(data);
            if (!headerIO.readHeader(headerStream))
            {
                FatalIOErrorInFunction(headerStream)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
            streamOptData = static_cast<IOstreamOption>(headerStream);
            labelWidth = headerStream.labelByteSize();
            scalarWidth = headerStream.scalarByteSize();
        }

        // Skip intermediate blocks
        for (label i = 1; i < blocki; ++i)
        {
            decomposedBlockData::skipBlockEntry(is);
        }

        // Read the block of interest
        decomposedBlockData::readBlockEntry(is, data);

        blockIsPtr.reset(new ICharStream(std::move(data)));
        blockIsPtr->name() = is.name();

        // Apply stream settings
        {
            auto& iss = blockIsPtr();
            iss.format(streamOptData.format());
            iss.version(streamOptData.version());
            iss.setLabelByteSize(labelWidth);
            iss.setScalarByteSize(scalarWidth);
        }
    }

    return blockIsPtr;
}


bool Foam::decomposedBlockData::readBlocks
(
    const label comm,
    autoPtr<ISstream>& isPtr,
    List<char>& localData,
    const UPstream::commsTypes  /* unused */
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "<null>")
            << " non-blocking comm:" << comm << endl;
    }

    // Read data on master and transmit. Always non-blocking
    bool ok = false;

    // The send buffers
    List<List<char>> procBuffers;

    // Some unique tag for this read/write grouping (as extra safety)
    const int messageTag = (UPstream::msgType() + 256);

    const auto startOfRequests = UPstream::nRequests();

    if (UPstream::master(comm))
    {
        auto& is = isPtr();
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, localData);

        // Read proc data and setup non-blocking sends
        procBuffers.resize(UPstream::nProcs(comm));
        for (const int proci : UPstream::subProcs(comm))
        {
            auto& slot = procBuffers[proci];

            decomposedBlockData::readBlockEntry(is, slot);

            // Send content (non-blocking)
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                slot.cdata_bytes(),
                slot.size_bytes(),
                messageTag,
                comm
            );
        }

        ok = is.good();
    }
    else if (UPstream::is_subrank(comm))
    {
        List<char>& slot = localData;

        // Probe for the message size
        auto [fromProci, num_bytes] =
            UPstream::probeMessage(UPstream::masterNo(), messageTag, comm);

        if (debug)
        {
            Pout<< "probed to receive " << num_bytes
                << " from " << fromProci << endl;
        }

        slot.resize_nocopy(num_bytes);

        // Receive content (can also be zero-sized)
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::masterNo(),
            slot.data_bytes(),
            slot.size_bytes(),
            messageTag,
            comm
        );
    }

    UPstream::waitRequests(startOfRequests);
    procBuffers.clear();

    // Sync the status
    Pstream::broadcast(ok, comm);

    return ok;
}


Foam::autoPtr<Foam::ISstream> Foam::decomposedBlockData::readBlocks
(
    const label comm,
    const fileName& fName,
    autoPtr<ISstream>& isPtr,
    IOobject& headerIO,
    const UPstream::commsTypes  /* unused */
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "<null>")
            << " non-blocking" << endl;
    }

    // Read data on master and transmit. Always non-blocking
    bool ok = false;
    List<char> localData;
    List<List<char>> procBuffers;
    autoPtr<ISstream> blockIsPtr;

    // Some unique tag for this read/write/probe grouping
    const int messageTag = (UPstream::msgType() + 256);

    const auto startOfRequests = UPstream::nRequests();

    if (UPstream::master(comm))
    {
        auto& is = *isPtr;
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, localData);

        // Move block data into a stream
        blockIsPtr.reset(new ICharStream(std::move(localData)));
        blockIsPtr->name() = fName;

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*blockIsPtr))
            {
                FatalIOErrorInFunction(*blockIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }

        // Read proc data and setup non-blocking sends
        procBuffers.resize(UPstream::nProcs(comm));
        for (const int proci : UPstream::subProcs(comm))
        {
            auto& slot = procBuffers[proci];

            decomposedBlockData::readBlockEntry(is, slot);

            // Send content - non-blocking mode
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                slot.cdata_bytes(),
                slot.size_bytes(),
                messageTag,
                comm
            );
        }

        ok = is.good();
    }
    else if (UPstream::is_subrank(comm))
    {
        List<char>& slot = localData;

        // Probe for the message size
        auto [fromProci, num_bytes] =
            UPstream::probeMessage(UPstream::masterNo(), messageTag, comm);

        if (debug)
        {
            Pout<< "probed to receive " << num_bytes
                << " from " << fromProci << endl;
        }

        slot.resize_nocopy(num_bytes);

        // Receive content (can also be zero-sized)
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::masterNo(),
            slot.data_bytes(),
            slot.size_bytes(),
            messageTag,
            comm
        );
    }

    UPstream::waitRequests(startOfRequests);
    procBuffers.clear();

    if (UPstream::is_subrank(comm))
    {
        // Move block data into a stream
        blockIsPtr.reset(new ICharStream(std::move(localData)));
        blockIsPtr->name() = fName;
    }

    // Broadcast master header info,
    // set stream properties from blockIsPtr on master

    int verValue(0);
    int fmtValue(0);
    unsigned labelWidth(0);
    unsigned scalarWidth(0);
    word headerName(headerIO.name());

    // The stream characteristics
    //
    // unsigned formatSizes
    // (
    //     ((static_cast<unsigned>(iss.format()) & 0xFF) << 16)
    //   | ((iss.labelByteSize() & 0xFF) << 8)
    //   | ((iss.scalarByteSize() & 0xFF))
    // );

    if (UPstream::master(comm))
    {
        auto& iss = blockIsPtr();
        verValue = iss.version().canonical();
        fmtValue = static_cast<int>(iss.format());
        labelWidth = iss.labelByteSize();
        scalarWidth = iss.scalarByteSize();
    }

    Pstream::broadcasts
    (
        comm,
        verValue,
        fmtValue,
        labelWidth,
        scalarWidth,
        headerName,
        headerIO.headerClassName(),
        headerIO.note()
        // Unneeded: headerIO.instance()
        // Unneeded: headerIO.local()
    );

    if (blockIsPtr)
    {
        auto& iss = *blockIsPtr;
        iss.version(IOstreamOption::versionNumber::canonical(verValue));
        iss.format(IOstreamOption::streamFormat(fmtValue));
        iss.setLabelByteSize(labelWidth);
        iss.setScalarByteSize(scalarWidth);
    }

    headerIO.rename(headerName);

    if (debug)
    {
        Info<< "reading ok:" << ok << endl;
    }

    return blockIsPtr;
}


void Foam::decomposedBlockData::gatherProcData
(
    const label comm,
    const UList<char>& localData,
    const labelUList& recvSizes,

    const labelRange& whichProcs,

    List<int>& sliceOffsets,
    DynamicList<char>& recvData,
    const UPstream::commsTypes commsType
)
{
    const label myRank = UPstream::myProcNo(comm);
    const label nProcs = UPstream::nProcs(comm);

    // Some unique tag for this read/write grouping (as extra safety)
    const int messageTag = (UPstream::msgType() + 256);

    int nSendBytes = 0;
    recvData.clear();

    // On master, calculate sizing/offsets and resize the recv buffer.
    // Do not need sliceSizes when nonBlocking
    List<int> sliceSizes;
    if (UPstream::master(comm))
    {
        sliceSizes.resize_nocopy(nProcs);
        sliceSizes = 0;
        sliceOffsets.resize_nocopy(nProcs+1);
        sliceOffsets = 0;

        int totalSize = 0;
        for (const label proci : whichProcs)
        {
            const auto nRecvBytes = static_cast<int>(recvSizes[proci]);

            sliceOffsets[proci] = totalSize;
            totalSize += nRecvBytes;

            sliceSizes[proci] = nRecvBytes;
        }

        // One beyond the end of the range
        const label endProci = whichProcs.end_value();

        sliceOffsets[endProci] = totalSize;
        recvData.resize_nocopy(totalSize);
    }
    else if (whichProcs.contains(myRank) && !localData.empty())
    {
        // Note: UPstream::gather limited to int
        nSendBytes = static_cast<int>(localData.size_bytes());
    }

    if (UPstream::commsTypes::nonBlocking == commsType)
    {
        if (UPstream::master(comm))
        {
            for (const label proci : whichProcs)
            {
                SubList<char> procSlice
                (
                    recvData,
                    sliceOffsets[proci+1]-sliceOffsets[proci],
                    sliceOffsets[proci]
                );

                if (procSlice.empty())
                {
                    continue;
                }
                else if (proci == UPstream::masterNo())
                {
                    // No self-communication, although masterNo is normally
                    // not contained in whichProcs range anyhow.
                    std::copy
                    (
                        localData.cbegin(),
                        localData.cbegin(procSlice.size()),
                        procSlice.begin()
                    );
                }
                else
                {
                    // Receive non-zero content
                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        procSlice.data_bytes(),
                        procSlice.size_bytes(),
                        messageTag,
                        comm
                    );
                }
            }
        }
        else if (whichProcs.contains(myRank) && !localData.empty())
        {
            // Send non-zero content
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::masterNo(),
                localData.cdata_bytes(),
                localData.size_bytes(),
                messageTag,
                comm
            );
        }

        // Waiting is done by the caller
    }
    else
    {
        // This is MPI_Gatherv() !! - but this path is unlikely to be used

        UPstream::mpiGatherv
        (
            localData.cdata(),
            nSendBytes,

            recvData.data(),
            sliceSizes,
            sliceOffsets,
            comm
        );
    }
}


bool Foam::decomposedBlockData::writeBlocks
(
    const label comm,
    autoPtr<OSstream>& osPtr,
    List<std::streamoff>& blockOffset,
    const UList<char>& localData,

    const labelUList& recvSizes,
    const UList<std::string_view>& procData,

    const UPstream::commsTypes commsType,
    const bool syncReturnState
)
{
    const label nProcs = UPstream::nProcs(comm);

    bool ok = true;

    // Recovery of blockOffset is optional
    if (UPstream::master(comm) && notNull(blockOffset))
    {
        blockOffset.resize(nProcs);
    }

    // Max proc data size to be received
    label maxNonLocalSize = 0;
    if (UPstream::master(comm) && procData.empty())
    {
        for (label proci = 1; proci < nProcs; ++proci)
        {
            maxNonLocalSize = Foam::max(maxNonLocalSize, recvSizes[proci]);
        }
    }

    if (debug)
    {
        Pout<< " stream:" << (osPtr ? osPtr->name() : "<null>")
            << " data:" << localData.size()
            << " proc-data:" << procData.size()
            << " max-size:" << maxNonLocalSize
            << " " << UPstream::commsTypeNames[commsType] << endl;
    }

    if (procData.size())
    {
        // --------
        // With pre-gathered proc data
        // --------

        if (UPstream::master(comm))
        {
            auto& os = osPtr();

            std::streamoff currOffset =
                decomposedBlockData::writeBlockEntry
                (
                    os,
                    UPstream::masterNo(),
                    localData
                );

            if (blockOffset.size() > UPstream::masterNo())
            {
                blockOffset[UPstream::masterNo()] = currOffset;
            }

            // Write all pre-gathered proc data.
            for (label proci = 1; proci < nProcs; ++proci)
            {
                currOffset =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        procData[proci]
                    );

                if (blockOffset.size() > proci)
                {
                    blockOffset[proci] = currOffset;
                }
            }

            ok = os.good();
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        // --------
        // Gather/write each rank, one at a time.
        // Note: This is often associated with maxMasterFileBufferSize == 0
        // --------

        // Some unique tag for this read/write grouping (as extra safety)
        const int messageTag = (UPstream::msgType() + 256);

        if (UPstream::master(comm))
        {
            auto& os = osPtr();

            std::streamoff currOffset =
                decomposedBlockData::writeBlockEntry
                (
                    os,
                    UPstream::masterNo(),
                    localData
                );

            if (blockOffset.size() > UPstream::masterNo())
            {
                blockOffset[UPstream::masterNo()] = currOffset;
            }

            // Could discard/recycle localData on master
            // (if we had taken ownership...)

            DynamicList<char> recvData(maxNonLocalSize);
            for (label proci = 1; proci < nProcs; ++proci)
            {
                recvData.resize_nocopy(recvSizes[proci]);

                if (!recvData.empty())
                {
                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        proci,
                        recvData.data_bytes(),
                        recvData.size_bytes(),
                        messageTag,
                        comm
                    );
                }

                currOffset =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        recvData
                    );

                if (blockOffset.size() > proci)
                {
                    blockOffset[proci] = currOffset;
                }
            }

            ok = os.good();
        }
        else if (UPstream::is_subrank(comm) && !localData.empty())
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                localData.cdata_bytes(),
                localData.size_bytes(),
                messageTag,
                comm
            );
        }
    }
    else
    {
        // --------
        // Gather/write ranks, packing together several smaller gathers
        // into a single buffer space
        // --------

        DynamicList<char> recvData;
        List<int> recvOffsets;  // Offsets into recvData

        // Offsets of combined ranks for communication.
        // Never includes master rank (handled separately)
        labelList procOffsets(nProcs, Foam::zero{});

        // Max combined data to be received (master only)
        label maxRecvCount = 0;

        if (UPstream::master(comm))
        {
            // Find out how many ranks can be received into
            // maxMasterFileBufferSize and the corresponding schedule

            off_t maxBufferSize
            (
                fileOperations::masterUncollatedFileOperation::
                maxMasterFileBufferSize
            );

            // Buffer must fit the largest off-processor size
            if (maxBufferSize < off_t(maxNonLocalSize))
            {
                maxBufferSize = off_t(maxNonLocalSize);
            }

            // Max combined proc data size to be received
            off_t maxCollected = 0;

            for (label proci = 1, nChunks = 0; proci < nProcs; /*nil*/)
            {
                procOffsets[nChunks] = proci;

                // At least one proc, regardless of maxBufferSize.
                // Also handles the corner case when the first proc has
                // size 0, but the next one is too large.

                for
                (
                    off_t total = 0;
                    (
                        proci < nProcs
                     && (!total || (total + recvSizes[proci] < maxBufferSize))
                    );
                    ++proci
                )
                {
                    total += recvSizes[proci];

                    if (maxCollected < total)
                    {
                        maxCollected = total;
                    }
                }

                procOffsets[++nChunks] = proci;
            }

            maxRecvCount = static_cast<label>(maxCollected);
        }

        if (debug && UPstream::master(comm))
        {
            OStringStream ranges;

            for (label nChunks = 1; nChunks < nProcs; ++nChunks)
            {
                const labelRange whichProcs
                (
                    procOffsets[nChunks-1],
                    procOffsets[nChunks]-procOffsets[nChunks-1]
                );

                if (whichProcs.start() >= nProcs || whichProcs.size() <= 0)
                {
                    break;
                }

                ranges << ' ' << whichProcs.min() << '-' << whichProcs.max();
            }

            Pout<< " write-schedule:" << ranges.str().c_str() << endl;
        }


        // Same schedule to be known by everyone
        UPstream::broadcast(procOffsets.data(), procOffsets.size(), comm);

        recvData.resize_nocopy(label(maxRecvCount));  // (master only)

        if (UPstream::master(comm))
        {
            auto& os = osPtr();

            std::streamoff currOffset =
                decomposedBlockData::writeBlockEntry
                (
                    os,
                    UPstream::masterNo(),
                    localData
                );

            if (blockOffset.size() > UPstream::masterNo())
            {
                blockOffset[UPstream::masterNo()] = currOffset;
            }
        }

        for (label nChunks = 1; nChunks < nProcs; ++nChunks)
        {
            const labelRange whichProcs
            (
                procOffsets[nChunks-1],
                procOffsets[nChunks]-procOffsets[nChunks-1]
            );

            if (whichProcs.start() >= nProcs || whichProcs.size() <= 0)
            {
                break;
            }

            const auto startOfRequests = UPstream::nRequests();

            // Setup non-blocking send/recv or MPI_Gatherv
            // - uses (UPstream::msgType()+256)
            gatherProcData
            (
                comm,
                localData,
                recvSizes,

                whichProcs,

                recvOffsets,
                recvData,
                commsType  // ie, blocking or non-blocking
            );

            // For sanity checks
            // const auto endOfRequests = UPstream::nRequests();

            if (UPstream::master(comm))
            {
                auto& os = osPtr();

                // Write received data
                auto currRequest = startOfRequests;
                for (const label proci : whichProcs)
                {
                    SubList<char> procSlice
                    (
                        recvData,
                        recvOffsets[proci+1]-recvOffsets[proci],
                        recvOffsets[proci]
                    );

                    if
                    (
                        (UPstream::commsTypes::nonBlocking == commsType)
                     && (proci != UPstream::masterNo())
                     && !procSlice.empty()
                    )
                    {
                        UPstream::waitRequest(currRequest);
                        ++currRequest;
                    }

                    std::streamoff currOffset =
                        decomposedBlockData::writeBlockEntry
                        (
                            os,
                            proci,
                            procSlice
                        );

                    if (blockOffset.size() > proci)
                    {
                        blockOffset[proci] = currOffset;
                    }
                }

                ok = os.good();
            }

            UPstream::waitRequests(startOfRequests);
        }
    }

    if (syncReturnState)
    {
        //- Enable to get synchronised error checking.
        //  Ensures that all procs are as slow as the master
        //  (which does all the writing)
        Pstream::broadcast(ok, comm);
    }

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<ISstream> isPtr;
    fileName objPath(fileHandler().filePath(false, *this, word::null));
    if (UPstream::master(comm_))
    {
        isPtr.reset(new IFstream(objPath));
        IOobject::readHeader(*isPtr);
    }

    return readBlocks(comm_, isPtr, contentData_, commsType_);
}


bool Foam::decomposedBlockData::writeData(Ostream& os) const
{
    IOobject io(*this);
    IOstreamOption streamOpt(os);

    int verValue(0);
    int fmtValue(0);
    // Unneeded: word masterName(name());
    fileName masterLocation(instance()/db().dbDir()/local());

    // Re-read my own data to find out the header information
    if (UPstream::master(comm_))
    {
        ISpanStream headerStream(contentData_);
        io.readHeader(headerStream);

        verValue = headerStream.version().canonical();
        fmtValue = static_cast<int>(headerStream.format());
    }

    // Broadcast header information
    Pstream::broadcasts
    (
        comm_,
        verValue,
        fmtValue,
        // Unneeded: masterName
        io.headerClassName(),
        io.note(),
        // Unneeded: io.instance()
        // Unneeded: io.local()
        masterLocation
    );

    streamOpt.version(IOstreamOption::versionNumber::canonical(verValue));
    streamOpt.format(IOstreamOption::streamFormat(fmtValue));

    if (UPstream::is_subrank(comm_))
    {
        decomposedBlockData::writeHeader
        (
            os,
            streamOpt,  // streamOpt for data
            io.headerClassName(),
            io.note(),
            masterLocation,
            name(),
            dictionary()
        );
    }

    // Write the character data
    if (isA<OSstream>(os))
    {
        // Serial stream - can output characters directly
        os.writeRaw(contentData_.cdata(), contentData_.size_bytes());
    }
    else
    {
        // Other cases are less fortunate, and no std::string_view
        os.writeQuoted(contentData_.cdata(), contentData_.size_bytes(), false);
    }

    if (UPstream::is_subrank(comm_))
    {
        IOobject::writeEndDivider(os);
    }

    return os.good();
}


bool Foam::decomposedBlockData::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm_))
    {
        // Note: always write binary. These are strings so readable anyway.
        //       They have already be tokenised on the sending side.

        osPtr.reset(new OFstream(objectPath(), IOstreamOption::BINARY));

        // Update meta-data for current state
        const_cast<regIOobject&>
        (
            static_cast<const regIOobject&>(*this)
        ).updateMetaData();

        decomposedBlockData::writeHeader
        (
            *osPtr,
            streamOpt,  // streamOpt for data
            static_cast<const IOobject&>(*this)
        );
    }

    const labelList recvSizes
    (
        UPstream::listGatherValues<label>(contentData_.size(), comm_)
    );

    List<std::streamoff> blockOffsets;  // Optional
    return writeBlocks
    (
        comm_,
        osPtr,
        blockOffsets,
        contentData_,
        recvSizes,
        UList<std::string_view>(),  // dummy proc data (nothing pre-gathered)
        commsType_
    );
}


// ************************************************************************* //
