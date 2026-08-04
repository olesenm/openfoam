/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
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
    Foam::CompactIOList

Description
    A List of objects of type \<T\> with automated input and output using
    a compact storage. Behaves like IOList except when binary output in
    case it writes a CompactListList.

    Useful for lists of small sublists e.g. faceList, cellList.

SourceFiles
    CompactIOList.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_CompactIOList_H
#define Foam_CompactIOList_H

#include "IOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class CompactIOList;

template<class T> Istream& operator>>(Istream&, CompactIOList<T>&);
template<class T> Ostream& operator<<(Ostream&, const CompactIOList<T>&);

/*---------------------------------------------------------------------------*\
                        Class CompactIOList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class CompactIOList
:
    public regIOobject,
    public List<T>
{
    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        //- Reads according to the header type.
        //  Return true if read.
        bool readIOcontents();

        //- Read the content size.
        //  Return the size if read, -1 otherwise
        label readIOsize();

        //- Has too many elements in it?
        bool overflows() const;


public:

    //- The underlying content type
    typedef List<T> content_type;

    //- Runtime type information
    TypeName("CompactList");


    // Constructors

        //- Default copy construct
        CompactIOList(const CompactIOList&) = default;

        //- Construct from IOobject. If not read, will be zero size.
        explicit CompactIOList(const IOobject& io);

        //- Construct from IOobject. If not read, will be zero size.
        CompactIOList(const IOobject& io, Foam::zero);

        //- Construct from IOobject. If not read, use given size.
        CompactIOList(const IOobject& io, const label len);

        //- Construct from IOobject. If not read, use given size/value.
        CompactIOList(const IOobject& io, const label len, const T& value);

        //- Construct from IOobject. If not read, use given size/zero-value.
        CompactIOList(const IOobject& io, const label len, Foam::zero);

        //- Construct from IOobject and List content
        CompactIOList(const IOobject& io, const UList<T>& content);

        //- Construct by transferring the List content
        CompactIOList(const IOobject& io, List<T>&& content);


    // Factory Methods

        //- Read and return content size, -1 if not read.
        //  The IOobject will not be registered
        static label readContentsSize(const IOobject& io);

        //- Read and return contents. The IOobject will not be registered
        static List<T> readContents(const IOobject& io);


    //- Destructor
    virtual ~CompactIOList() = default;


    // Member Functions

        //- Read as offsets/packed-values and transcribe into *this
        Istream& readCompact(Istream& is);

        //- Write as offsets/packed-values
        Ostream& writeCompact(Ostream& os) const;

        //- Write using stream options. Checks for overflow in binary
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;

        //- Write as plain or compact content (depends on stream format)
        virtual bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy or move assignment of entries
        using List<T>::operator=;

        //- Copy assignment of entries
        void operator=(const CompactIOList<T>& rhs)
        {
            List<T>::operator=(rhs);
        }

        //- Move assignment of entries
        void operator=(CompactIOList<T>&& rhs)
        {
            List<T>::operator=(std::move(static_cast<List<T>&>(rhs)));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactIOList.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
