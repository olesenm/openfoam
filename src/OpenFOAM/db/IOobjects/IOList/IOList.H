/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::IOList

Description
    A List of objects of type \<T\> with automated input and output.

SourceFiles
    IOList.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOList_H
#define Foam_IOList_H

#include "List.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOList
:
    public regIOobject,
    public List<T>
{
    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool readIOcontents();

        //- Read the content size.
        //  Return the size if read, -1 otherwise
        label readIOsize();

public:

    //- The underlying content type
    typedef List<T> content_type;

    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Default copy construct
        IOList(const IOList&) = default;

        //- Construct from IOobject. If not read, will be zero size.
        explicit IOList(const IOobject& io);

        //- Construct from IOobject. If not read, will be zero size.
        IOList(const IOobject& io, Foam::zero);

        //- Construct from IOobject. If not read, use given size.
        IOList(const IOobject& io, const label len);

        //- Construct from IOobject. If not read, use given size/value.
        IOList(const IOobject& io, const label len, const T& val);

        //- Construct from IOobject. If not read, use given size/zero-value.
        IOList(const IOobject& io, const label len, Foam::zero);

        //- Construct from IOobject. If not read, copy content.
        IOList(const IOobject& io, const UList<T>& content);

        //- Construct by transferring the content
        IOList(const IOobject& io, List<T>&& content);


    // Factory Methods

        //- Read and return content size, -1 if not read.
        //  The IOobject is never registered
        static label readContentsSize(const IOobject& io);

        //- Read and return contents. The IOobject is never registered
        static List<T> readContents(const IOobject& io);

        //- Write contents. The IOobject is never registered.
        //  Uses IOListRef internally.
        static void writeContents(const IOobject& io, const UList<T>& content);


    //- Destructor
    virtual ~IOList() = default;


    // Member Functions

        //- The writeData method for regIOobject write operation
        virtual bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy or move assignment of entries
        using List<T>::operator=;

        //- Copy assignment of entries
        void operator=(const IOList<T>& rhs)
        {
            List<T>::operator=(rhs);
        }

        //- Move assignment of entries
        void operator=(IOList<T>&& rhs)
        {
            List<T>::operator=(std::move(static_cast<List<T>&>(rhs)));
        }
};


/*---------------------------------------------------------------------------*\
                          Class IOListRef Declaration
\*---------------------------------------------------------------------------*/

//- A IOList wrapper for writing external List data
template<class T>
class IOListRef
:
    public regIOobject
{
    // Private Data

        //- Reference to the external content
        UList<T> contentRef_;


public:

    //- The underlying content type
    typedef List<T> content_type;


    //- Type is identical to IOList
    virtual const word& type() const
    {
        return IOList<T>::typeName;
    }


    // Generated Methods

        //- No default construct
        IOListRef() = delete;

        //- No copy construct
        IOListRef(const IOListRef&) = delete;

        //- No copy assignment
        void operator=(const IOListRef&) = delete;


    // Constructors

        //- Construct from IOobject and const data reference
        IOListRef(const IOobject& io, const UList<T>& content);


    //- Destructor
    virtual ~IOListRef() = default;


    // Member Functions

        //- Allow cast to const content
        operator const UList<T>&() const
        {
            return contentRef_;
        }

        //- The writeData method for regIOobject write operation
        //  Fatal if content is not set
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOList.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
