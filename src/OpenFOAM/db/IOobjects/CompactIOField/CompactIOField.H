/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::CompactIOField

Description
    A Field of objects of type \<T\> with automated input and output using
    a compact storage. Behaves like IOField except when binary output in
    case it writes a CompactListList.

    Useful for fields of small subfields e.g. in lagrangian

SourceFiles
    CompactIOField.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_CompactIOField_H
#define Foam_CompactIOField_H

#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class CompactIOField;

template<class T> Istream& operator>>(Istream&, CompactIOField<T>&);
template<class T> Ostream& operator<<(Ostream&, const CompactIOField<T>&);

/*---------------------------------------------------------------------------*\
                        Class CompactIOField Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class CompactIOField
:
    public regIOobject,
    public Field<T>
{
    // Private Member Functions

        //- Read if IOobject flags set and 'on-proc' is true.
        //- Reads according to the header type.
        //  Return true if read (only accurate when readOnProc == true).
        bool readIOcontents(bool readOnProc = true);

        //- Read the content size.
        //  Return the size if read, -1 otherwise
        label readIOsize(bool readOnProc = true);

        //- Has too many elements in it?
        bool overflows() const;

public:

    //- The underlying content type
    typedef Field<T> content_type;

    //- Runtime type information
    TypeName("FieldField");


    // Constructors

        //- Default copy construct
        CompactIOField(const CompactIOField&) = default;

        //- Construct from IOobject. If not read, will be zero size.
        explicit CompactIOField(const IOobject& io);

        //- Construct from IOobject, with local processor conditional reading
        CompactIOField(const IOobject& io, const bool readOnProc);

        //- Construct from IOobject. If not read, will be zero size.
        CompactIOField(const IOobject& io, Foam::zero);

        //- Construct from IOobject. If not read, use given size.
        CompactIOField(const IOobject& io, const label len);

        //- Construct from IOobject. If not read, use given size/value.
        CompactIOField(const IOobject& io, const label len, const T& value);

        //- Construct from IOobject. If not read, use given size/zero-value.
        CompactIOField(const IOobject& io, const label len, Foam::zero);

        //- Construct from IOobject and a List/Field content
        CompactIOField(const IOobject& io, const UList<T>& content);

        //- Construct by transferring the Field contents
        CompactIOField(const IOobject& io, Field<T>&& content);


    // Factory Methods

        //- Read and return content size, -1 if not read.
        //  The IOobject is never registered
        static label readContentsSize(const IOobject& io);

        //- Read and return contents. The IOobject is never registered
        static Field<T> readContents(const IOobject& io);


    //- Destructor
    virtual ~CompactIOField() = default;


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
        using Field<T>::operator=;

        //- Copy assignment of entries
        void operator=(const CompactIOField<T>& rhs)
        {
            Field<T>::operator=(rhs);
        }

        //- Move assignment of entries
        void operator=(CompactIOField<T>&& rhs)
        {
            Field<T>::operator=(std::move(static_cast<Field<T>&>(rhs)));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactIOField.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
