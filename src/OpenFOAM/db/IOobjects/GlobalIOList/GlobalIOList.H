/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::GlobalIOList

Description
    IOList with global data (so optionally read from master)

SourceFiles
    GlobalIOList.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_GlobalIOList_H
#define Foam_GlobalIOList_H

#include "IOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GlobalIOList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GlobalIOList
:
    public regIOobject,
    public List<Type>
{
public:

    //- The underlying content type
    typedef List<Type> content_type;

    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Default copy construct
        GlobalIOList(const GlobalIOList&) = default;

        //- Construct from IOobject. If not read, will be zero size.
        explicit GlobalIOList(const IOobject& io);

        //- Construct from IOobject. If not read, will be zero size.
        GlobalIOList(const IOobject& io, Foam::zero);

        //- Construct from IOobject. If not read, use given size.
        GlobalIOList(const IOobject& io, const label len);

        //- Construct from IOobject. If not read, use given size/value.
        GlobalIOList(const IOobject& io, const label len, const Type& val);

        //- Construct from IOobject. If not read, use given size/zero-value.
        GlobalIOList(const IOobject& io, const label len, Foam::zero);

        //- Construct from IOobject. If not read, copy content
        GlobalIOList(const IOobject& io, const UList<Type>& content);

        //- Construct by transferring the List content
        GlobalIOList(const IOobject& io, List<Type>&& content);


    // Factory Methods

        //- Read and return contents. The IOobject is never registered
        static List<Type> readContents(const IOobject& io);


    //- Destructor
    virtual ~GlobalIOList() = default;


    // Member Functions

        //- This object is global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //- either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }

        //- The readData method for regIOobject read operation
        virtual bool readData(Istream& is);

        //- The writeData method for regIOobject write operation
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy or move assignment of entries
        using List<Type>::operator=;

        //- Copy assignment of entries
        void operator=(const GlobalIOList<Type>& rhs)
        {
            List<Type>::operator=(rhs);
        }

        //- Move assignment of entries
        void operator=(GlobalIOList<Type>&& rhs)
        {
            List<Type>::operator=(std::move(static_cast<List<Type>&>(rhs)));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for GlobalIOList
template<class T>
struct is_globalIOobject<GlobalIOList<T>> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GlobalIOList.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
