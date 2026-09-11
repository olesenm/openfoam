/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::NullObject

Description
    Singleton null-object class and instance.

    Its contents occupy enough space to also be reinterpreted
    as another class with a null pointer or zero long for its first
    member, with additional zero parameters for safe casting to List etc.

SourceFiles
    nullObject.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_nullObject_H
#define Foam_nullObject_H

#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                         Class NullObject Declaration
\*---------------------------------------------------------------------------*/

class NullObject
{
    //- A %union of zero data types
    union zeros
    {
        void* ptr;
        unsigned long val;
    };


    // Private Data

        //- The zero data content
        zeros data_[8];


    // Constructors

        //- Private constructor for singleton only,
        //- using zero initialized union content.
        constexpr NullObject() : data_{} {}

        //- No copy construct (only use as singleton)
        NullObject(const NullObject&) = delete;

        //- No copy assignment (only use as singleton)
        void operator=(const NullObject&) = delete;


public:

    // Static Data

        //- A unique null object
        static const NullObject nullObject;


    // Member Functions

        //- A nullptr pointer content
        constexpr const void* pointer() const noexcept { return data_[0].ptr; }

        //- For zero-valued integer content
        static constexpr unsigned long value() noexcept { return 0; }

        //- Empty is true (no elements)
        static constexpr bool empty() noexcept { return true; }

        //- Size is zero (no elements)
        static constexpr label size() noexcept { return 0; }

        //- No-op method (for HashTable replacement)
        constexpr const NullObject& toc() const noexcept { return *this; }

        //- No-op method (for HashTable replacement)
        constexpr const NullObject& sortedToc() const noexcept { return *this; }


    // Member Operators

        //- Swallow assignment (cf, std::ignore)
        template<class T>
        constexpr const NullObject& operator=(const T&) const noexcept
        {
            return *this;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Globals

//- Pointer to the unique nullObject
extern const NullObject* nullObjectPtr;


// IOstream Operators

//- Read from Istream consumes no content, does not change NullObject
inline Istream& operator>>(Istream& is, const NullObject&) noexcept
{
    return is;
}

//- Write to Ostream emits no content
inline Ostream& operator<<(Ostream& os, const NullObject&) noexcept
{
    return os;
}


// Global Functions

//- Const pointer (of type T) to the nullObject
template<class T>
inline const T* NullObjectPtr() noexcept
{
    return reinterpret_cast<const T*>(nullObjectPtr);
}

//- Non-const pointer (of type T) to the nullObject.
//- Only use when nothing will be written into it!
template<class T>
inline T* NullObjectPtr_constCast() noexcept
{
    return reinterpret_cast<T*>(const_cast<NullObject*>(nullObjectPtr));
}


//- Const reference (of type T) to the nullObject
template<class T>
inline const T& NullObjectRef() noexcept
{
    return *reinterpret_cast<const T*>(nullObjectPtr);
}

//- Non-const reference (of type T) to the nullObject
//- Only use when nothing will be written into it!
template<class T>
inline T& NullObjectRef_constCast() noexcept
{
    return *reinterpret_cast<T*>(const_cast<NullObject*>(nullObjectPtr));
}


//- True if ptr is a pointer (of type T) to the nullObject
template<class T>
inline bool isNull(const T* ptr) noexcept
{
    return ptr == NullObjectPtr<T>();
}

//- True if obj is a reference (of type T) to the nullObject
template<class T>
inline bool isNull(const T& obj) noexcept
{
    return &obj == NullObjectPtr<T>();
}


//- True if ptr is not a pointer (of type T) to the nullObject
template<class T>
inline bool notNull(const T* ptr) noexcept
{
    return ptr != NullObjectPtr<T>();
}

//- True if obj is not a reference (of type T) to the nullObject
template<class T>
inline bool notNull(const T& obj) noexcept
{
    return &obj != NullObjectPtr<T>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
