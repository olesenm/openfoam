/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2026 OpenCFD Ltd.
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

#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneMesh>
void Foam::fvMeshDistribute::reorderZones
(
    const wordList& zoneNames,
    ZoneMesh& zones
)
{
    zones.clearAddressing();

    typedef typename ZoneMesh::zone_type zone_type;

    // Shift old ones to new position
    PtrList<zone_type>& oldPtrs = zones;
    PtrList<zone_type> newPtrs(zoneNames.size());

    forAll(zones, zonei)
    {
        if (!zones.get(zonei))
        {
            FatalErrorInFunction
                << "Problem with zone[" << zonei << "] " << zones.names()
                << exit(FatalError);
        }
    }

    forAll(zones, zonei)
    {
        // Take ownership
        auto ptr = oldPtrs.release(zonei);
        if (!ptr)
        {
            FatalErrorInFunction
                << "Problem with zone[" << zonei << "] == nullptr" << nl
                << exit(FatalError);
        }

        const label newIndex = zoneNames.find(ptr->name());
        // Check if not found?
        ptr->index() = newIndex;
        newPtrs.set(newIndex, std::move(ptr));
    }

    // Add empty zones for unknown ones
    forAll(newPtrs, i)
    {
        if (!newPtrs.get(i))
        {
            newPtrs.emplace_set(i, zoneNames[i], i, zones);
        }
    }

    // Transfer
    oldPtrs = std::move(newPtrs);
}


template<class GeoField>
void Foam::fvMeshDistribute::printIntFieldInfo(const fvMesh& mesh)
{
    typedef GeometricField
    <
        typename GeoField::value_type,
        fvPatchField,
        volMesh
    > excludeType;

    for (const GeoField& field : mesh.objectRegistry::csorted<GeoField>())
    {
        if (!isA<excludeType>(field))
        {
            Pout<< "Field:" << field.name() << " size:" << field.size()
                //<< " value:" << field
                << endl;
        }
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::printFieldInfo(const fvMesh& mesh)
{
    for (const GeoField& field : mesh.objectRegistry::csorted<GeoField>())
    {
        Pout<< "Field:" << field.name() << " size:" << field.size()
            //<< " value:" << field
            << endl;

        for (const auto& patchFld : field.boundaryField())
        {
            Pout<< "    " << patchFld.patch().index()
                << ' ' << patchFld.patch().name()
                << ' ' << patchFld.type()
                << ' ' << patchFld.size()
                << nl;
        }
    }
}


template<class T, class Mesh>
void Foam::fvMeshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvsPatchField, T>>& bflds
) const
{
    // Save whole boundary field

    typedef FieldField<fvsPatchField, T> ffType;
    typedef GeometricField<T, fvsPatchField, Mesh> fldType;

    const UPtrList<const fldType> flds
    (
        mesh_.objectRegistry::csorted<fldType>()
    );

    bflds.resize_null(flds.size());

    label fieldi = 0;
    for (const fldType& fld : flds)
    {
        bflds.set
        (
            fieldi,
            static_cast<const ffType&>
            (
                fld.boundaryField()
            ).clone(fld.internalField())
        );
        ++fieldi;
    }
}


template<class T, class Mesh>
void Foam::fvMeshDistribute::mapBoundaryFields
(
    const mapPolyMesh& map,
    const PtrList<FieldField<fvsPatchField, T>>& oldBflds
)
{
    // Map boundary field

    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    // typedef FieldField<fvsPatchField, T> ffType;
    typedef GeometricField<T, fvsPatchField, Mesh> fldType;

    UPtrList<fldType> flds
    (
        mesh_.objectRegistry::sorted<fldType>()
    );

    if (flds.size() != oldBflds.size())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    forAll(flds, fieldi)
    {
        auto& bfld = flds[fieldi].boundaryFieldRef();
        const auto& oldBfld = oldBflds[fieldi];

        // Pull from old boundary field into bfld.

        forAll(bfld, patchi)
        {
            fvsPatchField<T>& patchFld = bfld[patchi];
            label facei = patchFld.patch().start();

            forAll(patchFld, i)
            {
                label oldFacei = faceMap[facei++];

                // Find patch and local patch face oldFacei was in.
                forAll(oldPatchStarts, oldPatchi)
                {
                    label oldLocalI = oldFacei - oldPatchStarts[oldPatchi];

                    if (oldLocalI >= 0 && oldLocalI < oldBfld[oldPatchi].size())
                    {
                        patchFld[i] = oldBfld[oldPatchi][oldLocalI];
                    }
                }
            }
        }
    }
}


template<class T>
void Foam::fvMeshDistribute::saveInternalFields
(
    PtrList<Field<T>>& iflds
) const
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> fldType;

    const UPtrList<const fldType> fields
    (
        mesh_.objectRegistry::csorted<fldType>()
    );

    iflds.resize_null(fields.size());

    forAll(fields, fieldi)
    {
        // Bypass DynamicField level, save as Field
        iflds.set
        (
            fieldi,
            static_cast<const Field<T>&>(fields[fieldi]).clone()
        );
    }
}


template<class T>
void Foam::fvMeshDistribute::mapExposedFaces
(
    const mapPolyMesh& map,
    const PtrList<Field<T>>& oldFlds
)
{
    // Set boundary values of exposed internal faces

    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvsPatchField, surfaceMesh> fldType;

    UPtrList<fldType> flds
    (
        mesh_.objectRegistry::sorted<fldType>()
    );

    if (flds.size() != oldFlds.size())
    {
        FatalErrorInFunction
            << "problem"
            << abort(FatalError);
    }


    forAll(flds, fieldi)
    {
        auto& fld = flds[fieldi];
        const auto& oldInternal = oldFlds[fieldi];

        const bool oriented = fld.is_oriented();

        auto& bfld = fld.boundaryFieldRef();


        // Pull from old internal field into bfld.

        forAll(bfld, patchi)
        {
            fvsPatchField<T>& patchFld = bfld[patchi];

            forAll(patchFld, i)
            {
                const label faceI = patchFld.patch().start()+i;

                label oldFaceI = faceMap[faceI];

                if (oldFaceI < oldInternal.size())
                {
                    patchFld[i] = oldInternal[oldFaceI];

                    if (oriented && map.flipFaceFlux().contains(faceI))
                    {
                        patchFld[i] = flipOp()(patchFld[i]);
                    }
                }
            }
        }
    }
}


template<class GeoField, class PatchFieldType>
void Foam::fvMeshDistribute::initPatchFields
(
    const typename GeoField::value_type& initVal
)
{
    // Init patch fields of certain type
    // - field order is irrelevant

    for (GeoField& fld : mesh_.objectRegistry::objects<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        forAll(bfld, patchi)
        {
            if (isA<PatchFieldType>(bfld[patchi]))
            {
                bfld[patchi] == initVal;
            }
        }
    }
}


//template<class GeoField>
//void Foam::fvMeshDistribute::correctBoundaryConditions()
//{
//    // CorrectBoundaryConditions patch fields of certain type
//
//    for (GeoField& fld : mesh_.objectRegistry::sorted<GeoField>())
//    {
//        fld.correctBoundaryConditions();
//    }
//}


template<class GeoField>
void Foam::fvMeshDistribute::getFieldNames
(
    const fvMesh& mesh,
    HashTable<wordList>& allFieldNames,
    const word& excludeType,
    const bool syncPar
)
{
    wordList& list = allFieldNames(GeoField::typeName);
    list = mesh.sortedNames<GeoField>();

    // Use communicator associated with mesh
    const auto comm = mesh.comm();
    const int uniqTag = (UPstream::msgType() + 1907);

    if (!excludeType.empty())
    {
        const wordList& excludeList =
            allFieldNames.lookup(excludeType, wordList::null());

        if (!excludeList.empty())
        {
            label count = 0;
            forAll(list, i)
            {
                if (!excludeList.contains(list[i]))
                {
                    if (count != i)
                    {
                        list[count] = std::move(list[i]);
                    }
                    ++count;
                }
            }
            list.resize(count);
        }
    }


    // Check all procs have same names
    const auto& localNames = list;

    if (!syncPar || !UPstream::is_parallel(comm))
    {
        // No checking or not parallel
    }
    else if (UPstream::master(comm))
    {
        bool hasError = false;

        for (int proci : UPstream::subProcs(comm))
        {
            IPstream is(proci, uniqTag, comm);
            wordList procNames(is);

            if (procNames != localNames)
            {
                if (!hasError)
                {
                    hasError = true;

                    FatalErrorInFunction
                        << "When checking for equal " << GeoField::typeName
                        << " :" << nl
                        << "processor0 has:" << localNames << nl;
                }

                FatalError
                    << "processor" << proci << " has:" << procNames << nl;
            }
        }

        if (hasError)
        {
            FatalError
                << GeoField::typeName
                << " need to be synchronised on all ranks." << endl
                << exit(FatalError);
        }
    }
    else if (UPstream::is_subrank(comm))
    {
        OPstream os(UPstream::masterNo(), uniqTag, comm);
        os << localNames;
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::sendFields
(
    const label domain,
    const HashTable<wordList>& allFieldNames,
    const fvMeshSubset& subsetter,
    Ostream& toNbr
)
{
    // Send fields. Note order supplied so we can receive in exactly the same
    // order.
    // Note that field gets written as entry in dictionary so we
    // can construct from subdictionary.
    // (since otherwise the reading as-a-dictionary mixes up entries from
    // consecutive fields)
    // The dictionary constructed is:
    //  volScalarField
    //  {
    //      p {internalField ..; boundaryField ..;}
    //      k {internalField ..; boundaryField ..;}
    //  }
    //  volVectorField
    //  {
    //      U {internalField ...  }
    //  }

    // volVectorField {U {internalField ..; boundaryField ..;}}

    const wordList& fieldNames =
        allFieldNames.lookup(GeoField::typeName, wordList::null());

    toNbr << GeoField::typeName << token::NL << token::BEGIN_BLOCK << token::NL;

    for (const word& fieldName : fieldNames)
    {
        if (debug)
        {
            Pout<< "Subsetting " << GeoField::typeName
                << " field " << fieldName
                << " for domain:" << domain << endl;
        }

        // Send all fieldNames. This has to be exactly the same set as is
        // being received!
        const GeoField& fld =
            subsetter.baseMesh().lookupObject<GeoField>(fieldName);

        // Note: use subsetter to get sub field. Override default behaviour
        //       to warn for unset fields since they will be reset later on
        tmp<GeoField> tsubfld = subsetter.interpolate(fld, true);

        toNbr
            << fieldName << token::NL << token::BEGIN_BLOCK
            << tsubfld
            << token::NL << token::END_BLOCK << token::NL;
    }
    toNbr << token::END_BLOCK << token::NL;
}


template<class GeoField>
void Foam::fvMeshDistribute::receiveFields
(
    const label domain,
    const HashTable<wordList>& allFieldNames,
    fvMesh& mesh,
    PtrList<GeoField>& fields,
    const dictionary& allFieldsDict
)
{
    // Opposite of sendFields

    const wordList& fieldNames =
        allFieldNames.lookup(GeoField::typeName, wordList::null());

    const dictionary& fieldDicts =
        allFieldsDict.subDict(GeoField::typeName);


    if (debug)
    {
        Pout<< "Receiving:" << GeoField::typeName
            << " fields:" << fieldNames
            << " from domain:" << domain << endl;
    }

    fields.resize(fieldNames.size());

    label fieldi = 0;
    for (const word& fieldName : fieldNames)
    {
        if (debug)
        {
            Pout<< "Constructing type:"  << GeoField::typeName
                << " field:" << fieldName
                << " from domain:" << domain << endl;
        }

        fields.set
        (
            fieldi++,
            new GeoField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                fieldDicts.subDict(fieldName)
            )
        );
    }
}


// ************************************************************************* //
