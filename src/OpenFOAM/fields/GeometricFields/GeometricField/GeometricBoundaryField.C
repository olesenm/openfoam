/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
    Copyright (C) 2016-2026 OpenCFD Ltd.
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

#include "GeometricBoundaryField.H"
#include "globalMeshData.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
template<class CheckPatchFieldType>
bool Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::checkConsistency
(
    const scalar tol,
    const bool exitIfBad
) const
{
    auto& bfld = this->constCast();

    if (!bfld.size())
    {
        return true;
    }

    if (debug&2)
    {
        PoutInFunction
            << " Checking boundary consistency for field "
            << bfld[0].internalField().name() << endl;
    }

    // Store old values and states
    List<Field<Type>> oldFields(bfld.size());
    boolList oldUpdated(bfld.size());
    boolList oldManipulated(bfld.size());

    label nEvaluated(0);

    for (auto& pfld : bfld)
    {
        if (isA<CheckPatchFieldType>(pfld))
        {
            const label patchi = pfld.patch().index();
            oldFields[patchi] = pfld;
            oldUpdated[patchi] = pfld.updated();
            oldManipulated[patchi] = pfld.manipulatedMatrix();
            ++nEvaluated;
        }
    }

    if (!nEvaluated) return true;  // Early termination

    // Re-evaluate
    {
        const auto startOfRequests = UPstream::nRequests();

        nEvaluated = 0;

        for (auto& pfld : bfld)
        {
            if (isA<CheckPatchFieldType>(pfld))
            {
                pfld.initEvaluate(UPstream::commsTypes::nonBlocking);
                ++nEvaluated;
            }
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : bfld)
        {
            if (isA<CheckPatchFieldType>(pfld))
            {
                pfld.evaluate(UPstream::commsTypes::nonBlocking);
                if (--nEvaluated == 0) break;  // Early termination
            }
        }
    }


    // Check
    bool allOk(true);

    for (auto& pfld : bfld)
    {
        if (isA<CheckPatchFieldType>(pfld))
        {
            const label patchi = pfld.patch().index();
            auto& oldPfld = oldFields[patchi];

            bool localOk(true);

            if (allOk)
            {
                // Only check for first failed patch
                forAll(pfld, facei)
                {
                    if (tol < Foam::mag(pfld[facei]-oldPfld[facei]))
                    {
                        allOk = false;
                        localOk = false;
                        break;
                    }
                }
            }

            if (!localOk)
            {
                // Raise warning or error
                OSstream& err =
                (
                    exitIfBad
                  ? FatalErrorInFunction
                  : WarningInFunction
                );

                err << "Field "
                    << pfld.internalField().name()
                    << " is not evaluated?"
                    << " On patch " << pfld.patch().name()
                    << " type " << pfld.type()
                    << " : average of field = "
                    << average(oldPfld)
                    << ". Average of evaluated field = "
                    << average(pfld)
                    << ". Difference:" << average(pfld-oldPfld)
                    << ". Tolerance:" << tol << endl;

                if (exitIfBad)
                {
                    FatalError<< exit(FatalError);
                }
            }

            // Restore patch field values and states
            static_cast<Field<Type>&>(pfld) = std::move(oldPfld);
            pfld.setUpdated(oldUpdated[patchi]);
            pfld.setManipulated(oldManipulated[patchi]);
        }
    }

    if (debug&2)
    {
        PoutInFunction
            << " Result of checking for field "
            << bfld[0].internalField().name() << " : " << allOk << endl;
    }

    return allOk;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::readField
(
    const Internal& iField,
    const dictionary& dict
)
{
    // DebugInFunction << nl;

    // Clear the boundary field if already initialised
    this->resize_null(bmesh_.size());

    label nUnset = this->size();

    // 1. Handle explicit patch names. Note that there can be only one explicit
    //    patch name since is key of dictionary.

    for (const entry& dEntry : dict)
    {
        const auto* subdict = dEntry.dictPtr();

        if (subdict && dEntry.keyword().isLiteral())
        {
            const label patchi = bmesh_.findPatchID(dEntry.keyword());

            if (patchi != -1)
            {
                this->set
                (
                    patchi,
                    PatchField<Type>::New
                    (
                        bmesh_[patchi],
                        iField,
                        *subdict
                    )
                );
                --nUnset;
            }
        }
    }

    if (nUnset == 0)
    {
        return;
    }


    // 2. Patch-groups. (using non-wild card entries of dictionaries)
    // (patchnames already matched above)
    // Note: in reverse order of entries in the dictionary (last
    // patchGroups wins). This is so it is consistent with dictionary wildcard
    // behaviour
    for (auto iter = dict.crbegin(); iter != dict.crend(); ++iter)
    {
        const entry& dEntry = *iter;
        const auto* subdict = dEntry.dictPtr();

        if (subdict && dEntry.keyword().isLiteral())
        {
            const labelList patchIds
            (
                bmesh_.indices(dEntry.keyword(), true)  // use patchGroups
            );

            for (const label patchi : patchIds)
            {
                if (!this->set(patchi))
                {
                    this->set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            bmesh_[patchi],
                            iField,
                            *subdict
                        )
                    );
                }
            }
        }
    }


    // 3. Wildcard patch overrides
    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            if (bmesh_[patchi].type() == emptyPolyPatch::typeName)
            {
                this->set
                (
                    patchi,
                    PatchField<Type>::New
                    (
                        emptyPolyPatch::typeName,
                        bmesh_[patchi],
                        iField
                    )
                );
            }
            else
            {
                const auto* subdict = dict.findDict(bmesh_[patchi].name());

                if (subdict)
                {
                    this->set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            bmesh_[patchi],
                            iField,
                            *subdict
                        )
                    );
                }
            }
        }
    }


    // Check for any unset patches
    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            if (bmesh_[patchi].type() == cyclicPolyPatch::typeName)
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find patchField entry for cyclic "
                    << bmesh_[patchi].name() << endl
                    << "Is your field uptodate with split cyclics?" << endl
                    << "Run foamUpgradeCyclics to convert mesh and fields"
                    << " to split cyclics." << exit(FatalIOError);
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find patchField entry for "
                    << bmesh_[patchi].name() << exit(FatalIOError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    // DebugInFunction << nl;

    forAll(bmesh_, patchi)
    {
        this->set
        (
            patchi,
            PatchField<Type>::New
            (
                patchFieldType,
                bmesh_[patchi],
                iField
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const wordList& patchFieldTypes,
    const wordList& constraintTypes
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    // DebugInFunction << nl;

    if
    (
        patchFieldTypes.size() != this->size()
     || (constraintTypes.size() && (constraintTypes.size() != this->size()))
    )
    {
        FatalErrorInFunction
            << "Incorrect number of patch type specifications given" << nl
            << "    Number of patches = " << bmesh.size()
            << " number of patchTypes = " << patchFieldTypes.size()
            << " number of contraintTypes = " << constraintTypes.size()
            << abort(FatalError);
    }

    forAll(bmesh_, patchi)
    {
        const word& constraintType =
        (
            patchi < constraintTypes.size()
          ? constraintTypes[patchi]
          : word::null
        );

        this->set
        (
            patchi,
            PatchField<Type>::New
            (
                patchFieldTypes[patchi],
                constraintType,
                bmesh_[patchi],
                iField
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const UPtrList<PatchField<Type>>& pflds
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    if (pflds.empty())
    {
        // Accept an empty PtrList - ie, populate later
    }
    else
    {
        forAll(bmesh_, patchi)
        {
            this->set(patchi, pflds[patchi].clone(iField));
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const Internal& iField,
    const GeometricBoundaryField& btf
)
:
    FieldField<PatchField, Type>(btf, iField),
    bmesh_(btf.bmesh_)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const Internal& iField,
    const GeometricBoundaryField& btf,
    const labelUList& patchIDs,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    // DebugInFunction << nl;

    for (const label patchi : patchIDs)
    {
        this->set
        (
            patchi,
            PatchField<Type>::New
            (
                patchFieldType,
                bmesh_[patchi],
                iField
            )
        );
    }

    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            this->set(patchi, btf[patchi].clone(iField));
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const dictionary& dict
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    readField(iField, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::updateCoeffs()
{
    // DebugInFunction << nl;

    for (auto& pfld : *this)
    {
        pfld.updateCoeffs();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluate
(
    const UPstream::commsTypes commsType
)
{
    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const auto startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            pfld.initEvaluate(commsType);
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : *this)
        {
            pfld.evaluate(commsType);
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            if (schedEval.init)
            {
                pfld.initEvaluate(commsType);
            }
            else
            {
                pfld.evaluate(commsType);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType) << nl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class UnaryPredicate>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluate_if
(
    const UnaryPredicate& pred,
    const UPstream::commsTypes commsType
)
{
    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        label nEvaluated(0);
        const auto startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            if (pred(pfld))
            {
                pfld.initEvaluate(commsType);
                ++nEvaluated;
            }
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        if (!nEvaluated) return;  // Early termination

        for (auto& pfld : *this)
        {
            if (pred(pfld))
            {
                pfld.evaluate(commsType);
                if (--nEvaluated == 0) break;  // Early termination
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            if (pred(pfld))
            {
                if (schedEval.init)
                {
                    pfld.initEvaluate(commsType);
                }
                else
                {
                    pfld.evaluate(commsType);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType) << nl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluateSelected
(
    const labelUList& patchIDs
)
{
    const auto startOfRequests = UPstream::nRequests();

    for (const label patchi : patchIDs)
    {
        auto& pfld = (*this)[patchi];

        DebugInfo<< "Updating " << pfld.patch().name() << endl;

        pfld.initEvaluate(UPstream::commsTypes::nonBlocking);
    }

    // Wait for outstanding requests (non-blocking)
    UPstream::waitRequests(startOfRequests);

    for (const label patchi : patchIDs)
    {
        auto& pfld = (*this)[patchi];

        pfld.evaluate(UPstream::commsTypes::nonBlocking);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluateLocal
(
    const UPstream::commsTypes commsType
)
{
    // DebugInFunction << nl;

    if (!FieldBase::localBoundaryConsistency())
    {
        return;
    }

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const auto startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            pfld.initEvaluateLocal(commsType);
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : *this)
        {
            pfld.evaluateLocal(commsType);
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            if (schedEval.init)
            {
                pfld.initEvaluateLocal(commsType);
            }
            else
            {
                pfld.evaluateLocal(commsType);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType) << nl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class CoupledPatchType>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluateCoupled
(
    const UPstream::commsTypes commsType
)
{
    if constexpr (std::is_void_v<CoupledPatchType>)
    {
        this->evaluate_if
        (
            [](const auto& pfld) { return pfld.coupled(); },
            commsType
        );
    }
    else
    {
        this->evaluate_if
        (
            [](const auto& pfld) -> bool
            {
                const auto* cpp = isA<CoupledPatchType>(pfld.patch());
                return (cpp && cpp->coupled());
            },
            commsType
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::wordList
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::types() const
{
    const FieldField<PatchField, Type>& pff = *this;

    wordList list(pff.size());

    forAll(pff, patchi)
    {
        if (const auto* pfld = pff.get(patchi))
        {
            list[patchi] = pfld->type();
        }
        else
        {
            list[patchi] = "null";
        }
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::
boundaryInternalField() const
{
    auto tresult = tmp<GeometricBoundaryField<Type, PatchField, GeoMesh>>::New
    (
        DimensionedField<Type, GeoMesh>::null(),
        *this
    );

    auto& result = tresult.ref();

    forAll(result, patchi)
    {
        result[patchi] == this->operator[](patchi).patchInternalField();
    }

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::LduInterfaceFieldPtrsList<Type>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::interfaces() const
{
    LduInterfaceFieldPtrsList<Type> list(this->size());

    typedef LduInterfaceField<Type> InterfaceType;

    forAll(list, patchi)
    {
        list.set
        (
            patchi,
            dynamic_cast<const InterfaceType*>(this->get(patchi))
        );
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::lduInterfaceFieldPtrsList
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::
scalarInterfaces() const
{
    lduInterfaceFieldPtrsList list(this->size());

    typedef lduInterfaceField InterfaceType;

    forAll(list, patchi)
    {
        list.set
        (
            patchi,
            dynamic_cast<const InterfaceType*>(this->get(patchi))
        );
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);
    this->writeEntries(os);
    os.endBlock();

    os.check(FUNCTION_NAME);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::writeEntries
(
    Ostream& os
) const
{
    for (const auto& pfld : *this)
    {
        os.beginBlock(pfld.patch().name());
        os << pfld;
        os.endBlock();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::check
(
) const
{
    // Dummy op - template specialisations provide logic (usually call
    // to checkConsistency)
    return true;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const FieldField<PatchField, Type>& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const Type& val
)
{
    FieldField<PatchField, Type>::operator=(val);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    Foam::zero
)
{
    FieldField<PatchField, Type>::operator=(Foam::zero{});
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const FieldField<PatchField, Type>& bf
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const Type& val
)
{
    for (auto& pfld : *this)
    {
        pfld == val;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    os << static_cast<const FieldField<PatchField, Type>&>(bf);
    return os;
}


// ************************************************************************* //
