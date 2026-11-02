/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "lduMatrix.H"
#include "SubField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduMatrix::initMatrixInterfaces
(
    const bool add,
    const scalarField& interfaceCoeffsData,
    const labelUList& offsets,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt
) const
{
    const UPstream::commsTypes commsType = UPstream::defaultCommsType;

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                const label start = offsets[interfacei];
                const label len = offsets[interfacei+1] - offsets[interfacei];
                const SubField<scalar> coeffs(interfaceCoeffsData, len, start);

                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coeffs,
                    cmpt,
                    commsType
                );
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                const label start = offsets[interfacei];
                const label len = offsets[interfacei+1] - offsets[interfacei];
                const SubField<scalar> coeffs(interfaceCoeffsData, len, start);

                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coeffs,
                    cmpt,
                    UPstream::commsTypes::buffered
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << exit(FatalError);
    }
}


void Foam::lduMatrix::updateMatrixInterfaces
(
    const bool add,
    const scalarField& interfaceCoeffsData,
    const labelUList& offsets,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt,
    const int startRequest
) const
{
    const UPstream::commsTypes commsType = UPstream::defaultCommsType;

    if
    (
        commsType == UPstream::commsTypes::nonBlocking
     && UPstream::nPollProcInterfaces
    )
    {
        // Poll for completion and consume ready interfaces

        DynamicList<int> indices;  // (work array)

        const int maxPolling =
        (
            (UPstream::nPollProcInterfaces < 0)
          ? (UPstream::nRequests() - startRequest)
          : UPstream::nPollProcInterfaces
        );

        for (int cycle = 0; cycle < maxPolling; ++cycle)
        {
            bool allUpdated = true;

            forAll(interfaces, interfacei)
            {
                auto* intf = interfaces.get(interfacei);

                if (intf && !intf->updatedMatrix())
                {
                    if (intf->ready())
                    {
                        const label start = offsets[interfacei];
                        const label len = offsets[interfacei+1] - offsets[interfacei];
                        const SubField<scalar> coeffs(interfaceCoeffsData, len, start);

                        intf->updateInterfaceMatrix
                        (
                            result,
                            add,
                            mesh().lduAddr(),
                            interfacei,
                            psiif,
                            coeffs,
                            cmpt,
                            commsType
                        );
                    }
                    else
                    {
                        allUpdated = false;
                    }
                }
            }

            if (allUpdated)
            {
                break;  // Early exit
            }
        }
    }


    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        // Wait until sends/receives have finished.
        // - effectively a no-op (without waiting) if already completed.
        if (commsType == UPstream::commsTypes::nonBlocking)
        {
            UPstream::waitRequests(startRequest);
        }

        // Check/no-check for updatedMatrix() ?
        const bool noCheck = (commsType == UPstream::commsTypes::buffered);

        // Consume anything still outstanding
        forAll(interfaces, interfacei)
        {
            auto* intf = interfaces.get(interfacei);

            if (intf && (noCheck || !intf->updatedMatrix()))
            {
                const label start = offsets[interfacei];
                const label len = offsets[interfacei+1] - offsets[interfacei];
                const SubField<scalar> coeffs(interfaceCoeffsData, len, start);

                intf->updateInterfaceMatrix
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coeffs,
                    cmpt,
                    commsType
                );
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over patches according to the schedule
        forAll(patchSchedule, patchEvali)
        {
            const label interfacei = patchSchedule[patchEvali].patch;

            if (interfaces.set(interfacei))
            {
                if (patchSchedule[patchEvali].init)
                {
                    const label start = offsets[interfacei];
                    const label len = offsets[interfacei+1] - offsets[interfacei];
                    const SubField<scalar> coeffs(interfaceCoeffsData, len, start);

                    interfaces[interfacei].initInterfaceMatrixUpdate
                    (
                        result,
                        add,
                        mesh().lduAddr(),
                        interfacei,
                        psiif,
                        coeffs,
                        cmpt,
                        UPstream::commsTypes::buffered
                    );
                }
                else
                {
                    const label start = offsets[interfacei];
                    const label len = offsets[interfacei+1] - offsets[interfacei];
                    const SubField<scalar> coeffs(interfaceCoeffsData, len, start);

                    interfaces[interfacei].updateInterfaceMatrix
                    (
                        result,
                        add,
                        mesh().lduAddr(),
                        interfacei,
                        psiif,
                        coeffs,
                        cmpt,
                        UPstream::commsTypes::buffered
                    );
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                const label start = offsets[interfacei];
                const label len = offsets[interfacei+1] - offsets[interfacei];
                const SubField<scalar> coeffs(interfaceCoeffsData, len, start);

                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coeffs,
                    cmpt,
                    UPstream::commsTypes::buffered
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
