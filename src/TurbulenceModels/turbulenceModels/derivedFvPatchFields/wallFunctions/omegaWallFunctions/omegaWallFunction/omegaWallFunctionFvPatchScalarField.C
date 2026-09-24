/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
    Copyright (C) 2026 UKRI-STFC.
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

    Some parts here added have been inspired by the work done at the following
    repository https://github.com/SimFlowCFD/RapidCFD-dev

\*---------------------------------------------------------------------------*/

#include "omegaWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::omegaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::omegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<this_bctype>(bf[patchi]))
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void Foam::omegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<this_bctype>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            #ifdef FOAM_OFFLOAD
            const auto& pcells = mesh.lduAddr().patchSortCells(patchi);
            const auto& losortStart = mesh.lduAddr().patchSortStartAddr(patchi);

            labelRange cellRange(0, pcells.size());

            apply_policy(cellRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    [&](label i)
                    {
                        label delta = losortStart[i + 1] - losortStart[i];
                        weights[pcells[i]] += static_cast<scalar>(delta);
                    });
            });
            #else
            for (const auto& celli : faceCells)
            {
                ++weights[celli];
            }
            #endif
        }
    }

    cornerWeights_.setSize(bf.size());
    for (const auto& patchi : omegaPatches)
    {
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


Foam::omegaWallFunctionFvPatchScalarField&
Foam::omegaWallFunctionFvPatchScalarField::omegaPatch
(
    const label patchi
)
{
    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    return refConstCast<this_bctype>(omega.boundaryField()[patchi]);
}


void Foam::omegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void Foam::omegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    const scalarField& y = turbModel.y()[patchi];

    const labelUList& faceCells = patch.faceCells();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    // Calculate y-plus
    const auto yPlus = [&](const label facei) -> scalar
    {
        return
        (
            Cmu25*y[facei]*sqrt(k[faceCells[facei]])/nuw[facei]
        );
    };

    // Contribution from the viscous sublayer
    const auto omegaVis = [&](const label facei) -> scalar
    {
        return
        (
            6.0*nuw[facei]/(beta1_*sqr(y[facei]))
        );
    };

    // Contribution from the inertial sublayer
    const auto omegaLog = [&](const label facei) -> scalar
    {
        return
        (
            sqrt(k[faceCells[facei]])
          / (Cmu25*kappa*y[facei])
        );
    };

    #ifdef FOAM_OFFLOAD
    const auto& patchCells =
        patch.boundaryMesh().mesh().lduAddr().patchSortCells(patchi);
    const auto& sortStart =
        patch.boundaryMesh().mesh().lduAddr().patchSortStartAddr(patchi);
    const auto& sort =
        patch.boundaryMesh().mesh().lduAddr().patchSortAddr(patchi);

    const labelRange cellRange(patchCells.size());

    auto neiOp = [](scalar& result, const scalar value)
    {
        result += value;
    };
    #endif

    switch (blender_)
    {
        case blenderType::STEPWISE:
        {
            #ifdef FOAM_OFFLOAD
            auto neiFun = [&](const label, const label facei) -> scalar
            {
                if (yPlus(facei) > yPlusLam)
                {
                    return cornerWeights[facei]*omegaLog(facei);
                }

                return cornerWeights[facei]*omegaVis(facei);
            };

            apply_policy(cellRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    patchCellReduceFunctor
                    <
                        scalar,
                        decltype(neiFun),
                        decltype(neiOp)
                    >
                    (
                        omega0,
                        sortStart,
                        sort,
                        patchCells,
                        neiFun,
                        neiOp
                    )
                );
            });
            #else
            forAll(faceCells, facei)
            {
                if (yPlus(facei) > yPlusLam)
                {
                    omega0[faceCells[facei]] +=
                        cornerWeights[facei]
                      * omegaLog(facei);
                }
                else
                {
                    omega0[faceCells[facei]] +=
                        cornerWeights[facei]
                      * omegaVis(facei);
                }
            }
            #endif
            break;
        }

        case blenderType::BINOMIAL:
        {
            #ifdef FOAM_OFFLOAD
            auto neiFun = [&](const label, const label facei) -> scalar
            {
                return cornerWeights[facei] * pow
                    (
                        pow(omegaVis(facei), n_) + pow(omegaLog(facei), n_),
                        scalar(1)/n_
                    );
            };

            apply_policy(cellRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    patchCellReduceFunctor
                    <
                        scalar,
                        decltype(neiFun),
                        decltype(neiOp)
                    >
                    (
                        omega0,
                        sortStart,
                        sort,
                        patchCells,
                        neiFun,
                        neiOp
                    )
                );
            });
            #else
            forAll(faceCells, facei)
            {
                omega0[faceCells[facei]] +=
                    cornerWeights[facei]
                  * pow
                    (
                        pow(omegaVis(facei), n_) + pow(omegaLog(facei), n_),
                        scalar(1)/n_
                    );
            }
            #endif
            break;
        }

        case blenderType::MAX:
        {
            #ifdef FOAM_OFFLOAD
            auto neiFun = [&](const label, const label facei) -> scalar
            {
                return cornerWeights[facei]
                    * max(omegaVis(facei), omegaLog(facei));
            };

            apply_policy(cellRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    patchCellReduceFunctor
                    <
                        scalar,
                        decltype(neiFun),
                        decltype(neiOp)
                    >
                    (
                        omega0,
                        sortStart,
                        sort,
                        patchCells,
                        neiFun,
                        neiOp
                    )
                );
            });
            #else
            forAll(faceCells, facei)
            {
                // (PH:Eq. 27)
                omega0[faceCells[facei]] +=
                    cornerWeights[facei]
                  * max(omegaVis(facei), omegaLog(facei));
            }
            #endif
            break;
        }

        case blenderType::EXPONENTIAL:
        {
            #ifdef FOAM_OFFLOAD
            auto neiFun = [&](const label, const label facei) -> scalar
            {
                const scalar yPlusFace = yPlus(facei);
                const scalar Gamma = 0.01*pow4(yPlusFace)/(1 + 5*yPlusFace);
                const scalar invGamma = scalar(1)/(Gamma + ROOTVSMALL);

                return cornerWeights[facei]
                    * (
                        omegaVis(facei)*exp(-Gamma)
                    + omegaLog(facei)*exp(-invGamma)
                    );
            };

            apply_policy(cellRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    patchCellReduceFunctor
                    <
                        scalar,
                        decltype(neiFun),
                        decltype(neiOp)
                    >
                    (
                        omega0,
                        sortStart,
                        sort,
                        patchCells,
                        neiFun,
                        neiOp
                    )
                );
            });
            #else
            forAll(faceCells, facei)
            {
                // (PH:Eq. 31)
                const scalar yPlusFace = yPlus(facei);
                const scalar Gamma = 0.01*pow4(yPlusFace)/(1 + 5*yPlusFace);
                const scalar invGamma = scalar(1)/(Gamma + ROOTVSMALL);

                omega0[faceCells[facei]] +=
                    cornerWeights[facei]
                  * (
                        omegaVis(facei)*exp(-Gamma)
                      + omegaLog(facei)*exp(-invGamma)
                    );
            }
            #endif
            break;
        }

        case blenderType::TANH:
        {
            #ifdef FOAM_OFFLOAD
            auto neiFun = [&](const label, const label facei) -> scalar
            {
                const scalar omegaVisFace = omegaVis(facei);
                const scalar omegaLogFace = omegaLog(facei);

                const scalar b1 = omegaVisFace + omegaLogFace;

                const scalar b2 =
                    pow
                    (
                        pow(omegaVisFace, 1.2) + pow(omegaLogFace, 1.2),
                        1.0/1.2
                    );

                const scalar phiTanh = tanh(pow4(0.1*yPlus(facei)));

                return cornerWeights[facei]
                *(phiTanh*b1 + (1 - phiTanh)*b2);
            };

            apply_policy(cellRange.size(), [&](auto policy)
            {
                std::for_each
                (
                    policy,
                    cellRange.begin(),
                    cellRange.end(),
                    patchCellReduceFunctor
                    <
                        scalar,
                        decltype(neiFun),
                        decltype(neiOp)
                    >
                    (
                        omega0,
                        sortStart,
                        sort,
                        patchCells,
                        neiFun,
                        neiOp
                    )
                );
            });
            #else
            forAll(faceCells, facei)
            {
                // (KAS:Eqs. 33-34)
                const scalar omegaVisFace = omegaVis(facei);
                const scalar omegaLogFace = omegaLog(facei);
                const scalar b1 = omegaVisFace + omegaLogFace;
                const scalar b2 =
                    pow
                    (
                        pow(omegaVisFace, 1.2) + pow(omegaLogFace, 1.2),
                        1.0/1.2
                    );
                const scalar phiTanh = tanh(pow4(0.1*yPlus(facei)));

                omega0[faceCells[facei]] +=
                    cornerWeights[facei]
                   *(phiTanh*b1 + (1 - phiTanh)*b2);
            }
            #endif
            break;
        }
    }

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));

    #ifdef FOAM_OFFLOAD
    auto neiFun = [&](const label, const label facei) -> scalar
    {
        const label celli = faceCells[facei];

        if (!(blender_ == blenderType::STEPWISE) || yPlus(facei) > yPlusLam)
        {
            return cornerWeights[facei]
            *(nutw[facei] + nuw[facei])
            *magGradUw[facei]
            *Cmu25*sqrt(k[celli])
            /(kappa*y[facei]);
        }

        return Zero;
    };

    apply_policy(cellRange.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            patchCellReduceFunctor
            <
                scalar,
                decltype(neiFun),
                decltype(neiOp)
            >
            (
                omega0,
                sortStart,
                sort,
                patchCells,
                neiFun,
                neiOp
            )
        );
    });
    #else
    forAll(faceCells, facei)
    {
        if (!(blender_ == blenderType::STEPWISE) || yPlus(facei) > yPlusLam)
        {
            G0[faceCells[facei]] +=
                cornerWeights[facei]
               *(nutw[facei] + nuw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(k[faceCells[facei]])
               /(kappa*y[facei]);
        }
    }
    #endif
}


void Foam::omegaWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
    os.writeEntryIfDifferent<scalar>("beta1", 0.075, beta1_);
    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    wallFunctionBlenders(),
    initialised_(false),
    master_(-1),
    beta1_(0.075),
    wallCoeffs_(),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf),
    initialised_(false),
    master_(-1),
    beta1_(ptf.beta1_),
    wallCoeffs_(ptf.wallCoeffs_),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::BINOMIAL, scalar(2)),
    initialised_(false),
    master_(-1),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.075)),
    wallCoeffs_(dict),
    G_(),
    omega_(),
    cornerWeights_()
{
    // Apply zero-gradient condition on start-up
    this->extrapolateInternal();
}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const this_bctype& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(owfpsf, iF),
    wallFunctionBlenders(owfpsf),
    initialised_(false),
    master_(-1),
    beta1_(owfpsf.beta1_),
    wallCoeffs_(owfpsf.wallCoeffs_),
    G_(),
    omega_(),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::omegaWallFunctionFvPatchScalarField::G
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


Foam::scalarField& Foam::omegaWallFunctionFvPatchScalarField::omega
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void Foam::omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& omega = const_cast<FieldType&>(internalField());

    const labelUList& faceCells = patch().faceCells();

    labelRange faceRange(faceCells.size());

    apply_policy(faceCells.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            faceRange.begin(),
            faceRange.end(),
            [&](const label facei)
            {
                const label celli = faceCells[facei];

                G[celli] = G0[celli];
                omega[celli] = omega0[celli];
            }
        );
    });

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::omegaWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());
    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    const labelUList& faceCells = patch().faceCells();

#ifndef FOAM_OFFLOAD

    forAll(weights, facei)
    {
        const scalar w = weights[facei];

        if (w > tolerance_)
        {
            const label celli = faceCells[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

#else

    // only set the values if the weights are > tolerance

    const scalar tolerance = tolerance_;

    const auto& lduAddr = patch().boundaryMesh().mesh().lduAddr();
    const auto& patchCells = lduAddr.patchSortCells(patch().index());
    const auto& neiStart = lduAddr.patchSortStartAddr(patch().index());
    const auto& losort = lduAddr.patchSortAddr(patch().index());
    const label* const __restrict__ patchCellsPtr = patchCells.begin();
    const label* const __restrict__ neiStartPtr = neiStart.begin();
    const label* const __restrict__ losortPtr = losort.begin();

    const labelRange patchCellRange(weights.size());

    const scalar* const __restrict__ weightsPtr = weights.begin();
    const label* const __restrict__ faceCellsPtr = faceCells.begin();

    scalar* const __restrict__ GPtr = G.begin();
    scalar* const __restrict__ omegaPtr = omega.begin();
    scalar* const __restrict__ omegafPtr = omegaf.begin();

    const scalar* const __restrict__ G0Ptr = G0.begin();
    const scalar* const __restrict__ omega0Ptr = omega0.begin();

    apply_policy(patchCellRange.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            patchCellRange.begin(),
            patchCellRange.end(),
            [=](const label id)
            {
                const label celli = patchCellsPtr[id];

                scalar Gval = GPtr[celli];
                scalar omegaval = omegaPtr[celli];

                const label start = neiStartPtr[id];
                const label end = neiStartPtr[id + 1];

                for (label i = start; i < end; ++i)
                {
                    const label facei = losortPtr[i];
                    const scalar w = weightsPtr[facei];

                    if (w > tolerance)
                    {
                        Gval =
                            (1.0 - w)*Gval + w*G0Ptr[celli];

                        omegaval =
                            (1.0 - w)*omegaval + w*omega0Ptr[celli];

                        omegafPtr[facei] = omegaval;
                    }
                }

                GPtr[celli] = Gval;
                omegaPtr[celli] = omegaval;
            }
        );
    });

#endif

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintValues(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& fld = internalField();

    #ifdef FOAM_OFFLOAD
    const labelRange faceRange(weights.size());

    List<label> selected(weights.size(), -1);

    const scalar tolerance = tolerance_;

    const scalar* const __restrict__ weightsPtr = weights.begin();
    const label* const __restrict__ faceCellsPtr = faceCells.begin();
    label* const __restrict__ selectedPtr = selected.begin();

    apply_policy(faceRange.size(), [&](auto policy)
    {
        std::transform
        (
            policy,
            faceRange.begin(),
            faceRange.end(),
            selectedPtr,
            [=](const label facei) -> label
            {
                if (tolerance < weightsPtr[facei])
                {
                    return faceCellsPtr[facei];
                }

                return label(-1);
            }
        );
    });

    // Compaction step. Keep this on CPU unless you have a proper parallel
    // compact.
    constraintCells.clear();
    constraintValues.clear();

    forAll(selected, i)
    {
        const label celli = selected[i];

        if (celli != -1)
        {
            constraintCells.append(celli);
            constraintValues.append(fld[celli]);
        }
    }
    #else
    forAll(weights, facei)
    {
        // only set the values if the weights are > tolerance
        if (tolerance_ < weights[facei])
        {
            const label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintValues.append(fld[celli]);
        }
    }
    #endif

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << constraintCells.size()
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues(constraintCells, constraintValues);

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::omegaWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        omegaWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
