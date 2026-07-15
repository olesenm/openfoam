// Compressible pressure equations
PtrList<fvScalarMatrix> pEqnComps(phases.size());
PtrList<volScalarField> dmdts(fluid.dmdts());
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    const volScalarField& alpha = phase;
    volScalarField& rho = phase.thermoRef().rho();

    if (phase.compressible())
    {
        if (pimple.transonic())
        {
            surfaceScalarField phid
            (
                IOobject::groupName("phid", phase.name()),
                fvc::interpolate(phase.thermo().psi())*phase.phi()
            );

            pEqnComps.set
            (
                phasei,
                (
                    (
                        fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                      - fvc::Sp
                        (
                            fvc::ddt(alpha) + fvc::div(phase.alphaPhi()),
                            rho
                        )
                    )/rho
                  + correction
                    (
                        (alpha/rho)*
                        (
                            phase.thermo().psi()*fvm::ddt(p_rgh)
                          + fvm::div(phid, p_rgh)
                          - fvm::Sp(fvc::div(phid), p_rgh)
                        )
                    )
                ).ptr()
            );

            if (clearFaceFluxCorrection)
            {
                pEqnComps[phasei].faceFluxCorrectionPtr(nullptr);
            }

            pEqnComps[phasei].relax();
        }
        else
        {
            pEqnComps.set
            (
                phasei,
                (
                    (
                        fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                      - fvc::Sp
                        (
                            (fvc::ddt(alpha) + fvc::div(phase.alphaPhi())),
                            rho
                        )
                    )/rho
                  + (alpha*phase.thermo().psi()/rho)
                   *correction(fvm::ddt(p_rgh))
                ).ptr()
            );
        }
    }

    // Add option sources
    if (fvOptions.appliesToField(rho.name()))
    {
        tmp<fvScalarMatrix> optEqn = fvOptions(alpha, rho);
        if (pEqnComps.set(phasei))
        {
            pEqnComps[phasei] -= (optEqn&rho)/rho;
        }
        else
        {
            pEqnComps.set
            (
                phasei,
                fvm::Su(- (optEqn&rho)/rho, p_rgh).ptr()
            );
        }
    }

    // Add mass transfer
    if (dmdts.set(phasei))
    {
        if (pEqnComps.set(phasei))
        {
            pEqnComps[phasei] -= dmdts[phasei]/rho;
        }
        else
        {
            pEqnComps.set
            (
                phasei,
                fvm::Su(- dmdts[phasei]/rho, p_rgh)
            );
        }
    }
}
