// Construct the compressibility parts of the pressure equation
tmp<fvScalarMatrix> pEqnComp1, pEqnComp2;
if (phase1.compressible())
{
    if (pimple.transonic())
    {
        const surfaceScalarField phid1
        (
            IOobject::groupName("phid", phase1.name()),
            fvc::interpolate(psi1)*phi1
        );

        pEqnComp1 =
            (
                fvc::ddt(alpha1, rho1) + fvc::div(phase1.alphaRhoPhi())
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + correction
            (
                (alpha1/rho1)*
                (
                    psi1*fvm::ddt(p_rgh)
                  + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                )
            );

        pEqnComp1.ref().relax();
    }
    else
    {
        pEqnComp1 =
            (
                fvc::ddt(alpha1, rho1) + fvc::div(phase1.alphaRhoPhi())
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));
    }
}
if (phase2.compressible())
{
    if (pimple.transonic())
    {
        const surfaceScalarField phid2
        (
            IOobject::groupName("phid", phase2.name()),
            fvc::interpolate(psi2)*phi2
        );

        pEqnComp2 =
            (
                fvc::ddt(alpha2, rho2) + fvc::div(phase2.alphaRhoPhi())
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
            )/rho2
          + correction
            (
                (alpha2/rho2)*
                (
                    psi2*fvm::ddt(p_rgh)
                  + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
                )
            );

        pEqnComp2.ref().relax();
    }
    else
    {
        pEqnComp2 =
            (
                fvc::ddt(alpha2, rho2) + fvc::div(phase2.alphaRhoPhi())
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
            )/rho2
          + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
    }
}

// Add option sources
{
    if (fvOptions.appliesToField(rho1.name()))
    {
        tmp<fvScalarMatrix> optEqn1 = fvOptions(alpha1, rho1);
        if (pEqnComp1.valid())
        {
            pEqnComp1.ref() -= (optEqn1 & rho1)/rho1;
        }
        else
        {
            pEqnComp1 = fvm::Su(- (optEqn1 & rho1)/rho1, p_rgh);
        }
    }
    if (fvOptions.appliesToField(rho2.name()))
    {
        tmp<fvScalarMatrix> optEqn2 = fvOptions(alpha2, rho2);
        if (pEqnComp2.valid())
        {
            pEqnComp2.ref() -= (optEqn2 & rho2)/rho2;
        }
        else
        {
            pEqnComp2 = fvm::Su(- (optEqn2 & rho2)/rho2, p_rgh);
        }
    }
}

// Add mass transfer
{
    PtrList<volScalarField> dmdts(fluid.dmdts());
    if (dmdts.set(0))
    {
        if (pEqnComp1.valid())
        {
            pEqnComp1.ref() -= dmdts[0]/rho1;
        }
        else
        {
            pEqnComp1 = fvm::Su(- dmdts[0]/rho1, p_rgh);
        }
    }
    if (dmdts.set(1))
    {
        if (pEqnComp2.valid())
        {
            pEqnComp2.ref() -= dmdts[1]/rho2;
        }
        else
        {
            pEqnComp2 = fvm::Su(- dmdts[1]/rho2, p_rgh);
        }
    }
}
