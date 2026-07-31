tmp<fvScalarMatrix> pEqnComp1;
tmp<fvScalarMatrix> pEqnComp2;

// Construct the compressibility parts of the pressure equation
if (pimple.transonic())
{
    surfaceScalarField phid1
    (
        IOobject::groupName("phid", phase1.name()),
        fvc::interpolate(psi1)*phi1
    );
    surfaceScalarField phid2
    (
        IOobject::groupName("phid", phase2.name()),
        fvc::interpolate(psi2)*phi2
    );

    pEqnComp1 =
        (
            contErr1
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

    pEqnComp2 =
        (
            contErr2
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
    pEqnComp1 =
        (
            contErr1
          - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
        )/rho1
      + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));

    pEqnComp2 =
        (
            contErr2
          - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
        )/rho2
      + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
}
