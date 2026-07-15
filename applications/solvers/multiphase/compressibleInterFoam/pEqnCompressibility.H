tmp<fvScalarMatrix> p_rghEqnComp1;
tmp<fvScalarMatrix> p_rghEqnComp2;

if (pimple.transonic())
{
    #include "rhofs.H"

    surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
    surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

    p_rghEqnComp1 =
        pos(alpha1)
       *(
            (
                fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f)
              - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
            )/rho1
          - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
          + (alpha1/rho1)
           *correction
            (
                psi1*fvm::ddt(p_rgh)
              + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
            )
        );
    p_rghEqnComp1.ref().relax();

    p_rghEqnComp2 =
        pos(alpha2)
       *(
           (
               fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f)
             - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
           )/rho2
         - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
         + (alpha2/rho2)
          *correction
           (
               psi2*fvm::ddt(p_rgh)
             + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
           )
       );
    p_rghEqnComp2.ref().relax();
}
else
{
    #include "rhofs.H"

    p_rghEqnComp1 =
        pos(alpha1)
       *(
            (
                fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f)
              - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
            )/rho1
          - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
          + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh))
        );

    p_rghEqnComp2 =
        pos(alpha2)
       *(
           (
               fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f)
             - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
           )/rho2
         - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
         + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh))
        );
}
