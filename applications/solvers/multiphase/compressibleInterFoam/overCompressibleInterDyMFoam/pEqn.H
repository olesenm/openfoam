{
    volScalarField rAU("rAU", 1.0/UEqn.A());

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
    );

    MRF.makeRelative(phiHbyA);

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*faceMask*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    #include "pEqnCompressibility.H"

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            p_rghEqnComp1() + p_rghEqnComp2() + p_rghEqnIncomp,
            p_rgh.select(pimple.finalInnerIter())
        );

        if (pimple.finalNonOrthogonalIter())
        {
            p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);
            p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh;

            dgdt =
            (
                alpha1*(p_rghEqnComp2 & p_rgh)
              - alpha2*(p_rghEqnComp1 & p_rgh)
            );

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U =
                cellMask*
                (
                    HbyA
                  + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf)
                );

            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(fvc::absolute(phi, U)/mesh.magSf() - (n & Uf));
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    // Zero faces H-I for transport Eq after pEq
    phi *= faceMask;

    #include "pEqnDensityUpdate.H"
}
