// Update densities from change in p_rgh
mixture.thermo1().correctRho(psi1*(p_rgh - p_rgh_0));
mixture.thermo2().correctRho(psi2*(p_rgh - p_rgh_0));

rho = alpha1*rho1 + alpha2*rho2;

// Correct p_rgh for consistency with p and the updated densities
p_rgh = p - rho*gh;
p_rgh.correctBoundaryConditions();

K = 0.5*magSqr(U);
