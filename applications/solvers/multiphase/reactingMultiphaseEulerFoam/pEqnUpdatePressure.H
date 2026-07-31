// Update and limit the static pressure
p = max(p_rgh + rho*gh, pMin);

// Limit p_rgh
p_rgh = p - rho*gh;

// Update densities from change in p_rgh
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    phase.thermoRef().rho() += phase.thermo().psi()*(p_rgh - p_rgh_0);
}

// Correct p_rgh for consistency with p and the updated densities
rho = fluid.rho();
p_rgh = p - rho*gh;
p_rgh.correctBoundaryConditions();
