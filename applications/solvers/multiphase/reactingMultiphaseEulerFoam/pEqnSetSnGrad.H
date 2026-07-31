// Update the updateableSnGrad p_rgh BCs to ensure flux consistency
{
    surfaceScalarField::Boundary phib
    (
        phi.internalField(),
        phi.boundaryField()
    );
    phib = 0;
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        phib +=
            alphafs[phasei].boundaryField()*phase.phi()().boundaryField();
    }

    setSnGrad
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField() - phib
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );
}
