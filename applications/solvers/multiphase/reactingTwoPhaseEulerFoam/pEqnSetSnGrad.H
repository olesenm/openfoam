// Update the updateableSnGrad p_rgh BCs to ensure flux consistency
setSnGrad
(
    p_rgh.boundaryFieldRef(),
    (
        phiHbyA.boundaryField()
      - (
            alphaf1.boundaryField()*phi1.boundaryField()
          + alphaf2.boundaryField()*phi2.boundaryField()
        )
    )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
);
