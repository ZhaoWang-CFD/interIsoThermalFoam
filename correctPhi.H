CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    energy.vDot(),
    pimple
);

#include "continuityErrs.H"
