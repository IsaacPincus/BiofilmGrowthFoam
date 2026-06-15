// Before the time loop, declare storage for old alpha
volScalarField alpha1Old
(
    IOobject
    (
        "alpha1Old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1
);

label checkInterval = 100;  // compare every N steps
scalar stopTolerance = 1e-6;

