Info<< "Reading field C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading local biofilm concentration B\n" << endl;
volScalarField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Info<< "Reading local dead biofilm concentration Bdead\n" << endl;
// volScalarField Bdead
// (
//     IOobject
//     (
//         "Bdead",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );





Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity 'D'\n" << endl;
dimensionedScalar D(
    "D", 
    dimViscosity, // in this case, this is m^2/s since it's a dynamic viscosity
    transportProperties
);

Info<< "Reading biofilm growth rate mu\n" << endl;
dimensionedScalar mu
(
    "mu",
    dimensionSet(0, 0, -1, 0, 0, 0, 0), // 1/s
    transportProperties
);

Info<< "Reading Monod half-saturation constant Ks\n" << endl;
dimensionedScalar Ks
(
    "Ks",
    dimensionSet(0, -3, 0, 0, 1, 0, 0), // mol / m^3
    transportProperties
);

Info<< "Reading biomass yield rate Y\n" << endl;
dimensionedScalar Y
(
    "Y",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), // dimensionless
    transportProperties
);

// this calculates the surface fluxes fvc::flux(U)
#include "createPhi.H"

#include "createFvOptions.H"

/* end */