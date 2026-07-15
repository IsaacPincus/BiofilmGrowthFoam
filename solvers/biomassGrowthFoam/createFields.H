Info<< "Reading field C\n" << endl;
volScalarField C(
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
volVectorField U(
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
volScalarField B(
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

// these values must be present, but aren't actually convected
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field tau0\n" << endl;
volScalarField tau0
(
    IOobject
    (
        "tau0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field alpha.water\n" << endl;
volScalarField alphaWater
(
    IOobject
    (
        "alpha.water",       // <-- the on-disk field name; keep the dot here
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// Info<< "Reading field p\n" << endl;
// volScalarField p
// (
//     IOobject
//     (
//         "p",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );


// biomass diffusivity — recomputed each Picard pass, so no I/O
Info<< "Creating biomass diffusivity field DB\n" << endl;
volScalarField DB
(
    IOobject
    (
        "DB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DB", dimViscosity, 0.0),
    zeroGradientFvPatchScalarField::typeName
);


Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading biofilm spreading prefactor delta\n" << endl;
dimensionedScalar delta
(
    "delta",
    dimViscosity,   // [m^2/s] — same dims as D, but a distinct coefficient
    transportProperties
);

Info<< "Reading D(B) exponents a, b\n" << endl;
scalar a(transportProperties.get<scalar>("a"));   // M^a  (degeneracy at B=0)
scalar b(transportProperties.get<scalar>("b"));   // (1-M)^-b (singularity at B=1)

Info<< "Reading diffusivity 'D'\n" << endl;
dimensionedScalar D(
    "D", 
    dimViscosity, // in this case, this is m^2/s since it's a dynamic viscosity
    transportProperties
);

Info<< "Reading biofilm growth rate mu\n" << endl;
dimensionedScalar mu(
    "mu",
    dimensionSet(0, 0, -1, 0, 0, 0, 0), // 1/s
    transportProperties
);

Info<< "Reading Monod half-saturation constant Ks\n" << endl;
dimensionedScalar Ks(
    "Ks",
    dimensionSet(0, -3, 0, 0, 1, 0, 0), // mol / m^3
    transportProperties
);

Info<< "Reading biomass yield Y, in kg biomass/mol of limiting substrate\n" << endl;
dimensionedScalar Y(
    "Y",
    dimensionSet(1, 0, 0, 0, -1, 0, 0),  // [kg mol^-1] — kg biomass per mol substrate
    transportProperties
);

Info<< "Reading maximum cell density Binf\n" << endl;
dimensionedScalar Binf(
    "Binf",
    dimDensity,                          // [kg m^-3] — biomass density
    transportProperties
);


// for A equation
Info<< "Reading autoinducer concentration A\n" << endl;
volScalarField A(
    IOobject("A", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Reading AI-2 diffusivity DA\n" << endl;
dimensionedScalar DA("DA", dimViscosity, transportProperties);   // [m^2/s]

Info<< "Reading tracer switch parameters\n" << endl;
dimensionedScalar Cstar(
    "Cstar",
    Ks.dimensions(),          // same units as C  [mol/m^3]
    transportProperties
);
dimensionedScalar wSwitch(
    "wSwitch",
    Ks.dimensions(),          // transition width, same units as C
    transportProperties
);

// switch field: recomputed from C every step, no I/O
volScalarField fSwitch(
    IOobject("fSwitch", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("fSwitch", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// B is a volume fraction (dimensionless), so k_A± must carry [A]/s themselves
Info<< "Reading AI-2 production / uptake rates\n" << endl;
dimensionedScalar kAplus ("kAplus",  A.dimensions()/dimTime, transportProperties);
dimensionedScalar kAminus("kAminus", A.dimensions()/dimTime, transportProperties);

// this calculates the surface fluxes fvc::flux(U)
#include "createPhi.H"

#include "createFvOptions.H"





/* end */
