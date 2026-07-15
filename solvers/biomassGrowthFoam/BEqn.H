// --- numerical controls (fvSolution/biofilm) ---
const dictionary& biofilmDict = mesh.solutionDict().subDict("biofilm");

const label  maxPicard  = biofilmDict.getOrDefault<label>("maxPicard", 20);
const scalar picardTol  = biofilmDict.getOrDefault<scalar>("picardTol", 1e-6);
const scalar capFrac    = biofilmDict.getOrDefault<scalar>("capFrac", 1e-2);
const scalar DBmaxRatio = biofilmDict.getOrDefault<scalar>("DBmaxRatio", 1e8);
const scalar DBrelax    = biofilmDict.getOrDefault<scalar>("DBrelax", 1.0);

// --- derived, not read ---
const scalar Bceil = 1.0 - capFrac;              // B is a volume fraction: Bmax = 1
const dimensionedScalar DBmax("DBmax", DBmaxRatio*delta);
const dimensionedScalar Bzero("Bzero", B.dimensions(), 0.0);

// --- loop state ---
scalar picardResid = GREAT;
label  picardIter  = 0;

while (picardResid > picardTol && picardIter < maxPicard)
{
    B.storePrevIter();                     // snapshot for the convergence measure

    // capped biomass — used ONLY to evaluate the coefficient, never assigned to B
    volScalarField Bc(min(max(B, Bzero), Bceil));

    DB = delta*pow(Bc, a)/pow(1.0 - Bc, b);
    DB = min(DB, DBmax);                    // stiffness guard
    DB.correctBoundaryConditions();

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix BEqn
        (
            fvm::ddt(B)
          - fvm::laplacian(DB, B)
          ==
            fvm::SuSp(mu*C/(Ks + C), B)      // positive coeff → treated explicitly
        );

        BEqn.relax();
        BEqn.solve();
    }

    // boundedness
    B.max(0.0);
    B.min(Bceil);
    B.correctBoundaryConditions();

    // Picard residual: max normalized change over this sweep
    picardResid = gMax(mag(B.primitiveField() - B.prevIter().primitiveField()))
                / (gMax(mag(B.primitiveField())) + SMALL);

    picardIter++;
}

Info<< "D(B) Picard: " << picardIter << " sweeps, final residual "
    << picardResid << endl;
