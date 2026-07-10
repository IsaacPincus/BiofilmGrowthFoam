// --- Picard / capping controls (park these in fvSolution's SIMPLE subdict) ---
const label   maxPicard = 20;
const scalar  picardTol = 1e-6;
const scalar  capFrac   = 1e-2;   // stay this far below Bmax when evaluating D
const dimensionedScalar Bzero("Bzero", B.dimensions(), 0.0);
const dimensionedScalar Bceil = (1.0 - capFrac)*Bmax;

// optional hard ceiling on diffusivity, decoupled from capFrac (see note)
const dimensionedScalar DBmax("DBmax", DB.dimensions(), 1e8*delta.value());

scalar picardResid = GREAT;
label  picardIter  = 0;

while (picardResid > picardTol && picardIter < maxPicard)
{
    B.storePrevIter();                     // snapshot for the convergence measure

    // capped biomass — used ONLY to evaluate the coefficient, never assigned to B
    volScalarField Bc(min(max(B, Bzero), Bceil));
    volScalarField M(Bc/Bmax);

    DB = delta*pow(M, a)/pow(1.0 - M, b);
    DB = min(DB, DBmax);                    // stiffness guard
    DB.correctBoundaryConditions();

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix BEqn
        (
            fvm::ddt(B)
          - fvm::laplacian(DB, B)
          ==
            fvm::SuSp(mu*C/(k + C), B)      // positive coeff → treated explicitly
        );

        BEqn.relax();
        BEqn.solve();
    }

    // boundedness: clamp is a safety net, not a routine correction (see note)
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