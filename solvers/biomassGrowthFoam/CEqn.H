
while (simple.correctNonOrthogonal())
{
    C.max(0.); 
    fvScalarMatrix ConvectionEqn //(C=Concentration=rho*wf/mw)
    (
        fvm::ddt(C) 
        + fvm::div(phi,C) 
        - fvm::laplacian(D,C)
        + fvm::Sp(mu/Y*Binf*B/(Ks + C), C) // implicit treatment of source term
        // ==
        // fvOptions(C)
        // - mu*B*C/(Ks + C)/Y // biomass consumption, Monod kinetics
        //+ fvm::div(phic,C) //solid movement correction

        // // simplified version to test
        // fvm::ddt(C) 
        // + fvm::div(phi,C) 
        // - fvm::laplacian(Dion,C)
    );
    ConvectionEqn.relax();
    // fvOptions.constrain(ConvectionEqn);
    ConvectionEqn.solve();
    // fvOptions.correct(C);
    // Update boundary conditions
    C.correctBoundaryConditions();

    // Dion=Df*pow((1-epss),(n-1));
    // epsf=1-epss;
}
