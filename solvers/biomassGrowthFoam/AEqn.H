
fSwitch = 0.5*(1.0 - tanh((C - Cstar)/wSwitch));
fSwitch.correctBoundaryConditions();

while (simple.correctNonOrthogonal())
{
    A.max(0.); 
    fvScalarMatrix AEqn //(C=Concentration=rho*wf/mw)
    (
        fvm::ddt(A) 
        + fvm::div(phi,A) 
        - fvm::laplacian(DA,A)
        // creation term
        - kAplus*B
        // consumption term, can be either constant with kAplus > kAminus, or saturable
        + kAminus*fSwitch*B
        // - fvm::Sp(kAminus*fSwitch*B/(KA + A), A)     // saturable uptake, bounded at A=0

    );
    AEqn.relax();
    AEqn.solve();
    A.correctBoundaryConditions();
}
