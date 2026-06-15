// perhaps in the future this will just be calculated directly from B and Bdead

volScalarField& tau0 =
    const_cast<volScalarField&>
    (
        mesh.lookupObject<volScalarField>("tau0")
    );

// const surfaceScalarField& alphaPhi10 =
//     mesh.lookupObject<surfaceScalarField>
//     (
//         IOobject::groupName("alphaPhi0", alpha1.group())
//     );

// Advect tau0 field with phase-1 flux
{
    fvScalarMatrix tau0Eqn
    (
        fvm::ddt(tau0)
    //   + fvm::div(alphaPhi10, tau0)
        + fvm::div(phi,tau0)
    );
    tau0Eqn.solve();

    // Ensure physical bounds
    tau0 = max(tau0, dimensionedScalar("zero", dimViscosity/dimTime, 0.0));
}

