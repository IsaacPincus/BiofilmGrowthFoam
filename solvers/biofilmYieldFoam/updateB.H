// BvolScalarField& B =
//     const_cast<volScalarField&>
//     (
//         mesh.lookupObject<volScalarField>("B")
//     );

// const surfaceScalarField& alphaPhi10 =
//     mesh.lookupObject<surfaceScalarField>
//     (
//         IOobject::groupName("alphaPhi0", alpha1.group())
//     );

// Advect B field with phase-1 flux
{
    fvScalarMatrix BEqn
    (
        fvm::ddt(B)
    //   + fvm::div(alphaPhi10, B)
        + fvm::div(phi,B)
    );
    BEqn.solve();

    // // Ensure physical bounds
    // B = max(B, dimensionedScalar("zero", dimViscosity/dimTime, 0.0));
}

// // Advect Bdead field with phase-1 flux
// {
//     fvScalarMatrix BdeadEqn
//     (
//         fvm::ddt(Bdead)
//     //   + fvm::div(alphaPhi10, B)
//         + fvm::div(phi,Bdead)
//     );
//     BdeadEqn.solve();

//     // // Ensure physical bounds
//     // B = max(B, dimensionedScalar("zero", dimViscosity/dimTime, 0.0));
// }


