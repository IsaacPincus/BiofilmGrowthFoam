/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    convectionWithBiomassFoam

Description
    Solves the steady or transient transport equation for a scalar, and can include a reaction rate dependent constant

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    // Info<< "C dimensions: " << C.dimensions() << nl;
    // Info<< "B dimensions: " << B.dimensions() << nl;
    // Info<< "Ks dimensions: " << Ks.dimensions() << nl;
    // Info<< "mu dimensions: " << mu.dimensions() << nl;
    // Info<< "Y dimensions: " << Y.dimensions() << nl;
    // Info<< "phi dimensions: " << phi.dimensions() << nl;
    // Info<< "D dimensions: " << D.dimensions() << nl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // #include "CourantNo.H"
        // // #include "alphaCourantNo.H"
        // #include "setDeltaT.H"

        while (simple.correctNonOrthogonal())
        {

            C.max(0.); 
            fvScalarMatrix ConvectionEqn //(C=Concentration=rho*wf/mw)
            (
                fvm::ddt(C) 
                + fvm::div(phi,C) 
                - fvm::laplacian(D,C)
                + fvm::Sp(mu*B/(Ks + C)/Y, C) // implicit treatment of source term
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

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
