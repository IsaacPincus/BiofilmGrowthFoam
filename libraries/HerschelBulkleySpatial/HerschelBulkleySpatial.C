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

\*---------------------------------------------------------------------------*/

#include "HerschelBulkleySpatial.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(HerschelBulkleySpatial, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HerschelBulkleySpatial,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::HerschelBulkleySpatial::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());
    volScalarField tau0 = U_.mesh().lookupObject<volScalarField>("tau0"); 

    return
    (
        min
        (
            nu0_,
            (tau0 + k_*rtone*pow(tone*sr(), n_))
           /(max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::HerschelBulkleySpatial::HerschelBulkleySpatial
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    HerschelBulkleySpatialCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    k_("k", dimViscosity, HerschelBulkleySpatialCoeffs_),
    n_("n", dimless, HerschelBulkleySpatialCoeffs_),
    nu0_("nu0", dimViscosity, HerschelBulkleySpatialCoeffs_),
    tau0_
    (
        IOobject
        (
            "tau0",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh()
    ),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::HerschelBulkleySpatial::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    HerschelBulkleySpatialCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    HerschelBulkleySpatialCoeffs_.lookup("k") >> k_;
    HerschelBulkleySpatialCoeffs_.lookup("n") >> n_;
    // HerschelBulkleySpatialCoeffs_.lookup("tau0") >> tau0_;
    HerschelBulkleySpatialCoeffs_.lookup("nu0") >> nu0_;

    return true;
}


// ************************************************************************* //
