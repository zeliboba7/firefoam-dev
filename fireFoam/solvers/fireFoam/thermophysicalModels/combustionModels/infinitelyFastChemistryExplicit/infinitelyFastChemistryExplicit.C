/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "infinitelyFastChemistryExplicit.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(infinitelyFastChemistryExplicit, 0);
    addToRunTimeSelectionTable
    (
        combustionModel,
        infinitelyFastChemistryExplicit,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::infinitelyFastChemistryExplicit::infinitelyFastChemistryExplicit
(
    const dictionary& combustionProps,
    hsCombustionThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const surfaceScalarField& phi,
    const volScalarField& rho
)
:
    combustionModel(typeName, combustionProps, thermo, turbulence, phi, rho),
    C_(readScalar(coeffs_.lookup("C"))),
    singleMixture_
    (
        dynamic_cast<singleStepReactingMixture<gasThermoPhysics>&>(thermo)
    ),
    wFuelNorm_
    (
        IOobject
        (
            "wFuelNorm",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::infinitelyFastChemistryExplicit::~infinitelyFastChemistryExplicit()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::infinitelyFastChemistryExplicit::correct()
{
    singleMixture_.fresCorrect();

    const label fuelI = singleMixture_.fuelIndex();

    const volScalarField& YFuel = thermo_.composition().Y()[fuelI];

    const dimensionedScalar s = singleMixture_.s();

    if (thermo_.composition().contains("O2"))
    {
        const volScalarField& YO2 = thermo_.composition().Y("O2");
        wFuelNorm_ == rho_/(mesh_.time().deltaT()*C_)*min(YFuel, YO2/s.value());
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::infinitelyFastChemistryExplicit::R(volScalarField& Y) const
{
    const label specieI = thermo_.composition().species()[Y.name()];

    const label fNorm = singleMixture_.specieProd()[specieI];

    const volScalarField fres = singleMixture_.fres(specieI);

    const volScalarField wSpecie =
        wFuelNorm_*singleMixture_.specieStoichCoeffs()[specieI];
//      / max(fNorm*(Y - fres), SMALL);

    //return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
    return wSpecie + fvm::Sp(dimensionedScalar("zero", wSpecie.dimensions(), scalar(0.0)), Y);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::infinitelyFastChemistryExplicit::dQ() const
{
    const label fuelI = singleMixture_.fuelIndex();
    volScalarField& YFuel = thermo_.composition().Y(fuelI);

    return -singleMixture_.qFuel()*(R(YFuel) & YFuel);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::infinitelyFastChemistryExplicit::wFuelNorm() const
{
    return wFuelNorm_;
}


bool Foam::combustionModels::infinitelyFastChemistryExplicit::read
(
    const dictionary& combustionProps
)
{
    combustionModel::read(combustionProps);
    coeffs_.lookup("C") >> C_ ;

    return true;
}


// ************************************************************************* //
