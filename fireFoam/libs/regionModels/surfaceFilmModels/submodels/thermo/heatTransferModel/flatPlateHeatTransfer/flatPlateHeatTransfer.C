/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "flatPlateHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "thermoSingleLayerFM.H" //kvm

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(flatPlateHeatTransfer, 0);
        addToRunTimeSelectionTable
        (
            heatTransferModel,
            flatPlateHeatTransfer,
            dictionary
        );
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::surfaceFilmModels::flatPlateHeatTransfer::Nu
(
    const scalar Re,
    const scalar Pr
) const
{
    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Pr);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Pr);
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::flatPlateHeatTransfer::flatPlateHeatTransfer
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    heatTransferModel(typeName, owner, dict),
    L_(readScalar(coeffs_.lookup("L"))),
    htcConvFilm_
    (
        IOobject
        (
            "htcConvGas",
            owner.time().timeName(),
            owner.film(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        owner.film(),
        dimensionedScalar("zero", dimMass/pow3(dimTime)/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::flatPlateHeatTransfer::~flatPlateHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::flatPlateHeatTransfer::correct()
{
    const thermoSingleLayerFM& film = refCast<const thermoSingleLayerFM>(owner_); //kvm

    // retrieve fields from film model
    /*const scalarField& delta = film.delta();*/
    //const scalarField& T = film.T();
    //const scalarField& Tw = film.Tw();
    /*const scalarField& rho = film.rho();*/
    /*const scalarField& mu = film.mu();*/
    //const scalarField& magSf = film.magSf();
    const scalarField kappaInf = film.kappaPrimary();
    const vectorField dU = film.UPrimary() - film.Us();
    /*const scalarField& TInf = film.TPrimary();*/
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& muInf = film.muPrimary();
    const scalarField& cpInf = film.cpPrimary();

    forAll(htcConvFilm_, cellI)
    {
        // local temperature - impose upper/lower limit for stability
        /*const scalar Tloc = min(400.0, max(200.0, T[cellI])); */

        // Primary region thermal conductivity
        const scalar kappaInfc = kappaInf[cellI];

        // Primary region density [kg/m3]
        const scalar rhoInfc = rhoInf[cellI];

        // Primary region viscosity [Pa.s]
        const scalar muInfc = muInf[cellI];

        // Primary region heat capacity [J/kg/K]
        const scalar cpInfc = cpInf[cellI];

        // Prandtl number
        const scalar Pr = cpInfc*muInfc/(kappaInfc + ROOTVSMALL);

        // Reynolds number
        const scalar Re = rhoInfc*mag(dU[cellI])*L_/muInfc;

        // Nusselt number
        const scalar Nu = this->Nu(Re, Pr);

        // heat transfer coefficient [W/m2/K]
        const scalar htc = Nu*kappaInfc/(L_+ROOTVSMALL);

        // if Re -> 0, then the htc should revert to pure conduction, k/dx
        // TODO use the wall distance to represent dx
        /*const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");*/
        /*const scalarField& y = rasModel.y()[patchI];*/
        const scalar dx = 0.01; // for now, assume 1 cm
        const scalar htcCond=kappaInfc/dx;

        htcConvFilm_[cellI]=max(htc,htcCond);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::flatPlateHeatTransfer::h() const
{
    return htcConvFilm_;
}



// ************************************************************************* //
