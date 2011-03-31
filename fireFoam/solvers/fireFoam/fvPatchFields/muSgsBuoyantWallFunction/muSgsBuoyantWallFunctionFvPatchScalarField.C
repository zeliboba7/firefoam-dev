/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "muSgsBuoyantWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void muSgsBuoyantWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "muSgsBuoyantWallFunctionFvPatchScalarField::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muSgsBuoyantWallFunctionFvPatchScalarField::
muSgsBuoyantWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    beta_(1.0/293.0),
    magG_(9.80665),
    Tref_(293.0),
    Prt_(0.9)
{
    checkType();
//    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");
//    sgs.thermo().lookup("beta") >> beta_ ;
//    sgs.thermo().lookup("Tref") >> Tref_;

/*
    const IOdictionary& environmentalProperties =
        db().lookupObject<IOdictionary>
        (
            "environmentalProperties"
        );

    dimensionedVector g(environmentalProperties.lookup("g"));

    magG_ = mag(g).value();
*/
}


muSgsBuoyantWallFunctionFvPatchScalarField::
muSgsBuoyantWallFunctionFvPatchScalarField
(
    const muSgsBuoyantWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    beta_(ptf.beta_),
    magG_(ptf.magG_),
    Tref_(ptf.Tref_),
    Prt_(ptf.Prt_)
{}


muSgsBuoyantWallFunctionFvPatchScalarField::
muSgsBuoyantWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.0/293.0)),
    magG_(dict.lookupOrDefault<scalar>("magG", 9.80665)),
    Tref_(dict.lookupOrDefault<scalar>("Tref", 293.0)),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.9))
{
    checkType();
/*
    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");
    sgs.thermo().lookup("beta") >> beta_ ;
    sgs.thermo().lookup("Tref") >> Tref_;

    const IOdictionary& environmentalProperties =
        db().lookupObject<IOdictionary>
        (
            "environmentalProperties"
        );

    dimensionedVector g(environmentalProperties.lookup("g"));

    magG_ = mag(g).value();
*/
}


muSgsBuoyantWallFunctionFvPatchScalarField::
muSgsBuoyantWallFunctionFvPatchScalarField
(
    const muSgsBuoyantWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    beta_(tppsf.beta_),
    magG_(tppsf.magG_),
    Tref_(tppsf.Tref_),
    Prt_(tppsf.Prt_)
{
    checkType();
}


muSgsBuoyantWallFunctionFvPatchScalarField::
muSgsBuoyantWallFunctionFvPatchScalarField
(
    const muSgsBuoyantWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    beta_(tppsf.beta_),
    magG_(tppsf.magG_),
    Tref_(tppsf.Tref_),
    Prt_(tppsf.Prt_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void muSgsBuoyantWallFunctionFvPatchScalarField::evaluate //dummpy one
(
    const Pstream::commsTypes
)
{
}


void muSgsBuoyantWallFunctionFvPatchScalarField::evaluateInAlphaSgs
(
    const Pstream::commsTypes
)
{
    // Get info from the SGS model
    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");

    const basicThermo& thermo = db().lookupObject<basicThermo>
    (
        "thermophysicalProperties"
    );
    // Wall function constants

    // Field data
    const label patchI = patch().index();

//    const scalar Prt = sgs.Prt().value();
    const scalarField alphaEffw = sgs.alphaEff()().boundaryField()[patchI];
    const scalarField alphaSgsw = sgs.alphaSgs()().boundaryField()[patchI];
    const scalarField alphaw = sgs.alpha().boundaryField()[patchI];

    const scalarField muw = sgs.mu().boundaryField()[patchI];
    scalarField& muSgsw = *this;


    const scalarField& rhow = sgs.rho().boundaryField()[patchI];

    const fvPatchVectorField& Uw = sgs.U().boundaryField()[patchI];
    const vectorField U  = Uw.patchInternalField();
//    const scalarField magUp = max(mag(Uw.patchInternalField() - Uw), VSMALL);
    const scalarField magUp = max(mag(U - Uw), VSMALL);
    const scalarField magFaceGradU = max(mag(Uw.snGrad()),SMALL);


    const fvPatchScalarField& Tw = thermo.T().boundaryField()[patchI];
    const scalarField T = Tw.patchInternalField();
    const scalarField magGradTw = max(mag(Tw.snGrad()),SMALL);

    const scalarField Tc = 
        pow 
        (
            rhow/alphaw/magG_/beta_ * 
            pow 
            (    alphaEffw * magGradTw / rhow,
                 3
            ),
            0.25
        );

    const scalarField delta = 
        pow
        (        
            pow(alphaw,3)/magG_/beta_/alphaEffw/magGradTw,
            0.25
        ) / sqrt(rhow);

    const scalarField& ry = patch().deltaCoeffs();

    const scalarField yPlus = (1.0/ry)/delta;

    const scalarField thetaO = (Tw - Tref_)/Tc;

    // Populate boundary values
    forAll(muSgsw, faceI)
    {
        scalar Pr = muw[faceI]/alphaw[faceI];

        scalar Prat = Pr/Prt_;

        scalar logYPlus = log(yPlus[faceI]);

        scalar alpha = alphaw[faceI] / rhow[faceI]; 

        scalar DuDy = 1.0/Pr*alpha/sqr(delta[faceI])*
               (  Pr*delta[faceI]/alpha*magUp[faceI]  
                - 0.427*Prat*yPlus[faceI]
                * (  0.427*(logYPlus - 2.0) 
                   + 1.93
                   - thetaO[faceI]
                  )
                + 2.27*logYPlus
                - 1.28
               ) / 
               ( 0.49*logYPlus+1.28);

        scalar muEff = DuDy * muw[faceI] / magFaceGradU[faceI];

        muSgsw[faceI] = max(0.0, muEff - muw[faceI]);

        if (debug)
        {
            Info<< "    muEff           = " << muEff << nl
                << "    muw             = " << muw[faceI] << nl
                << "    muSgsw          = " << muSgsw[faceI] << nl
                << "    yPlus           = " << yPlus[faceI] << nl
                << "    Tw              = " << Tw[faceI] << nl
                << "    T               = " << T[faceI] << nl
                << "    Tc              = " << Tc[faceI] << nl
                << "    delta           = " << delta[faceI] << nl
                << "    ry              = " << ry[faceI] << nl
                << "    magUp           = " << magUp[faceI] << nl
//                << "    UcEff           = " << UcEff << nl
                << "    thetaO          = " << thetaO[faceI] << nl
//                << "    DuDyPlus        = " << DuDyPlus << nl
                << endl;
        }

    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    muSgsBuoyantWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
