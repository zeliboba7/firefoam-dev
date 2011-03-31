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

#include "alphaSgsBuoyantWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "muSgsBuoyantWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphaSgsBuoyantWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphaSgsBuoyantWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphaSgsBuoyantWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void alphaSgsBuoyantWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphaSgsBuoyantWallFunctionFvPatchScalarField::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphaSgsBuoyantWallFunctionFvPatchScalarField::
alphaSgsBuoyantWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    beta_(1.0/293.0),
    magG_(9.80665)
{
    checkType();
    //read();
}


alphaSgsBuoyantWallFunctionFvPatchScalarField::
alphaSgsBuoyantWallFunctionFvPatchScalarField
(
    const alphaSgsBuoyantWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    beta_(ptf.beta_),
    magG_(ptf.magG_)
{}


alphaSgsBuoyantWallFunctionFvPatchScalarField::
alphaSgsBuoyantWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.0/293.0)),
    magG_(dict.lookupOrDefault<scalar>("magG", 9.80665))
{
    checkType();
    //read();
}


alphaSgsBuoyantWallFunctionFvPatchScalarField::
alphaSgsBuoyantWallFunctionFvPatchScalarField
(
    const alphaSgsBuoyantWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    beta_(tppsf.beta_),
    magG_(tppsf.magG_)
{
    checkType();
}


alphaSgsBuoyantWallFunctionFvPatchScalarField::
alphaSgsBuoyantWallFunctionFvPatchScalarField
(
    const alphaSgsBuoyantWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    beta_(tppsf.beta_),
    magG_(tppsf.magG_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void alphaSgsBuoyantWallFunctionFvPatchScalarField::read()
{

    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");
    beta_ = readScalar(sgs.thermo().lookup("beta"));
    const IOdictionary& environmentalProperties =
        db().lookupObject<IOdictionary>
        (
            "environmentalProperties"
        );
    dimensionedVector g(environmentalProperties.lookup("g"));
    magG_ = mag(g).value();
}
*/

void alphaSgsBuoyantWallFunctionFvPatchScalarField::evaluate
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

    // Field data
    const label patchI = patch().index();

    const scalarField& alphaw = sgs.alpha().boundaryField()[patchI];
    scalarField& alphaSgsw = *this;

    const scalarField& rhow = sgs.rho().boundaryField()[patchI];

    const fvPatchScalarField& Tw = thermo.T().boundaryField()[patchI];

//    const scalarField Cpw = thermo.Cp()().boundaryField()[patchI];
    const scalarField Cpw = thermo.Cp(Tw, patchI);

    const scalarField T = Tw.patchInternalField();

    const scalarField magGradTw = max(mag(Tw.snGrad()), VSMALL);
//    const scalarField magGradTw = mag(Tw.snGrad());

    const scalarField& ry = patch().deltaCoeffs();


    // Populate boundary values
    forAll(alphaSgsw, faceI)
    {
        // Initial value from the previous iteration
        scalar qw =  (alphaw[faceI] + alphaSgsw[faceI])*Cpw[faceI]*magGradTw[faceI];

        scalar gBeta = pow(magG_*beta_,0.25);
        scalar aCoef = (Tw[faceI]-T[faceI])
                     * gBeta 
                     * pow(alphaw[faceI]*sqr(rhow[faceI])*pow(Cpw[faceI],3),0.25);
        scalar bCoef = 1.0 / ry[faceI] 
                     * gBeta 
                     * sqrt(rhow[faceI])
                     / pow(Cpw[faceI]*pow(alphaw[faceI],3),0.25);

//        if (yPlus[faceI] > 2.2)
//        {
            label iter = 0;
            scalar err = GREAT;

            do
            {
                scalar f =
                      aCoef*pow(qw, -0.75)
                    - 0.427*log(bCoef*pow(qw,0.25))
                    - 1.93;

                scalar df =
                    - 0.75*aCoef*pow(qw,-1.75)
                    - 0.427/4.0/qw;

                scalar qwNew = qw - f/df;
                err = mag((qw - qwNew)/qw);
                qw = qwNew;

            } while (qw > VSMALL && err > tolerance_ && ++iter < maxIters_);
//        }
//        else
//        {
//            qw = alphaw[faceI]*Cpw[faceI]*magGradTw[faceI];
//        }


        scalar alphaEff = qw/Cpw[faceI]/magGradTw[faceI];
        alphaSgsw[faceI] = max(0.0, alphaEff - alphaw[faceI]);

        if (debug)
        {
            Info<< "    alphaEff       = " << alphaEff << nl
                << "    alphaw         = " << alphaw[faceI] << nl
                << "    alphaSgsw      = " << alphaSgsw[faceI] << nl
//                << "    yPlus          = " << yPlus[faceI] << nl
                << "    Tw             = " << Tw[faceI] << nl
//                << "    Tc             = " << TcEff << nl
                << "    T              = "  << T[faceI] << nl
//                << "    delta          = "  << delta[faceI] << nl
                << "    ry             = "  << ry[faceI] << nl
                << "    magGradTw      = "  << magGradTw[faceI] << nl
                << "    qw             = "  << qw << nl
                << endl;
        }

    }

//  Grab the muSgs patch field using generic/base type
    const fvPatchScalarField& muSgsPatchField =
       patch().lookupPatchField<volScalarField, scalar>("muSgs"); 

    // Perform the type checking
    if (!isA<muSgsBuoyantWallFunctionFvPatchScalarField>(muSgsPatchField))
    {
        FatalErrorIn("alphaSgsBuoyantWallFunctionFvPatchScalarField::evaluate()")
            << "Invalid boundary condition for muSgs" << nl
            << "use muSgsBuoyantWallFunction" << nl
            << endl << abort(FatalError);
    }

//    muSgsBuoyantWallFunctionFvPatchScalarField& muSgsw = 
//       const_cast<muSgsBuoyantWallFunctionFvPatchScalarField&>(muSgsPatchField); 

    const muSgsBuoyantWallFunctionFvPatchScalarField& muSgsw = 
       refCast<const muSgsBuoyantWallFunctionFvPatchScalarField>(muSgsPatchField); 

    muSgsBuoyantWallFunctionFvPatchScalarField& muSgs = 
       const_cast<muSgsBuoyantWallFunctionFvPatchScalarField&>(muSgsw); 

    muSgs.evaluateInAlphaSgs();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphaSgsBuoyantWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
