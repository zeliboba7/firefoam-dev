/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "filmHeightInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "surfaceFilmModel.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    deltafName_("deltaf")
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    deltafName_(ptf.deltafName_)
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    deltafName_(dict.lookupOrDefault<word>("deltaf", "deltaf"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& fhivpvf
)
:
    fixedValueFvPatchVectorField(fhivpvf),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& fhivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fhivpvf, iF),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const fvPatchField<scalar>& deltafp =
        patch().lookupPatchField<volScalarField, scalar>(deltafName_);

    vectorField n = patch().nf();
    const scalarField& magSf = patch().magSf();

    static bool first=true;
    if(!first){
        first=false;
//###################Nusselt velocity######################
    if(1){
    const objectRegistry& obr=db().parent().lookupObject<objectRegistry>
        (
         "region0"
        );
    const regionModels::surfaceFilmModels::kinematicSingleLayer& film2 =
        refCast<const regionModels::surfaceFilmModels::kinematicSingleLayer>(obr.lookupObject<regionModels::surfaceFilmModels::surfaceFilmModel>
        (
         "surfaceFilmProperties"
        )
        );
    volVectorField gTan2=film2.gTan();
    const label patchI = patch().index();
    gTan2.boundaryField()[patchI];
    const fvPatchField<scalar>& mup =
        patch().lookupPatchField<volScalarField, scalar>("muf");
    scalarField nu(deltafp.size()); 
    nu=mup/rhop;
    scalarField gT(deltafp.size()); 
    gT=mag(gTan2.boundaryField()[patchI]);
    /*operator==(0.1*n);*/
    /*Info << "gT " << gT<<endl;*/
    /*Info << "deltafp " << deltafp<<endl;*/
    /*Info << "nu " << nu<<endl;*/
    /*Info << "n " << n<<endl;*/
    operator==(-(gT*deltafp*deltafp/3.0/nu)*n);
    }
//###################Nusselt velocity######################

    }
    else{
        operator==(deltafp*n*phip/(rhop*magSf*sqr(deltafp) + ROOTVSMALL));
        /*operator==(deltafp*n*phip/(rhop*magSf*sqr(deltafp) + SMALL));*/
    }
    

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::filmHeightInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "deltaf", "deltaf", deltafName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        filmHeightInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
