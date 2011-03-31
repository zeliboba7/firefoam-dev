/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "uniformDensityHydrostaticTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::
uniformDensityHydrostaticTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    rhoRefValue_(0.0),
    pRefValue_(p.size(), 0.0),
    pRefPoint_(vector::zero),
    phiName_("phi"),
    UName_("U"),
    rhoName_("rho")
{}


Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::
uniformDensityHydrostaticTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    rhoRefValue_(readScalar(dict.lookup("rhoRefValue"))),
    pRefValue_("pRefValue", dict, p.size()),
    pRefPoint_(dict.lookup("pRefPoint")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::
uniformDensityHydrostaticTotalPressureFvPatchScalarField
(
    const uniformDensityHydrostaticTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    rhoRefValue_(ptf.rhoRefValue_),
    pRefValue_(ptf.pRefValue_, mapper),
    pRefPoint_(ptf.pRefPoint_),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::
uniformDensityHydrostaticTotalPressureFvPatchScalarField
(
    const uniformDensityHydrostaticTotalPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    rhoRefValue_(tppsf.rhoRefValue_),
    pRefValue_(tppsf.pRefValue_),
    pRefPoint_(tppsf.pRefPoint_),
    phiName_(tppsf.phiName_),
    UName_(tppsf.UName_),
    rhoName_(tppsf.rhoName_)
{}


Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::
uniformDensityHydrostaticTotalPressureFvPatchScalarField
(
    const uniformDensityHydrostaticTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    rhoRefValue_(tppsf.rhoRefValue_),
    pRefValue_(tppsf.pRefValue_),
    pRefPoint_(tppsf.pRefPoint_),
    phiName_(tppsf.phiName_),
    UName_(tppsf.UName_),
    rhoName_(tppsf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    pRefValue_.autoMap(m);
}


void Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const uniformDensityHydrostaticTotalPressureFvPatchScalarField& tiptf =
        refCast
        <
            const uniformDensityHydrostaticTotalPressureFvPatchScalarField
        >(ptf);

    pRefValue_.rmap(tiptf.pRefValue_, addr);
}


void Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::
updateCoeffs(const vectorField& Up)
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& rho =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    scalarField p0 =
    (
        pRefValue_
      + rhoRefValue_*((g.value() & patch().Cf()) - (g.value() & pRefPoint_))
    );

    operator==(p0 - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));


    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::
updateCoeffs()
{
    updateCoeffs(patch().lookupPatchField<volVectorField, vector>(UName_));
}


void Foam::uniformDensityHydrostaticTotalPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("rhoRefValue") << rhoRefValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("pRefPoint") << pRefPoint_ << token::END_STATEMENT << nl;
    pRefValue_.writeEntry("pRefValue", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformDensityHydrostaticTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
