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

#include "surfaceFilmModel.H"
#include "fvMesh.H"
#include "directMappedFieldFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<
    Foam::regionModels::surfaceFilmModels::surfaceFilmModel::thermoModelType,
    2
>::names[] =
{
    "constant",
    "singleComponent"
};


const Foam::NamedEnum
<
    Foam::regionModels::surfaceFilmModels::surfaceFilmModel::thermoModelType,
    2
>
Foam::regionModels::surfaceFilmModels::surfaceFilmModel::thermoModelTypeNames_;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceFilmModel, 0);
defineRunTimeSelectionTable(surfaceFilmModel, mesh);

void surfaceFilmModel::constructMeshObjects()
{
    // construct pyrTemperature field if coupled to pyr model
    if (pyrCoupled_)
    {

        pyrTemperaturePtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pyrolysisTemperature",
                    time_.timeName(),
                    regionMesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                regionMesh()
            )
        );

        const volScalarField& pyrTemperature = pyrTemperaturePtr_();

        bool foundCoupledPatch = false;
        forAll(pyrTemperature.boundaryField(), patchI)
        {
            const fvPatchField<scalar>& fvp = pyrTemperature.boundaryField()[patchI];
            if (isA<directMappedFieldFvPatchField<scalar> >(fvp))
            {
                foundCoupledPatch = true;
                break;
            }
        }

        if (!foundCoupledPatch)
        {
            WarningIn("void pyrolysisModels::constructMeshObjects()")
                << "pyrCoupled flag set to true, but no "
                << directMappedFieldFvPatchField<scalar>::typeName
                << " patches found on " << pyrTemperature.name() << " field"
                << endl;
        }

    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool surfaceFilmModel::read()
{
    if (singleLayerRegion::read())
    {
        diagnostics_ = coeffs_.lookupOrDefault<Switch>("diagnostics",false);
        pyrCoupled_ = coeffs_.lookupOrDefault<Switch>("pyrolysisCoupled",false);
        thermoModel_ =
            thermoModelTypeNames_.read(coeffs_.lookup("thermoModel"));
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceFilmModel::surfaceFilmModel(const fvMesh& mesh)
:
    singleLayerRegion(mesh),
    pyrCoupled_(false),
    diagnostics_(false),
    g_(vector::zero),
    thermoModel_(tmConstant)
{}


surfaceFilmModel::surfaceFilmModel
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    singleLayerRegion(mesh, "surfaceFilm", modelType),
    pyrCoupled_(false),
    diagnostics_(false),
    g_(g),
    thermoModel_(tmConstant)
{
    if (active_)
    {
        read();
        constructMeshObjects();//kvm
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

surfaceFilmModel::~surfaceFilmModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void surfaceFilmModel::preEvolveRegion()
{
    if (pyrCoupled_)
    {
        pyrTemperaturePtr_->correctBoundaryConditions();  
        if(time_.outputTime()){
            pyrTemperaturePtr_->write();
        }
    }
}


Foam::scalar surfaceFilmModel::CourantNumber() const
{
    return ROOTVSMALL;
}


tmp<DimensionedField<scalar, volMesh> >
surfaceFilmModel::qRad() const
{
    FatalErrorIn("const volScalarField& noFilm::qRad() const")
        << "qRad field not available for " << type() << abort(FatalError);

    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "surfaceFilmModel::qRad()",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0)
        )
    );
}

tmp<DimensionedField<scalar, volMesh> >
surfaceFilmModel::qWall() const
{
    FatalErrorIn("const volScalarField& noFilm::qWall() const")
        << "qWall field not available for " << type() << abort(FatalError);

    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "surfaceFilmModel::qWall()",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0)
        )
    );
}

const volScalarField& surfaceFilmModel::omega() const
{
    FatalErrorIn("const volScalarField& noFilm::omega() const")
        << "omega field not available for " << type() << abort(FatalError);

    return reinterpret_cast<const volScalarField&>(null);
}


tmp<DimensionedField<scalar, volMesh> > surfaceFilmModel::Srho() const
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > surfaceFilmModel::Srho() const"
    )

    return tmp<DimensionedField<scalar, volMesh> >(NULL);
}


tmp<DimensionedField<scalar, volMesh> >
surfaceFilmModel::Srho(const label) const
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > surfaceFilmModel::Srho"
        "(const label) const"
    )

    return tmp<DimensionedField<scalar, volMesh> >(NULL);
}


tmp<DimensionedField<scalar, volMesh> > surfaceFilmModel::Sh() const
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > surfaceFilmModel::Sh() const"
    )

    return tmp<DimensionedField<scalar, volMesh> >(NULL);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
