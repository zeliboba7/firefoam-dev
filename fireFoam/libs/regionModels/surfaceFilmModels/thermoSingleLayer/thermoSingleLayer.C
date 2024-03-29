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

#include "thermoSingleLayer.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "directMappedFieldFvPatchField.H"
#include "mapDistribute.H"

// Sub-models
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"

#include "stdio.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayer, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList thermoSingleLayer::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchI)
    {
        if (bTypes[patchI] == directMappedFieldFvPatchField<scalar>::typeName)
        {
            bTypes[patchI] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayer::read()
{
    // no additional properties to read
    return kinematicSingleLayer::read();
}


void thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::resetPrimaryRegionSourceTerms()" << endl;
    }

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar("zero", hsSp_.dimensions(), 0.0);
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::resetPrimaryRegionSourceTerms()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayer::correctThermoFields()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::correctThermoFields()" << endl;
    }
    switch (thermoModel_)
    {
        case tmConstant:
        {
            rho_ == dimensionedScalar(coeffs_.lookup("rho0"));
            mu_ == dimensionedScalar(coeffs_.lookup("mu0"));
            sigma_ == dimensionedScalar(coeffs_.lookup("sigma0"));
            Cp_ == dimensionedScalar(coeffs_.lookup("Cp0"));
            kappa_ == dimensionedScalar(coeffs_.lookup("kappa0"));

            break;
        }
        case tmSingleComponent:
        {
            const liquid& liq = thermo_.liquids().properties()[liquidId_];
            forAll(rho_, cellI)
            {
                const scalar T = min(400.0,max(273.15,T_[cellI]));//TODO: add to code pack
                const scalar p = pPrimary_[cellI];
                rho_[cellI] = liq.rho(p, T);
                mu_[cellI] = liq.mu(p, T);
                sigma_[cellI] = liq.sigma(p, T);
                Cp_[cellI] = liq.cp(p, T);
                kappa_[cellI] = liq.K(p, T);
            }

            rho_.correctBoundaryConditions();
            mu_.correctBoundaryConditions();
            sigma_.correctBoundaryConditions();
            Cp_.correctBoundaryConditions();
            kappa_.correctBoundaryConditions();

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void thermoSingleLayer::"
                "correctThermoFields()"
            )   << "Unknown thermoType enumeration" << abort(FatalError);
        }
    }
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::correctThermoFields()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayer::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    forAll(T_.boundaryField(), patchI)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchI];
        if (isA<directMappedFieldFvPatchField<scalar> >(Tp))
        {
            hs_.boundaryField()[patchI] == hs(Tp, patchI);
        }
    }
}


void thermoSingleLayer::updateSurfaceTemperatures()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::updateSurfaceTemperatures()" << endl;
    }
    correctHsForMappedT();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchI = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchI];
        if(pyrCoupled_){
            UIndirectList<scalar>(Tw_, pp.faceCells()) =
            //T_.boundaryField()[patchI];
            pyrTemperaturePtr_->boundaryField()[patchI];
        }
        else{
            UIndirectList<scalar>(Tw_, pp.faceCells()) =
            T_.boundaryField()[patchI];
            //pyrTemperaturePtr_->boundaryField()[patchI];
        }
    }
    Tw_.correctBoundaryConditions();

    // Update heat transfer to wall (used in film/pyrolysis coupling)
    // heat flow out of film is positive
    qFilmToWall_ = htcw_->h()*(T_ - Tw_);
    forAll(qFilmToWall_,i){
        if(delta_[i]<1e-8){
            qFilmToWall_[i]=0.0;
        }
    }

    qFilmToWall_.correctBoundaryConditions();

    // Update heat transfer from gas phase (used in diagnostics)
    // heat flow out of film is positive
    qGasToFilm_ = htcs_->h()*(T_ - TPrimary_);
    forAll(qGasToFilm_,i){
        if(delta_[i]<1e-8){
            qGasToFilm_[i]=0.0;
        }
    }
    qGasToFilm_.correctBoundaryConditions();

    // Update film surface temperature
    Ts_ = T_;
    Ts_.correctBoundaryConditions();

    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::updateSurfaceTemperatures()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayer::transferPrimaryRegionThermoFields()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::transferPrimaryRegionThermoFields()" << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::transferPrimaryRegionThermoFields()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayer::transferPrimaryRegionSourceFields()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::transferPrimaryRegionSourceFields()" << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    // Retrieve the source fields from the primary region (spray) via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    //   kvm, this is coming from hsSpPrimary_
    hsSp_.correctBoundaryConditions();

    // Convert accummulated source terms into per unit area per unit time
    // Note: boundary values will still have original (neat) values
    const scalar deltaT = time_.deltaTValue();
    hsSp_.field() /= magSf()*deltaT;

    // Apply enthalpy source as difference between incoming and actual states
    hsSp_ -= rhoSp_*hs_;
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::transferPrimaryRegionSourceFields()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayer::updateSubmodels()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::updateSubmodels()" << endl;
    }

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassPCTrans_,
        primaryEnergyPCTrans_
    );

    //for diagnostics
    primaryMassPCTrans_.correctBoundaryConditions();
    primaryEnergyPCTrans_.correctBoundaryConditions();

    // Update radiation
    // push radiation source into qRad_ for use in standardPhaseChange (there must be a better way to do this!)
    radiation_->correct();
    const scalarField& qRad = radiation_->Shs();
    forAll(qRad_,i){
        qRad_[i]=qRad[i];
    }

    // Update kinematic sub-models
    kinematicSingleLayer::updateSubmodels();

    // Update source fields
    // kvm, why was this originally -= ?, A: because spray is added first
    hsSp_ += primaryEnergyPCTrans_/magSf()/time().deltaT();
    rhoSp_ += primaryMassPCTrans_/magSf()/time().deltaT();
    // vapor recoil pressure
    const dimensionedScalar rhov("rhov",dimDensity,1.0); //gas phase density (I need to get this from gas-phase instead)
    pSp_ -= pow(primaryMassPCTrans_/magSf()/time_.deltaT(),2)/2.0/rhov;

    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::updateSubmodels()" << endl;
        tabSubtract();
    }
}


tmp<fvScalarMatrix> thermoSingleLayer::q(volScalarField& hs) const
{
    dimensionedScalar Tstd("Tstd", dimTemperature, 298.15);

    return
    (
      - fvm::Sp(htcs_->h()/Cp_, hs) - htcs_->h()*(Tstd - TPrimary_)
      - fvm::Sp(htcw_->h()/Cp_, hs) - htcw_->h()*(Tstd - Tw_)
    );
}


void thermoSingleLayer::solveEnergy()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::solveEnergy()" << endl;
    }

    updateSurfaceTemperatures();

    dimensionedScalar hs0("SMALL", hs_.dimensions(), ROOTVSMALL);

    solve
    (
        fvm::ddt(deltaRho_, hs_)
      + fvm::div(phi_, hs_)
     ==
       -hsSp_ 
      + q(hs_) //convective heat transfer
      + radiation_->Shs() //radiative heat transfer
      - rhoSp_*hs_ 
    );

    correctThermoFields();
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::solveEnergy()" << endl;
        tabSubtract();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const bool readFields
)
:
    kinematicSingleLayer(modelType, mesh, g, false),
    thermo_(mesh.lookupObject<SLGThermo>("SLGThermo")),
    liquidId_(thermo_.liquidId(coeffs_.lookup("liquid"))),
    Cp_
    (
        IOobject
        (
            "Cp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("Cp", dimEnergy/dimMass/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar
        (
            "kappa",
            dimEnergy/dimTime/dimLength/dimTemperature,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    T_
    (
        IOobject
        (
            "Tf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Ts_
    (
        IOobject
        (
            "Tsf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    Tw_
    (
        IOobject
        (
            "Twf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        zeroGradientFvPatchScalarField::typeName
    ),
    hs_
    (
        IOobject
        (
            "hsf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimMass, 0.0),
        hsBoundaryTypes()
    ),
    qGasToFilm_
    (
     IOobject
     (
      "qGasToFilm",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qFilmToWall_
    (
     IOobject
     (
      "qWall",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),

    primaryMassPCTrans_
    (
        IOobject
        (
            "primaryMassPCTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    primaryEnergyPCTrans_
    (
        IOobject
        (
            "primaryEnergyPCTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    qRad_
    (
     IOobject
     (
      "qRad",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0)
    ),

    hsSp_
    (
        IOobject
        (
            "hsSp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    hsSpPrimary_
    (
        IOobject
        (
            hsSp_.name(), // must have same name as hSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar("zero", hsSp_.dimensions(), 0.0)
    ),

    TPrimary_
    (
        IOobject
        (
            "T", // same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimTemperature, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    YPrimary_(),

    htcs_
    (
        heatTransferModel::New(*this, coeffs().subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs().subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs())),
    radiation_(filmRadiationModel::New(*this, coeffs()))
{
    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar("zero", dimless, 0.0),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();
        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_);
        deltaRho_ == delta_*rho_;
        phi_ = fvc::interpolate(deltaRho_*U_) & regionMesh().Sf();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayer::addSources
(
    const label patchI,
    const label faceI,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchI,
        faceI,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::addSources()" << endl;
        Info<<tab.c_str()<< "    energy   = " << energySource << nl << endl;
    }

    hsSpPrimary_.boundaryField()[patchI][faceI] -= energySource; //hsSpPrimary_ is just a temporary placeholder in the gas phase mesh for spray boundary values
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::addSources()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayer::preEvolveRegion()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::preEvolveRegion()" << endl;
    }

    surfaceFilmModel::preEvolveRegion();//kvm, added to map T_pyrolysis to film model

//    correctHsForMappedT();

    kinematicSingleLayer::preEvolveRegion();

    // Update phase change
    primaryMassPCTrans_ == dimensionedScalar("zero", dimMass, 0.0);
    primaryEnergyPCTrans_ == dimensionedScalar("zero", dimEnergy, 0.0);
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::preEvolveRegion()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayer::evolveRegion()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayer::evolveRegion()" << endl;
    }

    //kvm, if updateSubmodels occurs after solveContinuity, then vaporization and separation don't get called until after deltaRho_ has been calculated, resulting in errors
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=0; oCorr<nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu(this->pu());

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp(this->pp());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Solve energy for hs_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    T_ == T(hs_);

#include "diagnostics.H"

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;


    // Update film wall and surface velocities
    updateSurfaceVelocities();

    // Update film wall and surface temperatures
    // only need to call once?
    // updateSurfaceTemperatures();

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();

    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayer::evolveRegion()" << endl;
        tabSubtract();
    }
}


const volScalarField& thermoSingleLayer::Cp() const
{
    return Cp_;
}


const volScalarField& thermoSingleLayer::kappa() const
{
    return kappa_;
}


const volScalarField& thermoSingleLayer::T() const
{
    return T_;
}


const volScalarField& thermoSingleLayer::Ts() const
{
    return Ts_;
}


const volScalarField& thermoSingleLayer::Tw() const
{
    return Tw_;
}

tmp<DimensionedField<scalar, volMesh> > thermoSingleLayer::qRad() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        qRad_
    );
}


const volScalarField& thermoSingleLayer::hs() const
{
    return hs_;
}


tmp<volScalarField> thermoSingleLayer::primaryMassTrans() const
{
    return primaryMassPCTrans_;
}


void thermoSingleLayer::info() const
{
    kinematicSingleLayer::info();

    Info<< indent << "min/max(T)         = " << min(T_).value() << ", "
        << max(T_).value() << nl;

    phaseChange_->info(Info);
}


tmp<DimensionedField<scalar, volMesh> > thermoSingleLayer::Srho() const
{
    tmp<DimensionedField<scalar, volMesh> > tSrho
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Srho",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    scalarField& Srho = tSrho();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchI = intCoupledPatchIDs()[i];
        const mapDistribute& distMap = mappedPatches_[filmPatchI].map();

        scalarField patchMass =
            primaryMassPCTrans_.boundaryField()[filmPatchI];
        distMap.distribute(patchMass);

        const label primaryPatchI = primaryPatchIDs()[i];
        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchI].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] = patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<DimensionedField<scalar, volMesh> > thermoSingleLayer::Srho
(
    const label i
) const
{
    const label vapId =
        thermo_.carrierId(thermo_.liquids().components()[liquidId_]);

    tmp<DimensionedField<scalar, volMesh> > tSrho
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Srho(" + Foam::name(i) + ")",
                time_.timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchI = intCoupledPatchIDs_[i];
            const mapDistribute& distMap = mappedPatches_[filmPatchI].map();

            scalarField patchMass =
                primaryMassPCTrans_.boundaryField()[filmPatchI];
            distMap.distribute(patchMass);

            const label primaryPatchI = primaryPatchIDs()[i];
            const unallocLabelList& cells =
                primaryMesh().boundaryMesh()[primaryPatchI].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] = patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<DimensionedField<scalar, volMesh> > thermoSingleLayer::Sh() const
{
    tmp<DimensionedField<scalar, volMesh> > tSh
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "thermoSingleLayer::Sh",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
/*
    phase change energy fed back into the film...

    scalarField& Sh = tSh();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchI = intCoupledPatchIDs_[i];
        const mapDistribute& distMap = mappedPatches_[filmPatchI].map();

        scalarField patchEnergy =
            primaryEnergyPCTrans_.boundaryField()[filmPatchI];
        distMap.distribute(patchEnergy);

        const label primaryPatchI = primaryPatchIDs()[i];
        const unallocLabelList& cells =
            primaryMesh().boundaryMesh()[primaryPatchI].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }
*/
    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
