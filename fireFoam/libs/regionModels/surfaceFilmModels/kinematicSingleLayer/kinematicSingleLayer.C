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

#include "kinematicSingleLayer.H"
#include "fvm.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"
#include "directMappedWallPolyPatch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kinematicSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmModel, kinematicSingleLayer, mesh);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool kinematicSingleLayer::read()
{
    if (surfaceFilmModel::read())
    {
        const dictionary& solution = this->solution().subDict("PISO");
        solution.lookup("momentumPredictor") >> momentumPredictor_;
        solution.lookup("nOuterCorr") >> nOuterCorr_;
        solution.lookup("nCorr") >> nCorr_;
        solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;

        coeffs_.lookup("Cf") >> Cf_;

        return true;
    }
    else
    {
        return false;
    }
}


void kinematicSingleLayer::correctThermoFields()
{
    if (thermoModel_ == tmConstant)
    {
        rho_ == dimensionedScalar(coeffs_.lookup("rho0"));
        mu_ == dimensionedScalar(coeffs_.lookup("mu0"));
        sigma_ == dimensionedScalar(coeffs_.lookup("sigma0"));
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::surfaceFilmModels::kinematicSingleLayer::"
            "correctThermo()"
        )   << "Kinematic surface film must use "
            << thermoModelTypeNames_[thermoModel_] << "thermodynamics" << endl;
    }
}


void kinematicSingleLayer::tabAdd(){
    tab+="\t";
    return;
}
void kinematicSingleLayer::tabSubtract(){
    tab.erase(tab.length()-1,tab.length());
    return;
}
void kinematicSingleLayer::resetPrimaryRegionSourceTerms()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::resetPrimaryRegionSourceTerms()" << endl;
    }

    rhoSpPrimary_ == dimensionedScalar("zero", rhoSp_.dimensions(), 0.0);
    USpPrimary_ == dimensionedVector("zero", USp_.dimensions(), vector::zero);
    pSpPrimary_ == dimensionedScalar("zero", pSp_.dimensions(), 0.0);
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::resetPrimaryRegionSourceTerms()" << endl;
        tabSubtract();
    }
}


void kinematicSingleLayer::transferPrimaryRegionThermoFields()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::"
            << "transferPrimaryRegionThermoFields()" << endl;
    }
    // Update fields from primary region via direct mapped
    // (coupled) boundary conditions
    UPrimary_.correctBoundaryConditions();
    pPrimary_.correctBoundaryConditions();
    rhoPrimary_.correctBoundaryConditions();
    muPrimary_.correctBoundaryConditions();
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::"
            << "transferPrimaryRegionThermoFields()" << endl;
        tabSubtract();
    }
}


void kinematicSingleLayer::transferPrimaryRegionSourceFields()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::"
            << "transferPrimaryRegionSourceFields()" << endl;
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    rhoSp_.correctBoundaryConditions();
    USp_.correctBoundaryConditions();
    pSp_.correctBoundaryConditions();

    // Convert accummulated source terms into per unit area per unit time
    // Note: boundary values will still have original (neat) values
    const scalar deltaT = time_.deltaTValue();
    rhoSp_.field() /= magSf()*deltaT; //at this point, rhoSp_ only contains impingement contribution
    USp_.field() /= magSf()*deltaT;
    pSp_.field() /= magSf()*deltaT;

    /* massImp_ is used for diagnostics */
    massImp_.field() = rhoSp_.field();
    massImp_.correctBoundaryConditions();
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::"
            << "transferPrimaryRegionSourceFields()" << endl;
        tabSubtract();
    }
}


tmp<volScalarField> kinematicSingleLayer::pu()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pu",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pPrimary_                  // pressure (mapped from primary region)
          - pSp_                           // accumulated particle impingement
          //    this term causes odd behavior on non-uniform meshes
          - fvc::laplacian(sigma_, delta_) // surface tension
        )
    );
}


tmp<volScalarField> kinematicSingleLayer::pp()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pp",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -rho_*gNormClipped() // hydrostatic effect only
        )
    );
}


void kinematicSingleLayer::updateSubmodels()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::updateSubmodels()" << endl;
    }

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    // Update source fields
    const dimensionedScalar deltaT = time().deltaT();
    rhoSp_ += cloudMassTrans_/magSf()/deltaT;
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::updateSubmodels()" << endl;
        tabSubtract();
    }
}


void kinematicSingleLayer::continuityCheck()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::continuityCheck()" << endl;
    }
    const volScalarField deltaRho0(deltaRho_);

    solveContinuity();

    if (debug)
    {
        const volScalarField mass(deltaRho_*magSf());
        const dimensionedScalar totalMass =
            fvc::domainIntegrate(mass)
          + dimensionedScalar("SMALL", dimMass*dimVolume, ROOTVSMALL);

        const scalar sumLocalContErr =
            (
                fvc::domainIntegrate(mag(mass - magSf()*deltaRho0))/totalMass
            ).value();

       const scalar globalContErr =
            (
                fvc::domainIntegrate(mass - magSf()*deltaRho0)/totalMass
            ).value();

        cumulativeContErr_ += globalContErr;

        Info<< "Surface film: " << type() << nl
            << "    time step continuity errors: sum local = "
            << sumLocalContErr << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_ << endl;
    }
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::continuityCheck()" << endl;
        tabSubtract();
    }
}


void kinematicSingleLayer::solveContinuity()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::solveContinuity()" << endl;
    }

    solve
    (
        fvm::ddt(deltaRho_)
      + fvc::div(phi_)
     ==
      - rhoSp_
    );

    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::solveContinuity()" << endl;
        tabSubtract();
    }
}


void kinematicSingleLayer::updateSurfaceVelocities()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::updateSurfaceVelocities()" << endl;
    }
    // Push boundary film velocity values into internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchI = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchI];
        UIndirectList<vector>(Uw_, pp.faceCells()) =
            U_.boundaryField()[patchI];
    }
    Uw_ -= nHat()*(Uw_ & nHat());
    Uw_.correctBoundaryConditions();

    // TODO: apply quadratic profile to determine surface velocity
    Us_ = U_;
    Us_.correctBoundaryConditions();
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::updateSurfaceVelocities()" << endl;
        tabSubtract();
    }
}


tmp<fvVectorMatrix> kinematicSingleLayer::tau(volVectorField& U) const
{
    // Calculate shear stress
    //TODO: kvm:  this is not correct, it should be Us_ - U(gas)
    volScalarField Cs("Cs", rho_*Cf_*mag(Us_ - U));
    volScalarField Cw
    (
        "Cw",
        mu_/(0.3333*(delta_ + dimensionedScalar("SMALL", dimLength, SMALL)))
    );
    /*Cw.min(1.0e+06);*/
    /*if (time().outputTime()){*/
        /*Cw.write();*/
    /*}*/
    Cw.min(5.0e+3); //kvm, reduced to help initiate inlet flow 

    return
    (
       - fvm::Sp(Cs, U) + Cs*Us_ // surface contribution
       - fvm::Sp(Cw, U) + Cw*Uw_ // wall contribution
    );
}


tmp<Foam::fvVectorMatrix> kinematicSingleLayer::solveMomentum
(
    const volScalarField& pu,
    const volScalarField& pp
)
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::solveMomentum()" << endl;
    }

    updateSurfaceVelocities();

    // Momentum
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(deltaRho_, U_)
      + fvm::div(phi_, U_)
     ==
      - USp_
      + tau(U_)
      + fvc::grad(sigma_)
      - fvm::SuSp(rhoSp_, U_)
    );

    fvVectorMatrix& UEqn = tUEqn();

    UEqn.relax();

    if (momentumPredictor_)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
              - fvc::interpolate(delta_)
              * (
                    regionMesh().magSf()
                  * (
                        fvc::snGrad(pu, "snGrad(p)")
                      + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
                      + fvc::snGrad(delta_)*fvc::interpolate(pp)
                    )
                  - (fvc::interpolate(rho_*gTan()) & regionMesh().Sf())
                )
            )
        );

        // Remove any patch-normal components of velocity
        U_ -= nHat()*(nHat() & U_);
        U_.correctBoundaryConditions();
    }

    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::solveMomentum()" << endl;
        tabSubtract();
    }
    return tUEqn;
}


void kinematicSingleLayer::solveThickness
(
    const volScalarField& pu,
    const volScalarField& pp,
    const fvVectorMatrix& UEqn
)
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::solveThickness()" << endl;
    }

    volScalarField rUA(1.0/UEqn.A());
    U_ = rUA*UEqn.H();

    surfaceScalarField deltarUAf(fvc::interpolate(delta_*rUA));
    surfaceScalarField rhof(fvc::interpolate(rho_));

    surfaceScalarField phiAdd
    (
        "phiAdd",
        regionMesh().magSf()
      * (
            fvc::snGrad(pu, "snGrad(p)")
          + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
        )
      - (fvc::interpolate(rho_*gTan()) & regionMesh().Sf())
    );
    constrainFilmField(phiAdd, 0.0);

    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(U_*rho_) & regionMesh().Sf())
      - deltarUAf*phiAdd*rhof
    );
    constrainFilmField(phid, 0.0);

    surfaceScalarField ddrhorUAppf
    (
        "deltaCoeff",
        fvc::interpolate(delta_)*deltarUAf*rhof*fvc::interpolate(pp)
    );
//    constrainFilmField(ddrhorUAppf, 0.0);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        // Film thickness equation
        fvScalarMatrix deltaEqn
        (
            fvm::ddt(rho_, delta_)
          + fvm::div(phid, delta_)
          - fvm::laplacian(ddrhorUAppf, delta_)
         ==
          - rhoSp_
        );

        deltaEqn.solve();

        if (nonOrth == nNonOrthCorr_)
        {
            phiAdd +=
                fvc::interpolate(pp)
              * fvc::snGrad(delta_)
              * regionMesh().magSf();

            phi_ == deltaEqn.flux();
        }
    }

    // Bound film thickness by a minimum of zero
    delta_.max(0.0);

    // Update U field
    U_ -= fvc::reconstruct(deltarUAf*phiAdd);

    // Remove any patch-normal components of velocity
    U_ -= nHat()*(nHat() & U_);

    U_.correctBoundaryConditions();

    // Continuity check
    continuityCheck();
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::solveThickness()" << endl;
        tabSubtract();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicSingleLayer::kinematicSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const bool readFields
)
:
    surfaceFilmModel(modelType, mesh, g),

    momentumPredictor_(solution().subDict("PISO").lookup("momentumPredictor")),
    nOuterCorr_(readLabel(solution().subDict("PISO").lookup("nOuterCorr"))),
    nCorr_(readLabel(solution().subDict("PISO").lookup("nCorr"))),
    nNonOrthCorr_(readLabel(solution().subDict("PISO").lookup("nNonOrthCorr"))),

    cumulativeContErr_(0.0),

    Cf_(readScalar(coeffs().lookup("Cf"))),

    rho_
    (
        IOobject
        (
            "rhof",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimDensity, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    mu_
    (
        IOobject
        (
            "muf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure*dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigma_
    (
        IOobject
        (
            "sigmaf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/sqr(dimTime), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    delta_
    (
        IOobject
        (
            "deltaf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    U_
    (
        IOobject
        (
            "Uf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Us_
    (
        IOobject
        (
            "Usf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_,
        zeroGradientFvPatchScalarField::typeName
    ),
    Uw_
    (
        IOobject
        (
            "Uwf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_,
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaRho_
    (
        IOobject
        (
            delta_.name() + "*" + rho_.name(),
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", delta_.dimensions()*rho_.dimensions(), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    phi_
    (
        IOobject
        (
            "phi",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimLength*dimMass/dimTime
    ),

    primaryMassTrans_
    (
        IOobject
        (
            "primaryMassTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    cloudMassTrans_
    (
        IOobject
        (
            "cloudMassTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    cloudDiameterTrans_
    (
        IOobject
        (
            "cloudDiameterTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimLength, -1.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    USp_
    (
        IOobject
        (
            "USpf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector
        (
            "zero", dimMass*dimVelocity/dimArea/dimTime, vector::zero
        ),
        this->mappedPushedFieldPatchTypes<vector>()
    ),
    pSp_
    (
        IOobject
        (
            "pSpf",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),
    rhoSp_
    (
        IOobject
        (
            "rhoSpf",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),
    massImp_
    (
        IOobject
        (
            "massImp",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    USpPrimary_
    (
        IOobject
        (
            USp_.name(), // must have same name as USp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedVector("zero", USp_.dimensions(), vector::zero)
    ),
    pSpPrimary_
    (
        IOobject
        (
            pSp_.name(), // must have same name as pSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar("zero", pSp_.dimensions(), 0.0)
    ),
    rhoSpPrimary_
    (
        IOobject
        (
            rhoSp_.name(), // must have same name as rhoSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar("zero", rhoSp_.dimensions(), 0.0)
    ),

    UPrimary_
    (
        IOobject
        (
            "U", // must have same name as U to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector("zero", dimVelocity, vector::zero),
        this->mappedFieldAndInternalPatchTypes<vector>()
    ),
    pPrimary_
    (
        IOobject
        (
            "p", // must have same name as p to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    rhoPrimary_
    (
        IOobject
        (
            "rho", // must have same name as rho to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimDensity, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    muPrimary_
    (
        IOobject
        (
            "mu", // must have same name as mu to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimPressure*dimTime, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    availableMass_(regionMesh().nCells(), 0.0),

    injection_(*this, coeffs_),

    addedMassTotal_(0.0)
{
    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctThermoFields();

        //slottedCircle();

        deltaRho_ == delta_*rho_;
        phi_ = fvc::interpolate(deltaRho_*U_) & regionMesh().Sf();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kinematicSingleLayer::~kinematicSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void kinematicSingleLayer::addSources
(
    const label patchI,
    const label faceI,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::addSources()" << endl;
        Info<<tab.c_str()<< "\nSurface film: " << type() << ": adding to film source:" << nl
            << "    mass     = " << massSource << nl
            << "    momentum = " << momentumSource << nl
            << "    pressure = " << pressureSource << endl;
    }

    rhoSpPrimary_.boundaryField()[patchI][faceI] -= massSource;
    USpPrimary_.boundaryField()[patchI][faceI] -= momentumSource;
    pSpPrimary_.boundaryField()[patchI][faceI] -= pressureSource;

    addedMassTotal_ += massSource;

    if (debug)
    {
            Info<<tab.c_str()<< "leaving kinematicSingleLayer::addSources()" << endl;
        tabSubtract();
    }
}


void kinematicSingleLayer::preEvolveRegion()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::preEvolveRegion()" << endl;
    }

    transferPrimaryRegionThermoFields();

    correctThermoFields();

    transferPrimaryRegionSourceFields();

    // Reset transfer fields
    availableMass_ = mass();
//    availableMass_ = netMass(); //kvm, using netMass yields very low vaporization rates
    cloudMassTrans_ == dimensionedScalar("zero", dimMass, 0.0);
    cloudDiameterTrans_ == dimensionedScalar("zero", dimLength, -1.0);
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::preEvolveRegion()" << endl;
        tabSubtract();
    }
}


void kinematicSingleLayer::evolveRegion()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::evolveRegion()" << endl;
    }

    //kvm, if updateSubmodels occurs after solveContinuity, then vaporization and separation don't get called until after deltaRho_ has been calculated, resulting in errors
    updateSubmodels();//kvm, this is where vaporization and separation mass are accounted

    // Solve continuity for deltaRho_
    solveContinuity();

    // Implicit pressure source coefficient - constant
    tmp<volScalarField> tpp(this->pp());

    for (int oCorr=0; oCorr<nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution - varies with delta_
        tmp<volScalarField> tpu(this->pu());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update film wall and surface velocities
    updateSurfaceVelocities();

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::evolveRegion()" << endl;
        tabSubtract();
    }
}


scalar kinematicSingleLayer::CourantNumber() const
{
    if (debug)
    {
        //tabAdd();
        Info<<tab.c_str()<< "kinematicSingleLayer::CourantNumber()" << endl;
    }
    scalar CoNum = 0.0;

    if (regionMesh().nInternalFaces() > 0)
    {
        const scalar deltaT = time_.deltaTValue();

        const surfaceScalarField SfUfbyDelta =
            regionMesh().surfaceInterpolation::deltaCoeffs()*mag(phi_);
        const surfaceScalarField rhoDelta = fvc::interpolate(rho_*delta_);
        const surfaceScalarField& magSf = regionMesh().magSf();

        forAll(rhoDelta, i)
        {
            if (rhoDelta[i] > ROOTVSMALL)
            {
                CoNum = max(CoNum, SfUfbyDelta[i]/rhoDelta[i]/magSf[i]*deltaT);
            }
        }
    }

    reduce(CoNum, maxOp<scalar>());

    Info<< "Film max Courant number: " << CoNum << endl;

    return CoNum;
    if (debug)
    {
        Info<<tab.c_str()<< "leaving kinematicSingleLayer::CourantNumber()" << endl;
        //tabSubtract();
    }
}


const volVectorField& kinematicSingleLayer::U() const
{
    return U_;
}


const volVectorField& kinematicSingleLayer::Us() const
{
    return Us_;
}


const volVectorField& kinematicSingleLayer::Uw() const
{
    return Uw_;
}


const surfaceScalarField& kinematicSingleLayer::phi() const
{
    return phi_;
}


const volScalarField& kinematicSingleLayer::rho() const
{
    return rho_;
}


const volScalarField& kinematicSingleLayer::T() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::T() const"
    )   << "T field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Ts() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::Ts() const"
    )   << "Ts field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Tw() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::Tw() const"
    )   << "Tw field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Cp() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::Cp() const"
    )   << "Cp field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::kappa() const
{
    FatalErrorIn
    (
        "const volScalarField& kinematicSingleLayer::kappa() const"
    )   << "kappa field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


tmp<volScalarField> kinematicSingleLayer::primaryMassTrans() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "kinematicSingleLayer::primaryMassTrans",
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
}


const volScalarField& kinematicSingleLayer::cloudMassTrans() const
{
    return cloudMassTrans_;
}


const volScalarField& kinematicSingleLayer::cloudDiameterTrans() const
{
    return cloudDiameterTrans_;
}


void kinematicSingleLayer::info() const
{
    Info<< "\nSurface film: " << type() << endl;

    Info<< indent << "added mass         = "
        << returnReduce<scalar>(addedMassTotal_, sumOp<scalar>()) << nl
        << indent << "current mass       = "
        << gSum((deltaRho_*magSf())()) << nl
        << indent << "min/max(mag(U))    = " << min(mag(U_)).value() << ", "
        << max(mag(U_)).value() << nl
        << indent << "min/max(delta)     = " << min(delta_).value() << ", "
        << max(delta_).value() << nl;

    injection_.info(Info);
}

tmp<DimensionedField<scalar, volMesh> > kinematicSingleLayer::qWall() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::qWall()",
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

tmp<DimensionedField<scalar, volMesh> > kinematicSingleLayer::qRad() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::qRad()",
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

const volScalarField& kinematicSingleLayer::omega() const
{
    FatalErrorIn("const volScalarField& kinematicSingleLayer::omega() const")
        << "omega field not available for " << type() << abort(FatalError);

    return reinterpret_cast<const volScalarField&>(null);
}

tmp<DimensionedField<scalar, volMesh> >  kinematicSingleLayer::Srho() const
{
    return tmp<DimensionedField<scalar, volMesh> > 
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::Srho",
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
}


tmp<DimensionedField<scalar, volMesh> > kinematicSingleLayer::Srho
(
    const label i
) const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::Srho(i)",
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
}


tmp<DimensionedField<scalar, volMesh> > kinematicSingleLayer::Sh() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "kinematicSingleLayer::Sh",
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
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void kinematicSingleLayer::slottedCircle(){
    Info << "Slotted Circle Initialization\n";

    static const vector center(0,0,0);
    static const scalar circleRadius=.01;
    static const scalar slotWidth=circleRadius*.0;
    static const volVectorField& cellCentres = regionMesh().C();

    forAll(delta_,i){
        scalar x=cellCentres[i][0];
        scalar y=cellCentres[i][1];
        scalar radius=sqrt(pow(x-center.x(),2)+pow(y-center.y(),2));
        if(radius<circleRadius){
            delta_[i]=0.0003;
        }
        else{
            delta_[i]=0.0;
        }
        /*cut out slot*/
        if(y<0&&x>-slotWidth&&x<slotWidth){
            delta_[i]=0.0;
        }
    }
    return;
}

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
