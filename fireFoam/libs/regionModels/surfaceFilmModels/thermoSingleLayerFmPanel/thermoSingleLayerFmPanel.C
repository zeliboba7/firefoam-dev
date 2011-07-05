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

#include "thermoSingleLayerFmPanel.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "directMappedWallPolyPatch.H"
#include "mathematicalConstants.H" 
#include "radiationConstants.H"

// Sub-models
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"
#include "stdio.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace regionModels
    {
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(thermoSingleLayerFmPanel, 0);
        addToRunTimeSelectionTable(surfaceFilmModel, thermoSingleLayerFmPanel, mesh);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

/*called from surfaceFilmModel::evolve()*/
bool thermoSingleLayerFmPanel::read()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayerFmPanel::read()" << endl;
    }
    // no additional properties to read
    if(thermoSingleLayerPw::read())
    {
        augmentedRadiation_ = coeffs_.lookupOrDefault<Switch>("augmentedRadiation",false);
        if(augmentedRadiation_){
            const dictionary& subdict=coeffs_.subDict("augmentedRadiationCoeffs");
            subdict.lookup("qRadConstant") >> qRadConstant_;
            subdict.lookup("qRadXMax") >> qRadXMax_; 
            subdict.lookup("qRadXMin") >> qRadXMin_; 
            subdict.lookup("qRadYMax") >> qRadYMax_; 
            subdict.lookup("qRadYMin") >> qRadYMin_; 
            subdict.lookup("qRadBegin") >> qRadBegin_; 
            subdict.lookup("qRadEnd") >> qRadEnd_; 
            subdict.lookup("qRadEmissivity") >> qRadEmissivity_; 
            subdict.lookup("qRadAbsorptivity") >> qRadAbsorptivity_; 
        }

        solveLumpedCapacitance_ = coeffs_.lookupOrDefault<Switch>("solveLumpedCapacitance",false);
        if(solveLumpedCapacitance_){
            const dictionary& subdict=coeffs_.subDict("solveLumpedCapacitanceCoeffs");
            subdict.lookup("thickness") >> (lc.thickness); //thickness of aluminum panel
            subdict.lookup("density") >> lc.density;
            subdict.lookup("cp") >> lc.cp;
            subdict.lookup("k") >> lc.k;
            subdict.lookup("Tinit") >> lc.Tinit;
        }
        
        xiangyang_ = coeffs_.lookupOrDefault<Switch>("XiangYang",false);
        perfectlyWettedInlet_ = coeffs_.lookupOrDefault<Switch>("perfectlyWettedInlet",false);
        if(perfectlyWettedInlet_){
            const dictionary& subdict=coeffs_.subDict("perfectlyWettedInletCoeffs");
            perfectlyWettedInletName_ = subdict.lookup("inletName");
            subdict.lookup("offsetDistance") >> perfectlyWettedInletDistance_ ;
        }
        if (debug)
        {
            Info<<tab.c_str()<< "leaving thermoSingleLayerFmPanel::read()" << endl;
            tabSubtract();
        }
        return true;
    }
    else
    {
        return false;
    }
}

void thermoSingleLayerFmPanel::updateSubmodels()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayerPw::updateSubmodels()" << endl;
    }

    thermoSingleLayerPw::updateSubmodels();

    /*partially wetted treatment*/
    if(perfectlyWettedInlet_){
        zeroContactAngleInlet();
    }
    thermoSingleLayerPw::updateContactLine();

    if(augmentedRadiation_){
        updateQRad();
    }

    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayerPw::updateSubmodels()" << endl;
        tabSubtract();
    }
}


void thermoSingleLayerFmPanel::updateSurfaceTemperatures()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayerFmPanel::updateSurfaceTemperatures2()" << endl;
    }

    if(!solveLumpedCapacitance_){
        thermoSingleLayerPw::updateSurfaceTemperatures();
    }
    else{
        // Push boundary film temperature values into internal field
        static label first=1;
        if(first){
            for (label i=0; i<intCoupledPatchIDs_.size(); i++)
            {
                label patchI = intCoupledPatchIDs_[i];
                const polyPatch& pp = regionMesh().boundaryMesh()[patchI];
                UIndirectList<scalar>(Tw_, pp.faceCells()) =
                    T_.boundaryField()[patchI];
            }
            first=0;
        }

        /*lumped-capacitance bc model for aluminum plate*/
        if(solveLumpedCapacitance_){
            const volScalarField qDotFilm = htcw_->h()*(T_ - Tw_);
            const dimensionedScalar dt = time().deltaT();

            #include "solveSolid.H"
        }

        qFilmToWall_ = htcw_->h()*(T_ - Tw_);
        //Info << qFilmToWall_<<endl;
        

        // Update film surface temperature
        Ts_ = T_;
        Ts_.correctBoundaryConditions();
    }
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayerFmPanel::updateSurfaceTemperatures2()" << endl;
        tabSubtract();
    }
}

tmp<fvScalarMatrix> thermoSingleLayerFmPanel::q
(
 volScalarField& hs
 ) const
{
    dimensionedScalar Tstd("Tstd", dimTemperature, 298.15);
    return
        (
         - fvm::Sp(htcs_->h()/Cp_, hs) - htcs_->h()*(Tstd - TPrimary_)
         - fvm::Sp(htcw_->h()/Cp_, hs) - htcw_->h()*(Tstd - Tw_)
         /* When the film becomes dry then qRad_ is to strong unless 
          * heat loss (conduction through solid, surface emission, etc) 
          * is accounted for */
      //+ omega_*qRad_ 
      //TODO:  I am changing qRad_ here, but in standardPhaseChange qRad_ is unmodified
      + qRadAbsorptivity_*qRad_ 
      //this term can case temp to be wildly unstable!!!     - qRadEmissivity_*sigmaSB_*pow(T_,4) //kvm, need to lineraize this term and use hs/Cp_ instead of T_
    );
}

void thermoSingleLayerFmPanel::zeroContactAngleInlet()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayerFmPanel::zeroContactAngleInlet()" << endl;
    }
    /*force contact angle to zero near inlet*/
    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();
    const volVectorField& cellCentres = regionMesh().C();
    static bool first=true;
    if(!contactAngleFromFile_){
        if(first){
            forAll(bm,patchI){

                if(bm[patchI].name()== perfectlyWettedInletName_){
                    const pointField& bCenters = bm[patchI].faceCentres();
                    forAll(cellCentres,i){
                        scalar distance=1e10;
                        forAll(bCenters,j){
                            scalar tmpDistance=sqrt(
                               +pow(bCenters[j][0]-cellCentres[i][0],2)
                               +pow(bCenters[j][1]-cellCentres[i][1],2)
                               +pow(bCenters[j][2]-cellCentres[i][2],2)
                               );
                            distance=(tmpDistance<distance)?(tmpDistance):(distance);
                        }
                        scalar distance1=perfectlyWettedInletDistance_.value();
                        scalar distance2=2.0*distance1;
                        if(distance<distance1){
                            contactAngle_[i]=0.0;
                        }
                        else if(distance<distance2){
                            contactAngle_[i]*=(distance2-distance)/(distance2-distance1);
                        }
                    }
                }
            }
            first=false;
        }
    }

    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayerFmPanel::zeroContactAngleInlet()" << endl;
        tabSubtract();
    }
}

void thermoSingleLayerFmPanel::updateQRad()
{
    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayerFmPanel::updateQRad()" << endl;
    }

    const volVectorField& cellCentres = regionMesh().C();
    forAll(qRad_,index){
        if(time().time().value() >= qRadBegin_.value() && time().time().value() <= qRadEnd_.value()){
            if(cellCentres[index][1]>qRadYMin_.value()&&cellCentres[index][1]<qRadYMax_.value()){
                if(cellCentres[index][0]>qRadXMin_.value()&&cellCentres[index][0]<qRadXMax_.value()){
                    //assume uniform heat flux for now
                    //qRad_[index]=qRadConstant_.value()-5.67e-8*pow((T_[index]-300.0),4);
                    qRad_[index]=qRadConstant_.value();
                }
            }
        }
        else{
            qRad_[index]=0.0;
        }

    }

    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayerFmPanel::updateQRad()" << endl;
        tabSubtract();
    }
}

void thermoSingleLayerFmPanel::wettedAreaInfo() const
{
    scalar wetArea=0.0;
    scalar dryArea=0.0;
    scalar totalArea=0.0;
    forAll(omega_,index){
        totalArea+=magSf()[index];
        if(omega_[index]==1){
            wetArea+=magSf()[index];
        }
        else if(omega_[index]==0){
            dryArea+=magSf()[index];
        }
    }

    const scalar totalWetArea=returnReduce(wetArea, sumOp<scalar>());
    const scalar totalDryArea=returnReduce(dryArea, sumOp<scalar>());
    Info << indent << "wettedAreaFraction  = " << totalWetArea/(totalDryArea+totalWetArea+ROOTVSMALL) << nl;
    Info << indent << "total area          = " << (totalDryArea+totalWetArea) << " m2"<< nl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayerFmPanel::thermoSingleLayerFmPanel
(
 const word& modelType,
 const fvMesh& mesh,
 const dimensionedVector& g
 )
:
    thermoSingleLayerPw(modelType, mesh, g),
    augmentedRadiation_(false),
    xiangyang_(false),
    perfectlyWettedInlet_(false),
    perfectlyWettedInletName_(""),
    perfectlyWettedInletDistance_(0.0),
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
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0)
    ),
    qRadConstant_(0.0),
    qRadXMax_(0.0),
    qRadXMin_(0.0),
    qRadYMax_(0.0),
    qRadYMin_(0.0),
    qRadBegin_(0.0),
    qRadEnd_(0.0),
    qRadEmissivity_(1.0),
    qRadAbsorptivity_(1.0),
    solveLumpedCapacitance_(false),
    lc()
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayerFmPanel::~thermoSingleLayerFmPanel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayerFmPanel::postEvolveRegion()
{

    if (debug)
    {
        tabAdd();
        Info<<tab.c_str()<< "thermoSingleLayerFmPanel::postEvolveRegion()" << endl;
    }

    if (time().outputTime())
    {
        //massPhaseChangeForPrimary_.write();
        htcw_->h()->write();
        htcs_->h()->write();
        Tw_.write();
        qRad_.write();
        qWall()->write();
#include "manualDecomposition.H"
    }
    if(diagnostics_){
        static dimensionedScalar qConvGas("qConvGas", dimEnergy, 0.0);
        qConvGas-=sum(omega_*htcs_->h()*(T_ - TPrimary_)*magSf()*time_.deltaT());   
        INFO << "qConvGas " << time_.value() << " " << qConvGas.value() << endl;
        static dimensionedScalar qPhaseChange("qPhaseChange", dimEnergy, 0.0);
        qPhaseChange-=sum(primaryEnergyPCTrans_);   
        INFO << "qPhaseChange " << time_.value() << " " << qPhaseChange.value() << endl;
        static dimensionedScalar qSensible("qSensible", dimEnergy, 0.0);
        //qSensible-=sum(delta_*rho_*(hs_-hs_.oldTime())*magSf());   
        qSensible-=sum((delta_*rho_*hs_-delta_.oldTime()*rho_.oldTime()*hs_.oldTime())*magSf());   
        INFO << "qSensible " << time_.value() << " " << qSensible.value() << endl;
        static dimensionedScalar qImp("qImp", dimEnergy, 0.0);
        qImp-=sum(rhoSp_*hs_*magSf()*time_.deltaT());   
        INFO << "qImp " << time_.value() << " " << qImp.value() << endl;
        static dimensionedScalar qImpSens("qImpSens", dimEnergy, 0.0);
        qImpSens-=sum(hsSp_*magSf()*time_.deltaT());   
        INFO << "qImpSens " << time_.value() << " " << qImpSens.value() << endl;
        INFO << "qTotal " << time_.value() << " " << qConvGas.value()+qPhaseChange.value()+qSensible.value()+qImpSens.value()+qImp.value() << endl;
//        static dimensionedScalar qAdv("qAdv", dimEnergy, 0.0);
//        qAdv-=sum(fvc::surfaceSum(phi_*fvc::interpolate(hs_)));   
//        INFO << "qAdv " << time_.value() << " " << qAdv.value() << endl;
    }
    if (xiangyang_){
        integrateSplashMass();
    }
    if (debug)
    {
        Info<<tab.c_str()<< "leaving thermoSingleLayerFmPanel::postEvolveRegion()" << endl;
        tabSubtract();
    }
}

tmp<DimensionedField<scalar, volMesh> >
thermoSingleLayerFmPanel::qWall() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
     qFilmToWall_
    );
}

tmp<DimensionedField<scalar, volMesh> >
thermoSingleLayerFmPanel::qRad() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        qRad_
    );
}

void thermoSingleLayerFmPanel::integrateSplashMass(){
    static const vector center(0,0,0);
    static scalarList binRadius(6,0.0);
    binRadius[0]=0.018;
    binRadius[1]=0.040;
    binRadius[2]=0.060;
    binRadius[3]=0.089;
    binRadius[4]=0.134;
    binRadius[5]=0.190;
    //const scalarList binRadius(.0,.018,.040,.060,.089,.134,.190);
    static const volVectorField& cellCentres = regionMesh().C();
    static scalarList cummulatedMass(binRadius.size(),0.0);

    forAll(delta_,i){
        scalar x=cellCentres[i][0];
        scalar y=cellCentres[i][1];
        scalar radius=sqrt(pow(x-center.x(),2)+pow(y-center.y(),2));
        for(label j=0;j<binRadius.size();j++){
            if(j==0){
                if(radius<binRadius[j]){
                    cummulatedMass[j]-=rhoSp_[i]*magSf()[i]*time().deltaT().value();
                }
            }
            else{
                if(binRadius[j-1]<=radius&&radius<binRadius[j]){
                    cummulatedMass[j]-=rhoSp_[i]*magSf()[i]*time().deltaT().value();
                }
            }

        }
    }
    static char buffer[256];
    sprintf(buffer,"SplashedMass %10.5f ",time().value());
    Info << buffer;
    forAll(cummulatedMass,j){
        sprintf(buffer," % 10.5g ",cummulatedMass[j]);
        Info << buffer;
    }
    sprintf(buffer," % 10.5g ",sum(cummulatedMass));
    Info << buffer;
    Info << endl;

}

void thermoSingleLayerFmPanel::info() const
{
    thermoSingleLayerPw::info();

    wettedAreaInfo();

    if(solveLumpedCapacitance_){
        Info<< indent << "min/max(Tw)         = " << min(Tw_).value() << ", " << max(Tw_).value() << nl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace surfaceFilmModels
} // end namespace regionModels
} // end namespace Foam


// ************************************************************************* //
