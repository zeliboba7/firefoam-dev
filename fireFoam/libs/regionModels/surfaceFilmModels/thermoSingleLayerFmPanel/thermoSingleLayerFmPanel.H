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

Class
    thermoSingleLayerFmPanel

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15 K

SourceFiles
    thermoSingleLayerFmPanel.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayerFmPanel_H
#define thermoSingleLayerFmPanel_H

#include "thermoSingleLayerPw.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// Forward declaration of classes
class heatTransferModel;
class phaseChangeModel;
class filmRadiationModel;

/*---------------------------------------------------------------------------*\
                      Class thermoSingleLayerFmPanel Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayerFmPanel
:
    public thermoSingleLayerPw
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayerFmPanel(const thermoSingleLayerFmPanel&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayerFmPanel&);


protected:

    // Protected data

                Switch augmentedRadiation_;
                Switch xiangyang_;
                Switch perfectlyWettedInlet_;
                string perfectlyWettedInletName_;
                dimensionedScalar perfectlyWettedInletDistance_;
        // Thermo properties

                //- Radiation Flux / [W/m^2]
                volScalarField qRad_;
                //- Film to Wall Flux / [W/m^2]
                volScalarField qFilmToWall_;

                //- Constant radiation source [W/m^2]
                dimensionedScalar qRadConstant_;
                //- radiation source applied to window xmin,ymin xmax,ymax
                dimensionedScalar qRadXMax_;
                dimensionedScalar qRadXMin_;
                dimensionedScalar qRadYMax_;
                dimensionedScalar qRadYMin_;
                dimensionedScalar qRadBegin_;
                dimensionedScalar qRadEnd_;
                dimensionedScalar qRadEmissivity_;
                dimensionedScalar qRadAbsorptivity_;

                Switch solveLumpedCapacitance_;

                class lumpedCapacitance
                {
                    public:
                    dimensionedScalar thickness;
                    dimensionedScalar density;
                    dimensionedScalar cp;
                    dimensionedScalar k;
                    dimensionedScalar Tinit;
                    //lumpedCapacitance(int,int,int,int,int)
                    lumpedCapacitance()
                        :
                            thickness(1.2e-3),
                            density(2700.),
                            cp(875.),
                            k(177.),
                            Tinit(300.)
                    {}
                };
                lumpedCapacitance lc;

    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();

        virtual void updateSubmodels();

        //- Set qRad
        virtual void updateQRad();

        //- compute/output wettedAreaInfo  
        virtual void wettedAreaInfo() const;

        //- Return the wall/surface heat transfer term for the enthalpy equation
        virtual tmp<fvScalarMatrix> q(volScalarField& h) const;
        
        void zeroContactAngleInlet();

public:

    //- Runtime type information
    TypeName("thermoSingleLayerFmPanel");


    // Constructors

        //- Construct from components
        thermoSingleLayerFmPanel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    //- Destructor
    virtual ~thermoSingleLayerFmPanel();

    // Member Functions

         // Source fields (read/write access)

            //- Post-evolve region
            virtual void postEvolveRegion();
            
            //qWall virtual const volScalarField& qWall() const;
            virtual tmp<DimensionedField<scalar, volMesh> > qWall() const;

            //qrad virtual const volScalarField& qRad() const;
            virtual tmp<DimensionedField<scalar, volMesh> > qRad() const;
            
            //- Correct the film surface and wall temperatures
            virtual void updateSurfaceTemperatures();

            //- Compute splashed mass for Xiangyang's experimental setup
            void integrateSplashMass();

       // I-O

            //- Provide some feedback
            virtual void info() const;
};


// * * * * *//notNeeded?   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end nam//notNeeded?  espace surfaceFilmModels
} // end nam//notNeeded?  espace regionModels
} // end nam//notNeeded?  espace Foam


// * * * * *//notNeeded?   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
