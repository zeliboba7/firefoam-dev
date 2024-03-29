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
    Foam::surfaceFilmModel

Description
    Base class for surface film models

SourceFiles
    surfaceFilmModelI.H
    surfaceFilmModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilmModel_H
#define surfaceFilmModel_H

#include "singleLayerRegion.H"

#include "dimensionedVector.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "labelList.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFilmModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilmModel
:
    public singleLayerRegion
{
public:

    // Data types

        //- Enumeration listing the possible thermo types
        enum thermoModelType
        {
            tmConstant,
            tmSingleComponent
        };

        //- Named enumeration for the thermoType
        static const NamedEnum<thermoModelType, 2> thermoModelTypeNames_;


private:

    // Private Member Functions

        //- Construct fields
        void constructMeshObjects();

        //- Disallow default bitwise copy construct
        surfaceFilmModel(const surfaceFilmModel&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceFilmModel&);


protected:

    // Protected data

        //- Flag to indicate whether pyrolysis region coupled to a film region
        Switch pyrCoupled_;
        //
        //- Flag to indicate whether diagnostics are turned on
        Switch diagnostics_;
        
        //- Pointer to film temperature field
        autoPtr<volScalarField> pyrTemperaturePtr_;

        //- Acceleration due to gravity [m/s2]
        const dimensionedVector& g_;

        //- Thermo type
        thermoModelType thermoModel_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("surfaceFilmModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             surfaceFilmModel,
             mesh,
             (
                const word& modelType,
                const fvMesh& mesh,
                const dimensionedVector& g
             ),
             (modelType, mesh, g)
         );

    // Constructors

        //- Construct null
        surfaceFilmModel(const fvMesh& mesh);

        //- Construct from type name, mesh and gravity vector
        surfaceFilmModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    // Selectors

        //- Return a reference to the selected surface film model
        static autoPtr<surfaceFilmModel> New
        (
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    //- Destructor
    virtual ~surfaceFilmModel();


    // Member Functions

        // Access

            //- Return the accleration due to gravity
            inline const dimensionedVector& g() const;

            //- Return the thermo type
            inline const thermoModelType& thermoModel() const;

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchI,
                const label faceI,
                const scalar massSource,
                const vector& momentumSource,
                const scalar pressureSource,
                const scalar energySource
            ) = 0;

            
        // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

        // Solution parameters

            //- Courant number evaluation
            virtual scalar CourantNumber() const;


        // Fields

            //- Return the film thickness [m]
            virtual const volScalarField& delta() const = 0;

            //- Return the film velocity [m/s]
            virtual const volVectorField& U() const = 0;

            //- Return the film surface velocity [m/s]
            virtual const volVectorField& Us() const = 0;

            //- Return the film wall velocity [m/s]
            virtual const volVectorField& Uw() const = 0;

            //- Return the film density [kg/m3]
            virtual const volScalarField& rho() const = 0;

            //- Return the film mean temperature [K]
            virtual const volScalarField& T() const = 0;

            //- Return the film surface temperature [K]
            virtual const volScalarField& Ts() const = 0;

            //- Return the film wall temperature [K]
            virtual const volScalarField& Tw() const = 0;

            //- Return the film specific heat capacity [J/kg/K]
            virtual const volScalarField& Cp() const = 0;

            //- Return the film thermal conductivity [W/m/K]
            virtual const volScalarField& kappa() const = 0;

            //- Return the film surface tension [N/m]
            virtual const volScalarField& sigma() const = 0;

            //- Return convective heat flux from film to wall [W/m2]
            virtual tmp<DimensionedField<scalar, volMesh> > qWall() const;

            //- Return radiative heat flux from gas to film [W/m2]
            virtual tmp<DimensionedField<scalar, volMesh> > qRad() const;

           //- Return the film wettness (0=dry, 1=wet)
           virtual const volScalarField& omega() const = 0;

            // Transfer fields - to the primary region

                //- Return mass transfer source - Eulerian phase only
                virtual tmp<volScalarField> primaryMassTrans() const = 0;

                //- Return the film mass available for transfer
                virtual const volScalarField& cloudMassTrans() const = 0;

                //- Return the parcel diameters originating from film
                virtual const volScalarField& cloudDiameterTrans() const = 0;


        // Source fields

            // Mapped into primary region

                //- Return total mass source - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Srho() const;

                //- Return mass source for specie i - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Srho
                (
                    const label i
                ) const;

                //- Return enthalpy source - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Sh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceFilmModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
