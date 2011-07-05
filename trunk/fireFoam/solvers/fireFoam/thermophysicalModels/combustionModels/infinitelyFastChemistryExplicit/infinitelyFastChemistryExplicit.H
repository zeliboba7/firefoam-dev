/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::combustionModels::infinitelyFastChemistryExplicit

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    infinitelyFastChemistryExplicit.C

\*---------------------------------------------------------------------------*/

#ifndef infinitelyFastChemistryExplicit_H
#define infinitelyFastChemistryExplicit_H

#include "combustionModel.H"
#include "singleStepReactingMixture.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class infinitelyFastChemistryExplicit Declaration
\*---------------------------------------------------------------------------*/

class infinitelyFastChemistryExplicit
:
    public combustionModel
{
    // Private data

        //- Model constant
        scalar C_;

        //- Reference to singleStepReactingMixture mixture
        singleStepReactingMixture<gasThermoPhysics>& singleMixture_;

        //- Normalised consumption rate of (fu - fres)
        volScalarField wFuelNorm_;


    // Private Member Functions

        //- Disallow copy construct
        infinitelyFastChemistryExplicit(const infinitelyFastChemistryExplicit&);

        //- Disallow default bitwise assignment
        void operator=(const infinitelyFastChemistryExplicit&);


public:

    //- Runtime type information
    TypeName("infinitelyFastChemistryExplicit");


    // Constructors

        //- Construct from components
        infinitelyFastChemistryExplicit
        (
            const dictionary& combustionProps,
            hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    //- Destructor
    virtual ~infinitelyFastChemistryExplicit();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix, i.e. source term for fuel equation
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //- Return normalised consumption rate of (fu - fres)
            virtual tmp<volScalarField> wFuelNorm() const;


    // I-O

            //- Update properties from given dictionary
            virtual bool read(const dictionary& combustionProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
