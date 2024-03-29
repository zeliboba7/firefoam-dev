/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::basicMultiComponentMixtureNew

Description
    Abstract base class for multi-component mixtures.

SourceFiles
    basicMultiComponentMixtureNew.C

\*---------------------------------------------------------------------------*/

#ifndef basicMultiComponentMixtureNew_H
#define basicMultiComponentMixtureNew_H

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basicMultiComponentMixtureNew Declaration
\*---------------------------------------------------------------------------*/

class basicMultiComponentMixtureNew
{

protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


public:

    // Constructors

        //- Construct from dictionary and mesh
        basicMultiComponentMixtureNew
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&
        );


    //- Destructor
    virtual ~basicMultiComponentMixtureNew()
    {}


    // Member functions

        //- Return the table of species
        const speciesTable& species() const
        {
            return species_;
        }

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        inline scalar fres(const scalar ft, const scalar stoicRatio) const;

        inline tmp<volScalarField> fres
        (
            const volScalarField& ft,
            const dimensionedScalar& stoicRatio
        ) const;


        // Per specie properties

            //- Number of moles []
            virtual scalar nMoles(const label specieI) const = 0;

            //- Molecular weight [kg/kmol]
            virtual scalar W(const label specieI) const = 0;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp(const label specieI, const scalar T) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv(const label specieI, const scalar T) const = 0;

            //- Enthalpy [J/kg]
            virtual scalar H(const label specieI, const scalar T) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs(const label specieI, const scalar T) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label specieI) const = 0;

            //- Entropy [J/(kg K)]
            virtual scalar S(const label specieI, const scalar T) const = 0;

            //- Internal energy [J/kg]
            virtual scalar E(const label specieI, const scalar T) const = 0;

            //- Gibbs free energy [J/kg]
            virtual scalar G(const label specieI, const scalar T) const = 0;

            //- Helmholtz free energy [J/kg]
            virtual scalar A(const label specieI, const scalar T) const = 0;


        // Per specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu(const label specieI, const scalar T) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa(const label specieI, const scalar T) const = 0;

            //- Thermal diffusivity [kg/m/s]
            virtual scalar alpha(const label specieI, const scalar T) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicMultiComponentMixtureNewI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
