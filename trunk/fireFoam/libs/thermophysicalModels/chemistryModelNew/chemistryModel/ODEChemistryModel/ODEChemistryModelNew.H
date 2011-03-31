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
    Foam::ODEChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    ODEChemistryModelNewI.H
    ODEChemistryModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ODEChemistryModelNew_H
#define ODEChemistryModelNew_H

#include "Reaction.H"
#include "ODE.H"
#include "volFieldsFwd.H"

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;


/*---------------------------------------------------------------------------*\
                      Class ODEChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class ODEChemistryModelNew
:
    public CompType,
    public ODE
{
    // Private Member Functions

        //- Disallow copy constructor
        ODEChemistryModelNew(const ODEChemistryModelNew&);

        //- Disallow default bitwise assignment
        ODEChemistryModelNew operator=(const ODEChemistryModelNew&);


protected:

    // Private data

        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reactions
        const PtrList<Reaction<ThermoType> >& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- Chemical source term [kg/m3/s]
        PtrList<scalarField> RR_;


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<scalarField>& RR();


public:

    //- Runtime type information
    TypeName("ODEChemistryModelNew");


    // Constructors

        //- Construct from components
        ODEChemistryModelNew
        (
            const fvMesh& mesh,
            const word& ODEmodelName,
            const word& thermoTypeName
        );



    //- Destructor
    virtual ~ODEChemistryModelNew();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;

        //- The number of reactions
        inline label nReaction() const;

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual tmp<scalarField> omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        //- Return the reaction rate for iReaction and the reference
        //  species and charateristic times
        virtual scalar omegaI
        (
            label iReaction,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        //- Calculates the reaction rates
        virtual void calculate();


        //- Update concentrations in reaction i given dt and reaction rate omeg
        // used by sequential solver
        void updateConcsInReactionI
        (
            const label i,
            const scalar dt,
            const scalar omeg,
            scalarField& c
        ) const;


        //- Update matrix RR for reaction i. Used by EulerImplicit
        void updateRRInReactionI
        (
            const label i,
            const scalar pr,
            const scalar pf,
            const scalar corr,
            const label lRef,
            const label rRef,
            simpleMatrix<scalar>& RR
        ) const;

        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms
            inline tmp<volScalarField> RR(const label i) const;

            //- Solve the reaction system for the given start time and time
            //  step and return the characteristic time
            virtual scalar solve(const scalar t0, const scalar deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            //- Return the heat release, i.e. enthalpy/sec [m2/s3]
            virtual tmp<volScalarField> dQ() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            virtual scalar solve
            (
                scalarField &c,
                const scalar T,
                const scalar p,
                const scalar t0,
                const scalar dt
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ODEChemistryModelNewI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ODEChemistryModelNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
