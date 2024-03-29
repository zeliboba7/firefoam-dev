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
    Foam::ReactingParcel

Description
    Reacting parcel class with one/two-way coupling with the continuous
    phase.

SourceFiles
    ReactingParcelI.H
    ReactingParcel.C
    ReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingParcel_H
#define ReactingParcel_H

#include "ThermoParcel.H"
#include "ReactingCloud.H"
#include "reactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                        Class ReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingParcel
:
    public reactingParcel,
    public ThermoParcel<ParcelType>
{
public:

    //- Class to hold reacting particle constant properties
    class constantProperties
    :
        public ThermoParcel<ParcelType>::constantProperties
    {
        // Private data

            //- Minimum pressure [Pa]
            const scalar pMin_;

            //- Constant volume flag - e.g. during mass transfer
            Switch constantVolume_;

            //- Vaporisation temperature [K]
            const scalar Tvap_;

            //- Boiling point [K]
            const scalar Tbp_;


    public:

        //- Constructor
        constantProperties(const dictionary& parentDict);

        // Access

            //- Return const access to the minimum pressure
            inline scalar pMin() const;

            //- Return const access to the constant volume flag
            inline Switch constantVolume() const;

            //- Return const access to the vaporisation temperature
            inline scalar Tvap() const;

            //- Return const access to the boiling point
            inline scalar Tbp() const;
   };


    //- Class used to pass reacting tracking data to the trackToFace function
    class trackData
    :
        public ThermoParcel<ParcelType>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            ReactingCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

            //- Interpolator for continuous phase pressure field
            const interpolation<scalar>& pInterp_;


    public:

        typedef ReactingCloud<ParcelType> cloudType;


        // Constructors

            //- Construct from components
            inline trackData
            (
                ReactingCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolation<scalar>& rhoInterp,
                const interpolation<vector>& UInterp,
                const interpolation<scalar>& muInterp,
                const interpolation<scalar>& TInterp,
                const interpolation<scalar>& CpInterp,
                const interpolation<scalar>& pInterp,
                const vector& g
            );


        // Member functions

            //- Return access to the owner cloud
            inline ReactingCloud<ParcelType>& cloud();

            //- Return const access to the constant properties
            inline const constantProperties& constProps() const;

            //- Return const access to the interpolator for continuous
            //  phase pressure field
            inline const interpolation<scalar>& pInterp() const;
    };


protected:

    // Protected data

        // Parcel properties

            //- Initial particle mass [kg]
            scalar mass0_;

            //- Mass fractions of mixture []
            scalarField Y_;


        // Cell-based quantities

            //- Pressure [Pa]
            scalar pc_;


    // Protected member functions

        //- Calculate Phase change
        template<class TrackData>
        void calcPhaseChange
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar Ts,           // Surface temperature
            const scalar nus,          // Surface kinematic viscosity
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label idPhase,       // id of phase involved in phase change
            const scalar YPhase,       // total mass fraction
            const scalarField& YComponents, // component mass fractions
            scalarField& dMassPC,      // mass transfer - local to particle
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        );

        //- Update mass fraction
        scalar updateMassFraction
        (
            const scalar mass0,
            const scalarField& dMass,
            scalarField& Y
        ) const;


public:

    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Runtime type information
        TypeName("ReactingParcel");


    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline ReactingParcel
        (
            ReactingCloud<ParcelType>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        inline ReactingParcel
        (
            ReactingCloud<ParcelType>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const scalarField& Y0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ReactingParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ReactingParcel(const ReactingParcel& p);

        //- Construct and return a clone
        autoPtr<ReactingParcel> clone() const
        {
            return autoPtr<ReactingParcel>(new ReactingParcel(*this));
        }


    // Member Functions

        // Access

            //- Return const access to initial mass
            inline scalar mass0() const;

            //- Return const access to mass fractions of mixture
            inline const scalarField& Y() const;

            //- Return the owner cell pressure
            inline scalar pc() const;


        // Edit

            //- Return access to initial mass
            inline scalar& mass0();

            //- Return access to mass fractions of mixture
            inline scalarField& Y();


        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct surface values due to emitted species
            template<class TrackData>
            void correctSurfaceValues
            (
                TrackData& td,
                const label cellI,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappa
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );


        // I-O

            //- Read
            static void readFields(Cloud<ParcelType>& c);

            //- Write
            static void writeFields(const Cloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
