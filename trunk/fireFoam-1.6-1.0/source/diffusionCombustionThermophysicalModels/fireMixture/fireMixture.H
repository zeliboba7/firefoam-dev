/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original developer
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    This file is part of a preliminary version of the FireFOAM code which is 
    currently under development at FM Global. FM Global designates this file
    WITHOUT ANY WARRANTY to be used for development purpose only. 

License
    This file is based on OpenFOAM.

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
    fireMixture

Description

SourceFiles
    fireMixture.C

\*---------------------------------------------------------------------------*/

#ifndef fireMixture_H
#define fireMixture_H

#include "basicMultiComponentMixture.H"
//#include "interpolationLookUpTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fireMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class fireMixture
:
    public basicMultiComponentMixture
{
public:
    // Public data types
        //- Enumeration of pdf calculation methods
        enum pdfMethod
        {
            delta,          // 0- delta function
            lookUpTable,    // 1- lookup table
            recursive       // 2-recursive beta pdf integration
        };

private:
    // Private data

        static const int nSpecies_ = 6;
        static const char* specieNames_[6];

        dimensionedScalar stoicRatio_;

        ThermoType fuel_;
        ThermoType oxidant_;
        ThermoType products_;

        mutable ThermoType mixture_;

        //- Mixture fraction
        volScalarField& ft_;

        //- Regress variable
        volScalarField& b_;

        //- Mixture fraction
        volScalarField& ftVar_;

        //- Fuel mass fraction
        mutable volScalarField& fu_;

        //- Oxidizer mass fraction
        mutable volScalarField& ox_;

        //- Product mass fraction
        mutable volScalarField& pr_;

        //- Look Up table
        //interpolationLookUpTable<scalar> LookUpTable_;

        //- Construct as copy (not implemented)
        fireMixture(const fireMixture<ThermoType>&);

        //- Flag of how to do the pdf integration 
        label pdfMethod_;  

        //- Order of polynomial
        static const int nPoly_=21; 

        //- Polynomial coefficient
        scalar coeff_[nPoly_];


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        fireMixture(const dictionary&, const fvMesh&);


    // Member functions

        //interpolationLookUpTable<scalar> LookUpTable()
        //{
        //    return LookUpTable_;
        //}

        const dimensionedScalar& stoicRatio() const
        {
            return stoicRatio_;
        }

        const ThermoType& mixture
        (
            const scalar, 
            const scalar, 
            const scalar, 
            scalar&, 
            scalar&, 
            scalar&
        ) const;

        scalar getMassFraction
        (
            const scalar ,
            const scalar ,
            const scalar 
        ) const;

        const ThermoType& cellMixture(const label celli) const
        {
            return mixture
            (
                ft_[celli], 
                b_[celli], 
                ftVar_[celli], 
                fu_[celli], 
                ox_[celli], 
                pr_[celli]
            );
        }

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                b_.boundaryField()[patchi][facei],
                ftVar_.boundaryField()[patchi][facei],
                fu_.boundaryField()[patchi][facei],
                ox_.boundaryField()[patchi][facei],
                pr_.boundaryField()[patchi][facei]
            );
        }

/*
        const ThermoType& cellReactants(const label celli) const
        {
            return mixture(ft_[celli], 1, ftVar_[celli]);
        }

        const ThermoType& patchFaceReactants
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                1,
                ftVar_.boundaryField()[patchi][facei]
            );
        }

        const ThermoType& cellProducts(const label celli) const
        {
            return mixture(ft_[celli], 0, ftVar_[celli]);
        }

        const ThermoType& patchFaceProducts
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                0,
                ftVar_.boundaryField()[patchi][facei]
            );
        }
*/
        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "fireMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
