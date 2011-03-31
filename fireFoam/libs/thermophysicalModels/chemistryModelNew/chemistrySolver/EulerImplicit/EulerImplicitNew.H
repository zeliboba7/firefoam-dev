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
    Foam::EulerImplicitNew

Description
    An Euler implicit solver for chemistry

SourceFiles
    EulerImplicitNew.C

\*---------------------------------------------------------------------------*/

#ifndef EulerImplicitNew_H
#define EulerImplicitNew_H

#include "chemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                       Class EulerImplicitNew Declaration
\*---------------------------------------------------------------------------*/

template<class ODEChemistryType>
class EulerImplicitNew
:
    public chemistrySolver<ODEChemistryType>
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;


        // Model constants

            //- Chemistry timescale
            scalar cTauChem_;

            //- Equilibrium rate limiter flag (on/off)
            Switch eqRateLimiter_;


public:

    //- Runtime type information
    TypeName("EulerImplicitNew");


    // Constructors

        //- Construct from components
        EulerImplicitNew
        (
            const fvMesh& mesh,
            const word& ODEmodelName,
            const word& thermoType
        );


    //- Destructor
    virtual ~EulerImplicitNew();


    // Member Functions

        //- Update the concentrations and return the chemical time
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

#ifdef NoRepository
#   include "EulerImplicitNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
