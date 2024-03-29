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
    Foam::odeNew

Description
    An ODE solver for chemistry

SourceFiles
    odeNew.C

\*---------------------------------------------------------------------------*/

#ifndef odeNew_H
#define odeNew_H

#include "chemistrySolver.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                            Class odeNew Declaration
\*---------------------------------------------------------------------------*/

template<class ODEChemistryType>
class odeNew
:
    public chemistrySolver<ODEChemistryType>
{
    // Private data

        dictionary coeffsDict_;
        const word solverName_;
        autoPtr<ODESolver> odeNewSolver_;

        // ModeNewl constants

            scalar eps_;


public:

    //- Runtime type information
    TypeName("ODE");


    // Constructors

        //- Construct from components
        odeNew
        (
            const fvMesh& mesh,
            const word& ODEmodeNewlName,
            const word& thermoType
        );


    //- Destructor
    virtual ~odeNew();


    // Member Functions

        virtual scalar solve
        (
            scalarField& c,
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
#   include "odeNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
