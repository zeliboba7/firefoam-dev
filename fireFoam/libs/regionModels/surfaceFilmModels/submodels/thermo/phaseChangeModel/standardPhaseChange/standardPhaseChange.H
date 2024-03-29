/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::standardPhaseChange

Description
    Standard phase change model with modification for boiling

SourceFiles
    standardPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef standardPhaseChange_H
#define standardPhaseChange_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class standardPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class standardPhaseChange
:
    public phaseChangeModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        standardPhaseChange(const standardPhaseChange&);

        //- Disallow default bitwise assignment
        void operator=(const standardPhaseChange&);


protected:

    // Protected data

        //- Boiling temperature / [K]
        const scalar Tb_;

        //- Minimum film height for model to be active
        const scalar deltaMin_;

        //- Length scale / [m]
        const scalar L_;


    // Protected member functions

        //- Return Sherwood number as a function of Reynolds and Schmidt numbers
        scalar Sh(const scalar Re, const scalar Sc) const;


public:

    //- Runtime type information
    TypeName("standardPhaseChange");


    // Constructors

        //- Construct from surface film model
        standardPhaseChange
        (
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~standardPhaseChange();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
