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

\*---------------------------------------------------------------------------*/

#include "hCombustionThermo.H"
#include "hPsiMixtureThermo.H"

#include "makeCombustionThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "perfectGas.H"
#include "specieThermo.H"
#include "janafThermo.H"
#include "sutherlandTransport.H"

#include "fireMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeCombustionThermo
(
    hCombustionThermo,
    hPsiMixtureThermo,
    fireMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
