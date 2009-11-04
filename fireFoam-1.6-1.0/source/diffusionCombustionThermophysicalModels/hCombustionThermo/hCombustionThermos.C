/*---------------------------------------------------------------------------*\

     This file is part of a preliminary version of the FireFOAM, LES solver for 
     fire applications, which is currently under development at FM Global and 
     OpenCFD. This file is based on OpenFOAM, the Open Source CFD Toolbox.
     
     Copyright (c) 2009 OpenCFD / FM Global        
  
 License 
       
     FireFOAM is free software; you can redistribute it and/or modify it 
     under the terms of the GNU General Public License as published by the 
     Free Software Foundation; either version 2 of the License, or (at your 
     option) any later version. 
  
     FireFOAM is distributed in the hope that it will be useful, but WITHOUT 
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
     for more details. 
  
     You should have received a copy of the GNU General Public License 
     along with FireFOAM; if not, write to the Free Software Foundation, 
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
