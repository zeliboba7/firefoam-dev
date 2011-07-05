/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "basicReactingParcel.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeParcelDispersionModels.H"
#include "makeParcelDragModels.H"
#include "makeReactingParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(basicReactingParcel);

    // Kinematic sub-models
    makeParcelDispersionModels(basicReactingParcel);
    makeParcelDragModels(basicReactingParcel);
    makeReactingParcelInjectionModels(basicReactingParcel);
    makeParcelPatchInteractionModels(basicReactingParcel);

    // Thermo sub-models
    makeParcelHeatTransferModels(basicReactingParcel);

    // Reacting sub-models
    makeReactingParcelCompositionModels(basicReactingParcel);
    makeReactingParcelPhaseChangeModels(basicReactingParcel);
    makeReactingParcelSurfaceFilmModels(basicReactingParcel);
};


// ************************************************************************* //
