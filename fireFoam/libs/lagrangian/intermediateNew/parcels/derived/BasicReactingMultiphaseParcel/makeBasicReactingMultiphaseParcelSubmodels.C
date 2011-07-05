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

#include "basicReactingMultiphaseParcel.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeParcelDispersionModels.H"
#include "makeParcelDragModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP variant
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(basicReactingMultiphaseParcel);

    // Kinematic sub-models
    makeParcelDispersionModels(basicReactingMultiphaseParcel);
    makeParcelDragModels(basicReactingMultiphaseParcel);
    makeReactingMultiphaseParcelInjectionModels(basicReactingMultiphaseParcel);
    makeParcelPatchInteractionModels(basicReactingMultiphaseParcel);

    // Thermo sub-models
    makeParcelHeatTransferModels(basicReactingMultiphaseParcel);

    // Reacting sub-models
    makeReactingMultiphaseParcelCompositionModels(basicReactingMultiphaseParcel);
    makeReactingParcelPhaseChangeModels(basicReactingMultiphaseParcel);

    // Reacting multiphase sub-models
    makeReactingMultiphaseParcelDevolatilisationModels
    (
        basicReactingMultiphaseParcel
    );
    makeReactingParcelSurfaceFilmModels
    (
        basicReactingMultiphaseParcel
    );
    makeReactingMultiphaseParcelSurfaceReactionModels
    (
        basicReactingMultiphaseParcel
    );
};


// ************************************************************************* //
