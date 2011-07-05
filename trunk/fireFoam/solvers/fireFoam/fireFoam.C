/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Application
    fireFoam

Description
    Transient Solver for Fires and turbulent diffusion flames

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hsCombustionThermo.H"
#include "turbulenceModel.H"
#include "combustionModel.H"
#include "basicReactingCloud.H"
#include "surfaceFilmModel.H"
#include "radiationModel.H"
#include "solidChemistryModel.H"
#include "pyrolysisModel.H"
#include "MULES.H"

#include "singleStepReactingMixture.H"
#include "thermoPhysicsTypes.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    //#include "printVersion.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createRadiationModel.H"
    #include "createPyrolysisModel.H"
    #include "createClouds.H"
    #include "createSurfaceFilmModel.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    #include "readPyrolysisTimeControls.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readPISOControls.H"
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        //#include "setMultiRegionDeltaT.H"
        #include "setDeltaT.H"
        #include "readMultivarMULEControls.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        surfaceFilm.evolve();

        if (solvePrimaryRegion)
        {
             #include "rhoEqn.H"

             // --- Pressure-velocity PIMPLE corrector loop
             for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
             {
                 #include "UEqn.H"
                 #include "YhsEqn.H"

                 // --- PISO loop
                 for (int corr=0; corr<nCorr; corr++)
                 {
                     //#include "pEqn.H"
                     #include "p_rghEqn.H"
                 }

                 turbulence->correct();

                 if (oCorr == nOuterCorr-1)
                 {
                     #include "infoOutput.H"
                 }
             }

             //turbulence->correct();

             rho = thermo.rho();

//             #include "infoOutput.H" 

             pyrolysis.evolve();

             runTime.write();
        }
        else
        {
            pyrolysis.evolve();
            runTime.write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s";
        if(runTime.outputTime()){
            Info<< " +";
        }
        Info<< nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
