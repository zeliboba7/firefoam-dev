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
  
 Application 
     fireFoam 
  
 Description 
     Transient solver for buoyant, turbulent diffusion flame 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "LESModel.H"

#include "radiationModel.H"
#include "MULES.H"

#include "OFstream.H"

#include "specieThermo.H"
#include "janafThermo.H"
#include "perfectGas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"

    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"

    //#include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readPISOControls.H"
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        #include "readLowMachControls.H"
        #include "readMixingControls.H"
        #include "readMultivarMULEControls.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "computeB.H"

        // Pressure-velocity PIMPLE corrector loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            #include "rhoEqn.H"
            #include "UEqn.H"

            #include "ftEqn.H"
            #include "hEqn.H"

            // --- PISO loop
            for (int corr=0; corr<nCorr; corr++)
            {
                #include "pEqn.H"
            }
        }

        turbulence->correct();

        rho = thermo.rho();

        runTime.write();

        //#include "infoOutput.H"     

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
