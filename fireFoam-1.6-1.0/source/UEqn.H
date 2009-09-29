/*---------------------------------------------------------------------------*\
    This file is part of a preliminary version of the FireFOAM code which is 
    currently under development at FM Global. FM Global designates this file
    WITHOUT ANY WARRANTY to be used for development purpose only. 
\*---------------------------------------------------------------------------*/


    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn.relax();

    if (oCorr == nOuterCorr - 1)
    {
        solve(UEqn == -fvc::grad(pd) - fvc::grad(rho)*gh, mesh.solver("UFinal"));
    }
    else
    {
        solve(UEqn == -fvc::grad(pd) - fvc::grad(rho)*gh);
    }


