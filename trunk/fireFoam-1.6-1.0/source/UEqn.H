
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


