/*---------------------------------------------------------------------------*\
    This file is part of a preliminary version of the FireFOAM code which is 
    currently under development at FM Global. FM Global designates this file
    WITHOUT ANY WARRANTY to be used for development purpose only. 
\*---------------------------------------------------------------------------*/

{
  if (multiVariate&&!isMULES)
  {

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi,h)
      - fvm::laplacian(turbulence->alphaEff(), h)   
     ==
        dpdt
      + fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p))
      - p*fvc::div(phi/fvc::interpolate(rho))
      + radiation->Sh(thermo)
    );

    hEqn.relax();
    hEqn.solve();
  }
  else
  {
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)   
     ==
        dpdt
      + fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p))
      - p*fvc::div(phi/fvc::interpolate(rho))
      + radiation->Sh(thermo)
    );

    hEqn.relax();
    hEqn.solve();
  }

  # include "computeHRR.H"

  thermo.correct();

  radiation->correct();
}
