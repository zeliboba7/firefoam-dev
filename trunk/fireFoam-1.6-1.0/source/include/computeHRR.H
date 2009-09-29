  hc = thermo.hc(); 
  hs = thermo.hs(); 

  if (multiVariate&&!isMULES)
  {
      HRR =
          -  fvc::ddt(rho, hc)
          - mvConvection->fvcDiv(phi, hc)
          + fvc::laplacian(turbulence->alphaEff(), hc);

      outHRR <<  runTime.value() << "  "     
             <<  fvc::domainIntegrate(HRR).value() << "  "
             <<  endl;
  }
  else  
  {
      HRR =
          -  fvc::ddt(rho, hc)
          - fvc::div(phi, hc)
          + fvc::laplacian(turbulence->alphaEff(), hc);

      outHRR <<  runTime.value() << "  "     
             <<  fvc::domainIntegrate(HRR).value() << "  "
             <<  endl;
  }
  
  // hrr is used in fvm radiation solver
  hrr=HRR;


