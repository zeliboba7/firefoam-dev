
bool closedVolume = false;

rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

phi =
    fvc::interpolate(rho)
   *(
       (fvc::interpolate(U) & mesh.Sf())
     + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

phi -= fvc::interpolate(rUA*rho*gh)*fvc::snGrad(rho)*mesh.magSf();

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pdEqn
    (
        fvc::ddt(psi)*pRef
      + fvc::ddt(psi, rho)*gh
      + fvc::div(phi)
      - fvm::laplacian(rho*rUA, pd)
    );

    if (!lowMach)
    {  
        pdEqn += fvm::ddt(psi,pd);
    }

    closedVolume = pd.needReference();

    if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
    {
        pdEqn.solve(mesh.solver(pd.name() + "Final"));
    }
    else
    {
        pdEqn.solve(mesh.solver(pd.name()));
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi += pdEqn.flux();
    }
}

if (lowMach)
{
  p = rho*gh + pRef;
}
else
{
  p = pd + rho*gh + pRef;
}

dpdt = fvc::ddt(p);

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rUA*(fvc::grad(pd) + fvc::grad(rho)*gh);

U.correctBoundaryConditions();


// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(thermo.psi()*p))
        /fvc::domainIntegrate(thermo.psi());
    rho = thermo.rho();
}
