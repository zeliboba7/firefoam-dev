/*---------------------------------------------------------------------------*\
    This file is part of a preliminary version of the FireFOAM code which is 
    currently under development at FM Global. FM Global designates this file
    WITHOUT ANY WARRANTY to be used for development purpose only. 
\*---------------------------------------------------------------------------*/

if (composition.contains("ftVar"))
{
    volScalarField& ftVar = composition.Y("ftVar");

    if (ftVarModel == "LESsimilarity")
    {
        ftVar = LESsimilarityC*turbulence->delta()*turbulence->delta()
        *magSqr(fvc::grad(ft));
    }
    else if (ftVarModel == "lengthScale")
    {
        ftVar = 1.0/lengthScaleC*turbulence->k()/turbulence->epsilon()
        *turbulence->alphaEff()/rho*magSqr(fvc::grad(ft));
    }
    else if (ftVarModel == "transport")
    {
        if (multiVariate&&!isMULES)
        {
            solve
            (
                fvm::ddt(rho, ftVar)
              + mvConvection->fvmDiv(phi, ftVar)
              - fvm::laplacian(turbulence->alphaEff(), ftVar)
              ==
                2.*turbulence->alphaEff()*magSqr(fvc::grad(ft))
              - lengthScaleC*fvm::SuSp(2.0*rho*turbulence->epsilon()/
                turbulence->k(),ftVar)
            );
        }
        else
        {
            solve
            (
                fvm::ddt(rho, ftVar)
              + fvm::div(phi, ftVar)
              - fvm::laplacian(turbulence->alphaEff(), ftVar)
              ==
                2.*turbulence->alphaEff()*magSqr(fvc::grad(ft))
              - lengthScaleC*fvm::SuSp(2.0*rho*turbulence->epsilon()/
                turbulence->k(),ftVar)
            );
        }
    }
    else
    {
        FatalError
            << args.executable() << " : Unknown ftVar model " << ftVarModel
            << abort(FatalError);
    }

    Info<< "max(ftVar) = " << max(ftVar).value() << endl;
    Info<< "min(ftVar) = " << min(ftVar).value() << endl;

}
