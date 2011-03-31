tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_hs)")
    )
);
{
    combustion->correct();
    dQ = combustion->dQ();
    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
    	    fvScalarMatrix R = combustion->R(Yi);

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->alphaEff(), Yi)
            );

            if (oCorr == nOuterCorr-1)
            {
                YiEqn.relax(1);   
            }
            else
            {
                YiEqn.relax(); 
            }

            solve
            (
                YiEqn
             == 
                parcels.Srho(i)
              + R,
                mesh.solver("Yi")
            );

            //Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
        Info << "max of " << Y[i].name() << " = " << max(Y[i]).value() << endl;
        Info << "min of " << Y[i].name() << " = " << min(Y[i]).value() << endl;
    }

    Y[inertIndex] = scalar(1) - Yt;
    //Y[inertIndex].max(0.0);

    radiation->correct();

    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        DpDt
      + parcels.Sh()         
      + dQ                    
      + radiation->Shs(thermo)
    );

    if (oCorr == nOuterCorr-1)
    {
        hsEqn.relax(1);
    }
    else
    {
        hsEqn.relax();
    }

    hsEqn.solve();

    thermo.correct();
}
