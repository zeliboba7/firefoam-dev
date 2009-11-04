
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,ft_b_h)")
    )
);

if (composition.contains("ft"))
{
    if (isMULES)
    {
        MULES::explicitSolve
        (
            rho,
            ft,
            phi,
            fvc::flux
            (
                phi,
                ft
            )(),
            zeroField(),
            zeroField(),
            1, 0
        );

        solve
        (
            fvm::ddt(rho, ft) - fvc::ddt(rho, ft)
          - fvm::laplacian(turbulence->alphaEff(), ft)
        );
    }
    else
    {
        if (multiVariate)
        {
            solve
            (
                fvm::ddt(rho, ft)
              + mvConvection->fvmDiv(phi, ft)
              - fvm::laplacian(turbulence->alphaEff(), ft)                
            );
        }
        else
        {
            solve
            (
                fvm::ddt(rho, ft)
              + fvm::div(phi, ft)
              - fvm::laplacian(turbulence->alphaEff(), ft)               
            );
        }
    }

    Info<< "max(ft) = " << max(ft).value() << endl;
    Info<< "min(ft) = " << min(ft).value() << endl;
}

#       include "ftVarEqn.H"


