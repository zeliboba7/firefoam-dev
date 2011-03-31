
        ft = (fu*s-O2+YO2Inf)/(s+YO2Inf);

        HRR_fu = ( - mvConvection->interpolate(phi,fu)*phi
                   + (
                         fvc::interpolate
                         (
                             turbulence->alphaEff()*fvc::grad(fu)
                         )
                         & mesh.Sf()
                     )                         
                 ) * qFuel;


        if (mag(g).value() > 0.0) 
        {    
            flameHeight = dimensionedScalar("zero",dimLength,0.0);

            forAll (ft.internalField(),cellI)
            {
                if (ft[cellI] >= ftSt)
                { 
                    flameHeight[cellI] = mag( cellCentres[cellI] & (g/mag(g)).value() );
                }
            }

            outFlameHeight <<  runTime.value() << "  "     
                    << max(flameHeight).value() << "  "
                    <<  endl;
        }


        volScalarField hc = thermo.hc();  
        h = hs + hc;

        phiFt = mvConvection->interpolate(phi,ft)*phi
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(ft)*mesh.magSf();

        phiHc = phi*mvConvection->interpolate(phi,hc)
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(hc)*mesh.magSf();

        phiHs = phi*mvConvection->interpolate(phi,hs)
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(hs)*mesh.magSf();

        phiH  = phi*mvConvection->interpolate(phi,h)
                 - fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(h)*mesh.magSf();


        UT = U * T;

        B = turbulence->R();


    




