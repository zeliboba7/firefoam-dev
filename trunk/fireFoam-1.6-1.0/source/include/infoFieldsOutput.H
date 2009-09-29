    //OFstream outHRR(fileName);
    OFstream outHRR("outHRR_"+runTime.timeName());
    outHRR << "# Time    " << "HRR    " << endl;

    // compute stoichiometric mixture fraction
    scalar ftSt = 1.0 / ( 1.0 + stoicRatio.value() );
    Info << "stoichiometric mixture fraction is = " << ftSt << endl;

    // compute heat of combustion
    specieThermo<janafThermo<perfectGas> >  fuel(thermo.lookup("fuel"));
    specieThermo<janafThermo<perfectGas> >  oxidant(thermo.lookup("oxidant"));
    specieThermo<janafThermo<perfectGas> >  products(thermo.lookup("burntProducts"));

    dimensionedScalar hoC
    (  
        "hoC",
        dimVelocity*dimVelocity, 
        - products.hc()*(1.0+stoicRatio.value())/products.W()
        + fuel.hc()/fuel.W() 
        + oxidant.hc()*stoicRatio.value()/oxidant.W()
    );

    Info << "heat of combustion of fuel is " << hoC << endl;

    // create fields for HRR estimation 
    volScalarField HRR
    (
        IOobject
        (
            "HRR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(rho,h)
    );
    
    // filed used by fvm radiation solver
    volScalarField hrr 
    (
        IOobject
        (
            "hrr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::ddt(rho,h)
    );

