    singleStepReactingMixture<gasThermoPhysics>& singleMixture
    (
        dynamic_cast<singleStepReactingMixture<gasThermoPhysics>&>
        (thermo)
    );

    // stoichiometric O2 to fuel ratio
    scalar s(singleMixture.s().value());

    // stoichiometric air to fuel ratio
    scalar stoicRatio(singleMixture.stoicRatio().value());  

    // heat of combustion [J/kg]
    scalar qFuel(singleMixture.qFuel().value());  

    scalar fuelIndex(singleMixture.fuelIndex());

    const volScalarField& O2 = thermo.composition().Y("O2");
    const volScalarField& fu = Y[fuelIndex]; 

    scalar YO2Inf = 0.23301; //hardcode for now

    // compute stoichiometric mixture fraction
    scalar ftSt = 1.0 / ( 1.0 + stoicRatio );
    Info << "stoichiometric mixture fraction is = " << ftSt << endl;

    // create fileds for mixture fraction
    volScalarField ft
    (
        IOobject
        (
            "ft",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fu*s-O2+YO2Inf)/(s+YO2Inf)
    );

    // create fileds for patch integration of HRR
    surfaceScalarField HRR_fu
    ( 
        IOobject
        (
            "HRR_fu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(fu)*qFuel
    );

    // for outputing flame height
    OFstream outFlameHeight("outFlameHeight_"+runTime.timeName());
    outFlameHeight << "# Time    " << "FlameHeight  "  << endl;

    volScalarField flameHeight
    (
        IOobject
        (
            "flameHeight",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensionedScalar("zero",dimLength,0.0)
    );

    const pointField& cellCentres = mesh.cellCentres();


    // for check ft conservation
    surfaceScalarField phiFt
    ( 
        IOobject
        (
            "phiFt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    // create fileds for surface integration of sensible enthalpy flux
    surfaceScalarField phiHs
    (
        IOobject
        (
            "phiHs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(hs) 
    );

    // create fileds for surface integration of chemical enthalpy flux
    const volScalarField hc = thermo.hc();  
    surfaceScalarField phiHc
    (
        IOobject
        (
           "phiHc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );

    // create fileds for surface integration of total enthalpy flux
    volScalarField h = hc + hs;
    surfaceScalarField phiH
    (
        IOobject
        (
            "phiH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );

    // create UT field for average purpose
    volVectorField UT
    (
        IOobject
        (
            "UT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U*T
    );

    // create B field for average purpose
    volSymmTensorField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(U)
    );


