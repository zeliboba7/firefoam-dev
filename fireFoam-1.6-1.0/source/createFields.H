
    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<hCombustionThermo> pThermo
    (
        hCombustionThermo::New(mesh)
    );
    hCombustionThermo& thermo = pThermo();

    basicMultiComponentMixture& composition = thermo.composition();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    dimensionedScalar stoicRatio
    (
        thermo.lookup("stoichiometricAirFuelMassRatio")
    );

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();

    const volScalarField& psi = thermo.psi();

    volScalarField& b = composition.Y("b");
    volScalarField& ft = composition.Y("ft");

    volScalarField& fu = composition.Y("fu");
    volScalarField& ox = composition.Y("ox");
    volScalarField& pr = composition.Y("pr");

    IOdictionary fireFoamProperties
    (
        IOobject
        (
            "fireFoamProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar pRef(fireFoamProperties.lookup("pRef"));

    // read ignition flag
    Switch ignited(fireFoamProperties.lookup("ignited"));

    const scalar ignRampTime(readScalar(fireFoamProperties.lookup("ignRampTime")));

    const volScalarField& T = thermo.T();

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Store the chemical enthalpy for use in evaluating the heat-release rate
    volScalarField hc = thermo.hc();
    hc.oldTime();

//    volScalarField hs = thermo.hs();
//    hs.oldTime();

#   include "compressibleCreatePhi.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::LESModel> turbulence
    (
        compressible::LESModel::New(rho, U, phi, thermo)
    );

#   include "infoFieldsOutput.H"

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt = fvc::ddt(p);


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());


    Info<< "Creating field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//    if (lowMach)
//    {
//        p = rho*gh + pRef;
//    }
//    else
//    {
        p = pd + rho*gh + pRef;
//    }


    Info<< "Creating radiation model\n" << endl;
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );


    dimensionedScalar initialMass = fvc::domainIntegrate(rho);

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    if (composition.contains("ft"))
    {
        fields.add(composition.Y("ft"));
    }

    fields.add(h);

/*
//tmp
//    volScalarField alphaOuptput = turbulence->alpha();
    volScalarField alphaOutput
    (
        IOobject
        (
            "alphaOutput",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->alpha()
    );


    volScalarField cpOutput
    (
        IOobject
        (
            "cpOutput",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.Cp()
    );

*/
