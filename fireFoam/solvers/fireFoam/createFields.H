Info<< "Reading thermophysical properties\n" << endl;

autoPtr<hsCombustionThermo> pThermo
(
    hsCombustionThermo::New(mesh)
);

hsCombustionThermo& thermo = pThermo();

//SLGThermo slgThermo(mesh, thermo);  

basicMultiComponentMixtureNew& composition = thermo.composition();

PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

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


volScalarField& p = thermo.p();
volScalarField& hs = thermo.hs();
const volScalarField& T = thermo.T();

const volScalarField& psi = thermo.psi();

Info<< "Reading field p_rgh\n" << endl;

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New(rho, U, phi, thermo)
);

IOdictionary combustionProperties
(
    IOobject
    (
        "combustionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Creating combustion model\n" << endl;
//autoPtr<combustionModel::combustionModel> combustion
autoPtr<combustionModel> combustion
(
    combustionModel::combustionModel::New
    (
        combustionProperties,
        thermo,
        turbulence(),
        phi,
        rho
    )
);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimMass/pow3(dimTime)/dimLength, 0.0)
);


#include "infoFieldsOutput.H"


Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt =
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);


Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());

surfaceScalarField ghf("gh", g & mesh.Cf());

//p += rho*gh;
//thermo.correct();

p = p_rgh + rho*gh;
thermo.correct();
rho = thermo.rho();
p_rgh = p - rho*gh;

dimensionedScalar initialMass = fvc::domainIntegrate(rho);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y,i)
{
    fields.add(Y[i]);
}

fields.add(hs);

IOdictionary additionalControlsDict
(
    IOobject
    (
        "additionalControls",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Switch solvePrimaryRegion
(
    additionalControlsDict.lookup("solvePrimaryRegion")
);


