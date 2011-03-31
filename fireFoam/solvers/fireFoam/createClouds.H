Info<< "\nConstructing reacting cloud" << endl;
//basicReactingCloud parcels
thermoReactingCloud parcels
(
    "reactingCloud1",
    rho,
    U,
    g,
    //slgThermo
    thermo
);
