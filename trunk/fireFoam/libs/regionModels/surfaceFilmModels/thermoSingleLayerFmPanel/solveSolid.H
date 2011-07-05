{
    static volScalarField rhos
        (
         IOobject
         (
          "rhos",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.density,
        zeroGradientFvPatchScalarField::typeName
        );
    static volScalarField cps
        (
         IOobject
         (
          "cps",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.cp,
        zeroGradientFvPatchScalarField::typeName
        );
    static volScalarField ks
        (
         IOobject
         (
          "ks",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.k,
        zeroGradientFvPatchScalarField::typeName
        );

    tmp<fvScalarMatrix> TEqn
        (
         fvm::ddt(rhos*cps, Tw_)
         - fvm::laplacian(ks, Tw_)
         == qDotFilm/lc.thickness
        );

    TEqn().relax();
    TEqn().solve();

}
