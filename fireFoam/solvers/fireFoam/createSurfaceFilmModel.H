Info<< "\nConstructing surface film model" << endl;

autoPtr<regionModels::surfaceFilmModels::surfaceFilmModel>
    tsurfaceFilm
    (
        regionModels::surfaceFilmModels::surfaceFilmModel::New(mesh, g)
    );

regionModels::surfaceFilmModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();

