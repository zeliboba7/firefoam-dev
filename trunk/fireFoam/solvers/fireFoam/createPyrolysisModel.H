Info<< "Creating pyrolysis model\n" << endl;

autoPtr<regionModels::pyrolysisModels::pyrolysisModel> tpyrolysis
(
    regionModels::pyrolysisModels::pyrolysisModel::New(mesh)
);
regionModels::pyrolysisModels::pyrolysisModel& pyrolysis = tpyrolysis();

//scalar tau = runTime.deltaT().value();
