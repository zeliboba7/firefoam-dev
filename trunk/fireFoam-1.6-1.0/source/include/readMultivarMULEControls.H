    dictionary dictMultiVariate =
    mesh.solutionDict().subDict("multiVariate");

    dictionary dictMULES =
    mesh.solutionDict().subDict("MULES");

    Switch multiVariate(dictMultiVariate.lookup("bMultiVariate"));

    Switch isMULES(dictMULES.lookup("bMULES"));

