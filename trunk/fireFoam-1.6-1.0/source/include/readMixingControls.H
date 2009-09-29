    dictionary dictMixing = mesh.solutionDict().subDict("mixing");

    scalar LESsimilarityC
    (
        readScalar(dictMixing.lookup("LESsimilarityC"))
    );

    scalar lengthScaleC
    (
        readScalar(dictMixing.lookup("lengthScaleC"))
    );

    word ftVarModel
    (
        dictMixing.lookup("ftVarModel")
    );
