    dictionary dictcompressibleVersion =
    mesh.solutionDict().subDict("compressibleVersion");

//    word lowMach
    Switch lowMach
    (
        dictcompressibleVersion.lookup("bLowMach")
    );
