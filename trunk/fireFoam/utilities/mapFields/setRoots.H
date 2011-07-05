    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("source dir");

    argList::validOptions.insert("sourceTime", "scalar");
    argList::validOptions.insert("sourceRegion", "word");
    argList::validOptions.insert("targetRegion", "word");
    argList::validOptions.insert("parallelSource", "");
    argList::validOptions.insert("parallelTarget", "");
    argList::validOptions.insert("consistent", "");
    argList::validOptions.insert("mapMethod", "word");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath(args.additionalArgs()[0]);
    fileName rootDirSource = casePath.path();
    fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;

    word sourceRegion = fvMesh::defaultRegion;
    if (args.optionFound("sourceRegion"))
    {
        sourceRegion = args.option("sourceRegion");
        Info<< "Source region: " << sourceRegion << endl;
    }

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
    word targetRegion = fvMesh::defaultRegion;
    if (args.optionFound("targetRegion"))
    {
        targetRegion = args.option("targetRegion");
        Info<< "Target region: " << targetRegion << endl;
    }

    bool parallelSource = args.optionFound("parallelSource");
    bool parallelTarget = args.optionFound("parallelTarget");
    bool consistent = args.optionFound("consistent");


    meshToMesh::order mapOrder = meshToMesh::INTERPOLATE;
    if (args.optionFound("mapMethod"))
    {
        const word mapMethod(args.option("mapMethod"));
        if (mapMethod == "mapNearest")
        {
            mapOrder = meshToMesh::MAP;
        }
        else if (mapMethod == "interpolate")
        {
            mapOrder = meshToMesh::INTERPOLATE;
        }
        else if (mapMethod == "cellPointInterpolate")
        {
            mapOrder = meshToMesh::CELL_POINT_INTERPOLATE;
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Unknown mapMethod " << mapMethod << ". Valid options are: "
                << "mapNearest, interpolate and cellPointInterpolate"
                << exit(FatalError);
        }

        Info<< "Mapping method: " << mapMethod << endl;
    }
