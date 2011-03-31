constant:
    blockMesh
    setSet:
        faceSet f0 new patchToFace movingWall
        faceZoneSet f0 new setToFaceZone f0

0.005:
    extrudeToRegionMesh17x
