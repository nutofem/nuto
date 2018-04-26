add_sources(cell/SimpleAssembler.cpp
    constraints/ConstraintCompanion.cpp
    constraints/Constraints.cpp
    dofs/DofNumbering.cpp
    elements/ElementShapeFunctions.cpp
    elements/SpectralShapeFunctions.cpp

    mesh/MeshFem.cpp
    mesh/MeshFemDofConvert.cpp
    mesh/MeshGmsh.cpp
    mesh/UnitMeshFem.cpp

    solver/Solve.cpp

    tools/AdaptiveSolve.cpp
    tools/CellStorage.cpp
    tools/GlobalFractureEnergyIntegrator.cpp
    tools/NodalValueMerger.cpp
    tools/QuasistaticSolver.cpp
    tools/TimeDependentProblem.cpp

    interpolation/InterpolationTriangleLinear.cpp

    integrationtypes/IntegrationCompanion.cpp    
    integrationtypes/IntegrationType0DBoundary.cpp
    integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.cpp
    integrationtypes/IntegrationType3D6NGauss1Ip.cpp
    integrationtypes/IntegrationType3D6NGauss2x3Ip.cpp
    integrationtypes/IntegrationTypeEnum.cpp
    integrationtypes/IntegrationTypeTensorProduct.cpp
    integrationtypes/IntegrationTypeTetrahedron.cpp
    integrationtypes/IntegrationTypeTriangle.cpp
    interpolation/InterpolationCompanion.cpp
    )
