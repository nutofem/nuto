add_sources(cell/SimpleAssembler.cpp
    constraints/ConstraintCompanion.cpp
    constraints/Constraints.cpp
    dofs/DofNumbering.cpp
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

    interpolation/InterpolationTrussLinear.cpp
    interpolation/InterpolationTrussQuadratic.cpp
    interpolation/InterpolationTruss3rdOrder.cpp
    interpolation/InterpolationTruss4thOrder.cpp
    interpolation/InterpolationTriangleLinear.cpp
    interpolation/InterpolationTriangleQuadratic.cpp
    interpolation/InterpolationTriangle3rdOrder.cpp
    interpolation/InterpolationTriangle4thOrder.cpp
    interpolation/InterpolationTetrahedronLinear.cpp
    interpolation/InterpolationTetrahedronQuadratic.cpp
    interpolation/InterpolationQuadLinear.cpp
    interpolation/InterpolationQuadSerendipity.cpp
    interpolation/InterpolationBrickLinear.cpp
    interpolation/InterpolationBrickQuadratic.cpp
    interpolation/InterpolationPrismLinear.cpp
    interpolation/InterpolationPrismQuadratic.cpp
    interpolation/InterpolationPyramidLinear.cpp

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
