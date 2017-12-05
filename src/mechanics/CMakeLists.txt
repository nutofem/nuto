set(MechanicsConstraintSources
    constraintsPde/Constraints.cpp
    constraintsPde/ConstraintCompanion.cpp
    )

set(MechanicsDofsSources
    dofs/DofNumbering.cpp
    )

list(APPEND MechanicsElementsSources
    elements/ElementShapeFunctions.cpp
    )

set(MechanicsIntegrationtypesSources
    integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.cpp
    integrationtypes/IntegrationType0DBoundary.cpp
    integrationtypes/IntegrationType2D3NGauss13Ip.cpp
    integrationtypes/IntegrationType2D3NGauss16Ip.cpp
    integrationtypes/IntegrationType2D3NGauss1Ip.cpp
    integrationtypes/IntegrationType2D3NGauss3Ip.cpp
    integrationtypes/IntegrationType2D3NGauss4Ip.cpp
    integrationtypes/IntegrationType2D3NGauss6Ip.cpp
    integrationtypes/IntegrationType2D3NGauss12Ip.cpp
    integrationtypes/IntegrationType2D3NGauss12IpDetail.cpp
    integrationtypes/IntegrationType3D4NGauss1Ip.cpp
    integrationtypes/IntegrationType3D4NGauss4Ip.cpp
    integrationtypes/IntegrationType3D6NGauss1Ip.cpp
    integrationtypes/IntegrationType3D6NGauss2x3Ip.cpp
    integrationtypes/IntegrationTypeTensorProduct.cpp
    integrationtypes/IntegrationTypeBase.cpp
    integrationtypes/IntegrationTypeEnum.cpp
    )

set(MechanicsMeshSources
    mesh/MeshFem.cpp
    mesh/MeshFemDofConvert.cpp
    mesh/UnitMeshFem.cpp
    )

set(MechanicsStructuresSources
    structures/StructureBase.cpp
    structures/StructureBaseConstitutiveLaw.cpp
    structures/StructureBaseElement.cpp
    structures/StructureBaseEnum.cpp
    structures/StructureBaseGroup.cpp
    structures/StructureBaseIntegrationTypes.cpp
    structures/StructureBaseLoad.cpp
    structures/StructureBaseLogger.cpp
    structures/StructureBaseNode.cpp
    structures/StructureOutputBase.cpp
    structures/StructureOutputBlockMatrix.cpp
    structures/StructureOutputBlockVector.cpp
    structures/unstructured/Structure.cpp
    structures/unstructured/StructureElement.cpp
    structures/unstructured/StructureGroup.cpp
    structures/unstructured/StructureInterpolationType.cpp
    structures/unstructured/StructureNode.cpp
    structures/Assembler.cpp
    )



set(MechanicsTimeIntegrationSources
    timeIntegration/TimeControl.cpp
    )

set(MechanicsToolsSources
    tools/GlobalFractureEnergyIntegrator.cpp
    )

set(MechanicsSources
    ${MechanicsConstitutiveSources}
    ${MechanicsConstraintSources}
    ${MechanicsDofsSources}
    ${MechanicsElementsSources}
    ${MechanicsIntegrationtypesSources}
    ${MechanicsMeshSources}
    ${MechanicsStructuresSources}
    ${MechanicsTimeIntegrationSources}
    ${MechanicsToolsSources}
    )

create_nuto_module(Mechanics "${MechanicsSources}")
target_link_libraries(Mechanics Base Math Boost::filesystem Ann::Ann)

if(ENABLE_MPI)
    message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
    target_link_libraries(Mechanics ${MPI_LIBRARIES})
    message(STATUS "Creating feti library...")
    add_subdirectory(feti)
endif()
