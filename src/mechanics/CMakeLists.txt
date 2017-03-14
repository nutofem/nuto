set(NuToMechanicsConstitutiveSources
    constitutive/ConstitutiveBase.cpp
    constitutive/ConstitutiveEnum.cpp
    constitutive/inputoutput/ConstitutiveIOBase.cpp
    constitutive/inputoutput/ConstitutiveIOMap.cpp
    constitutive/inputoutput/EngineeringStrain.cpp
    constitutive/inputoutput/EngineeringStress.cpp
    constitutive/inputoutput/EquivalentStrain.cpp
    constitutive/laws/AdditiveBase.cpp
    constitutive/laws/AdditiveInputExplicit.cpp
    constitutive/laws/AdditiveInputImplicit.cpp
    constitutive/laws/AdditiveOutput.cpp
    constitutive/laws/EngineeringStressHelper.cpp
    constitutive/laws/FibreMatrixBondStressSlip.cpp
    constitutive/laws/GradientDamageEngineeringStress.cpp
    constitutive/laws/GradientDamageFatigueEngineeringStress.cpp
    constitutive/laws/PhaseField.cpp
    constitutive/laws/HeatConduction.cpp
    constitutive/laws/LinearDampingEngineeringStress.cpp
    constitutive/laws/LinearElasticEngineeringStress.cpp
    constitutive/laws/LocalDamageModel.cpp
    constitutive/laws/MisesPlasticityEngineeringStress.cpp
    constitutive/laws/MoistureTransport.cpp
    constitutive/laws/ShrinkageCapillaryStrainBased.cpp
    constitutive/laws/ShrinkageCapillaryStressBased.cpp
    constitutive/laws/ThermalStrains.cpp
    constitutive/staticData/DataAdditiveInputImplicit.cpp
    constitutive/staticData/DataMisesPlasticity.cpp
    constitutive/staticData/DataMoistureTransport.cpp
    constitutive/staticData/IPAdditiveInputExplicit.cpp
    constitutive/staticData/IPAdditiveInputImplicit.cpp
    constitutive/staticData/IPAdditiveOutput.cpp
    )

set(NuToMechanicsConstraintSources
    constraints/ConstraintBase.cpp
    constraints/ConstraintEquationTerm.cpp
    constraints/ConstraintLinear.cpp
    constraints/ConstraintLinearDisplacementsPeriodic2D.cpp
    constraints/ConstraintLinearEquation.cpp
    constraints/ConstraintLinearNodeDisplacements1D.cpp
    constraints/ConstraintLinearNodeDisplacements2D.cpp
    constraints/ConstraintLinearNodeDisplacements3D.cpp
    constraints/ConstraintLinearNodeGroupDisplacements1D.cpp
    constraints/ConstraintLinearNodeGroupDisplacements2D.cpp
    constraints/ConstraintLinearNodeGroupDisplacements3D.cpp
    constraints/ConstraintLinearNodeGroupRotations2D.cpp
    constraints/ConstraintLinearNodeGroupTemperature.cpp
    constraints/ConstraintLinearNodeRelativeHumidity.cpp
    constraints/ConstraintLinearNodeRotations2D.cpp
    constraints/ConstraintLinearNodeTemperature.cpp
    constraints/ConstraintLinearNodeWaterVolumeFraction.cpp
    constraints/ConstraintNode.cpp
    constraints/ConstraintNodeGroup.cpp
    constraints/ConstraintNonlinear.cpp
    )

set(NuToMechanicsDofSubMatrixStorageSources
    dofSubMatrixStorage/DofStatus.cpp
    dofSubMatrixStorage/BlockStorageBase.cpp
    dofSubMatrixStorage/BlockSparseMatrix.cpp
    dofSubMatrixStorage/BlockScalar.cpp
    dofSubMatrixStorage/BlockFullVector.cpp
    dofSubMatrixStorage/BlockFullMatrix.cpp
    )

list(APPEND NuToMechanicsElementsSources
    elements/ContinuumElement.cpp
    elements/ContinuumBoundaryElement.cpp
    elements/ContinuumBoundaryElementConstrainedControlNode.cpp
    elements/ContinuumContactElement.cpp
    elements/ElementBase.cpp
    elements/Element1DInXD.cpp
    elements/Element2DInterface.cpp
    elements/ContinuumElementIGA.cpp
    elements/ElementOutputBase.cpp
    elements/ElementShapeFunctions.cpp
    elements/IpDataEnum.cpp
    elements/IPData.cpp
    )

set(NuToMechanicsGroupsSources
    groups/Group.cpp
    groups/GroupBase.cpp
    )

set(NuToMechanicsIntegrationtypesSources
    integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.cpp
    integrationtypes/IntegrationType0DBoundary.cpp
    integrationtypes/IntegrationType1D2NGauss1Ip.cpp
    integrationtypes/IntegrationType1D2NGauss2Ip.cpp
    integrationtypes/IntegrationType1D2NGauss3Ip.cpp
    integrationtypes/IntegrationType1D2NGauss4Ip.cpp
    integrationtypes/IntegrationType1D2NGauss5Ip.cpp
    integrationtypes/IntegrationType1D2NLobatto3Ip.cpp
    integrationtypes/IntegrationType1D2NLobatto4Ip.cpp
    integrationtypes/IntegrationType1D2NLobatto5Ip.cpp
    integrationtypes/IntegrationType1D2NLobatto6Ip.cpp
    integrationtypes/IntegrationType2D3NGauss13Ip.cpp
    integrationtypes/IntegrationType2D3NGauss16Ip.cpp
    integrationtypes/IntegrationType2D3NGauss1Ip.cpp
    integrationtypes/IntegrationType2D3NGauss3Ip.cpp
    integrationtypes/IntegrationType2D3NGauss4Ip.cpp
    integrationtypes/IntegrationType2D3NGauss6Ip.cpp
    integrationtypes/IntegrationType2D3NGauss12Ip.cpp
    integrationtypes/IntegrationType2D3NGauss12IpDetail.cpp
    integrationtypes/IntegrationType2D4NGauss1Ip.cpp
    integrationtypes/IntegrationType2D4NGauss4Ip.cpp
    integrationtypes/IntegrationType2D4NGauss9Ip.cpp
    integrationtypes/IntegrationType2D4NLobatto9Ip.cpp
    integrationtypes/IntegrationType2D4NLobatto16Ip.cpp
    integrationtypes/IntegrationType2D4NLobatto25Ip.cpp
    integrationtypes/IntegrationType3D4NGauss1Ip.cpp
    integrationtypes/IntegrationType3D4NGauss4Ip.cpp
    integrationtypes/IntegrationType3D6NGauss2x3Ip
    integrationtypes/IntegrationType3D8NGauss1Ip.cpp
    integrationtypes/IntegrationType3D8NGauss2x2x2Ip.cpp
    integrationtypes/IntegrationTypeBase.cpp
    integrationtypes/IntegrationTypeEnum.cpp
    )

set(NuToMechanicsInterpolationtypesSources
    interpolationtypes/InterpolationBase.cpp
    interpolationtypes/InterpolationBaseFEM.cpp
    interpolationtypes/InterpolationBaseIGA.cpp
    interpolationtypes/InterpolationType.cpp
    interpolationtypes/InterpolationTypeEnum.cpp
    interpolationtypes/Interpolation1D.cpp
    interpolationtypes/Interpolation1DTruss.cpp
    interpolationtypes/Interpolation1DInterface.cpp
    interpolationtypes/Interpolation2D.cpp
    interpolationtypes/Interpolation2DTriangle.cpp
    interpolationtypes/Interpolation2DQuad.cpp
    interpolationtypes/Interpolation3D.cpp
    interpolationtypes/Interpolation3DTetrahedron.cpp
    interpolationtypes/Interpolation3DBrick.cpp
    interpolationtypes/Interpolation3DPrism.cpp
    interpolationtypes/Interpolation1DIGA.cpp
    interpolationtypes/Interpolation2DIGA.cpp
    )

set(NuToMechanicsIGASources
    IGA/BSplineCurve.cpp
    IGA/BSplineSurface.cpp
    )

set(NuToMechanicsLoadsSources
    loads/LoadBase.cpp
    loads/LoadNode.cpp
    loads/LoadNodeGroup.cpp
    loads/LoadNodeForces1D.cpp
    loads/LoadNodeForces2D.cpp
    loads/LoadNodeForces3D.cpp
    loads/LoadNodeGroupForces1D.cpp
    loads/LoadNodeGroupForces2D.cpp
    loads/LoadNodeGroupForces3D.cpp
    loads/LoadNodeHeatFlux1D.cpp
    loads/LoadSurfaceBase2D.cpp
    loads/LoadSurfaceBase3D.cpp
    loads/LoadSurfaceConstDirection2D.cpp
    loads/LoadSurfaceConstDirection3D.cpp
    loads/LoadSurfacePressure2D.cpp
    loads/LoadSurfacePressure3D.cpp
    loads/LoadSurfacePressureFunction2D.cpp
    )

set(NuToMechanicsMeshSources
    mesh/ConvertToInterpolationType.cpp
    mesh/ElementPrismCreate.cpp
    mesh/ImportFromGmsh.cpp
    mesh/MeshCompanion.cpp
    mesh/MeshGenerator.cpp
    )

set(NuToMechanicsNodesSources
    nodes/NodeBase.cpp
    nodes/NodeDof.cpp
    nodes/NodeEnum.cpp
    )

set(NuToMechanicsSectionSources
    sections/SectionBase.cpp
    sections/SectionTruss.cpp
    sections/SectionPlane.cpp
    sections/SectionVolume.cpp
    sections/SectionSpring.cpp
    sections/SectionFibreMatrixBond.cpp
    )

set(NuToMechanicsStructuresSources
    structures/StructureBase.cpp
    structures/StructureBaseConstitutiveLaw.cpp
    structures/StructureBaseConstraints.cpp
    structures/StructureBaseElement.cpp
    structures/StructureBaseEnum.cpp
    structures/StructureBaseGroup.cpp
    structures/StructureBaseIntegrationTypes.cpp
    structures/StructureBaseLoad.cpp
    structures/StructureBaseLogger.cpp
    structures/StructureBaseNode.cpp
    structures/StructureBaseSection.cpp
    structures/StructureOutputBase.cpp
    structures/StructureOutputBlockMatrix.cpp
    structures/StructureOutputBlockVector.cpp
    structures/unstructured/Structure.cpp    
    structures/unstructured/StructureElement.cpp
    structures/unstructured/StructureGroup.cpp
    structures/unstructured/StructureInterpolationType.cpp
    structures/unstructured/StructureNode.cpp
    )



set(NuToMechanicsTimeIntegrationSources
    timeIntegration/NewmarkBase.cpp
    timeIntegration/NewmarkDirect.cpp    
    #timeIntegration/JumpDirect.cpp
    timeIntegration/NystroemBase.cpp
    timeIntegration/NystroemQinZhu.cpp
    timeIntegration/ResultBase.cpp
    timeIntegration/ResultElementIpData.cpp
    timeIntegration/ResultGroupNodeForce.cpp
    timeIntegration/ResultGroupNodeDof.cpp
    timeIntegration/ResultNodeAcceleration.cpp
    timeIntegration/ResultNodeDisp.cpp
    timeIntegration/ResultNodeDof.cpp
    timeIntegration/ResultTime.cpp
    timeIntegration/RungeKutta2.cpp
    timeIntegration/RungeKutta3.cpp
    timeIntegration/RungeKutta38.cpp
    timeIntegration/RungeKutta4.cpp
    timeIntegration/RungeKuttaCashKarp.cpp
    timeIntegration/RungeKuttaDormandPrince.cpp
    timeIntegration/RungeKuttaBase.cpp
    timeIntegration/TimeDependencyBase.cpp
    timeIntegration/TimeIntegrationBase.cpp
    timeIntegration/ImplicitExplicitBase.cpp
    timeIntegration/ImplEx.cpp
    timeIntegration/ImplExCallback.cpp
    timeIntegration/VelocityVerlet.cpp
    )


set(NuToMechanicsSources
    ${NuToMechanicsConstitutiveSources}
    ${NuToMechanicsConstraintSources}
    ${NuToMechanicsCrackSources}
    ${NuToMechanicsDofSubMatrixStorageSources}
    ${NuToMechanicsElementsSources}
    ${NuToMechanicsGroupsSources}
    ${NuToMechanicsIntegrationtypesSources}
    ${NuToMechanicsInterpolationtypesSources}
    ${NuToMechanicsIGASources}
    ${NuToMechanicsLoadsSources}
    ${NuToMechanicsMeshSources}
    ${NuToMechanicsNodesSources}
    ${NuToMechanicsSectionSources}
    ${NuToMechanicsStructuresSources}
    ${NuToMechanicsTimeIntegrationSources}
    )


sources_to_objects("${NuToMechanicsSources}" NuToMechanicsObjects)

# build static library
add_library(NuToMechanics ${NuToMechanicsObjects})
target_link_libraries(NuToMechanics ${ANN_LIBRARIES})

target_link_libraries(NuToMechanics jsoncpp)

if(ENABLE_MPI)
    message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
    target_link_libraries(NuToMechanics ${MPI_LIBRARIES})
endif()

# install library
install(TARGETS NuToMechanics
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )

if(ENABLE_MPI)
    message(STATUS "Creating feti library...")
    add_subdirectory(feti)
endif()

