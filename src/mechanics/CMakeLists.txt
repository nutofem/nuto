set(MechanicsConstitutiveSources
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
    constitutive/laws/Creep.cpp
    constitutive/laws/EngineeringStressHelper.cpp
    constitutive/laws/FibreMatrixBondStressSlip.cpp
    constitutive/laws/GradientDamageEngineeringStress.cpp
    constitutive/laws/GradientDamageFatigueEngineeringStress.cpp
    constitutive/laws/PhaseField.cpp
    constitutive/laws/HeatConduction.cpp
    constitutive/laws/LinearDampingEngineeringStress.cpp
    constitutive/laws/LinearDielectric.cpp
    constitutive/laws/LinearElasticEngineeringStress.cpp
    constitutive/laws/LinearElasticInhomogeneous.cpp
    constitutive/laws/LinearElasticAnisotropic.cpp
    constitutive/laws/LinearPiezoelectric.cpp
    constitutive/laws/LocalDamageModel.cpp
    constitutive/laws/MisesPlasticityEngineeringStress.cpp
    constitutive/laws/MoistureTransport.cpp
    constitutive/laws/ShrinkageCapillaryStrainBased.cpp
    constitutive/laws/ShrinkageCapillaryStressBased.cpp
    constitutive/laws/ThermalStrains.cpp
    constitutive/staticData/DataAdditiveInputImplicit.cpp
    constitutive/staticData/DataCreep.cpp
    constitutive/staticData/DataMisesPlasticity.cpp
    constitutive/staticData/DataMoistureTransport.cpp
    constitutive/staticData/IPAdditiveInputExplicit.cpp
    constitutive/staticData/IPAdditiveInputImplicit.cpp
    constitutive/staticData/IPAdditiveOutput.cpp
    )

set(MechanicsConstraintSources
    constraints/Constraints.cpp
    constraints/ConstraintCompanion.cpp
    )

set(MechanicsDofSubMatrixStorageSources
    dofSubMatrixStorage/DofStatus.cpp
    dofSubMatrixStorage/BlockStorageBase.cpp
    dofSubMatrixStorage/BlockSparseMatrix.cpp
    dofSubMatrixStorage/BlockScalar.cpp
    dofSubMatrixStorage/BlockFullVector.cpp
    dofSubMatrixStorage/BlockFullMatrix.cpp
    )

list(APPEND MechanicsElementsSources
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

set(MechanicsGroupsSources
    groups/Group.cpp
    groups/GroupBase.cpp
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

set(MechanicsInterpolationtypesSources
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

set(MechanicsIGASources
    IGA/BSplineCurve.cpp
    IGA/BSplineSurface.cpp
    )

set(MechanicsLoadsSources
    loads/LoadNode.cpp
    loads/LoadNodeGroup.cpp
    loads/LoadNodeForces1D.cpp
    loads/LoadNodeForces2D.cpp
    loads/LoadNodeForces3D.cpp
    loads/LoadNodeScalarSource.cpp
    loads/LoadNodeGroupForces1D.cpp
    loads/LoadNodeGroupForces2D.cpp
    loads/LoadNodeGroupForces3D.cpp
    loads/LoadSurfaceBase2D.cpp
    loads/LoadSurfaceBase3D.cpp
    loads/LoadSurfaceConstDirection2D.cpp
    loads/LoadSurfaceConstDirection3D.cpp
    loads/LoadSurfacePressure2D.cpp
    loads/LoadSurfacePressure3D.cpp
    loads/LoadSurfacePressureFunction2D.cpp
    )

set(MechanicsMeshSources
    mesh/ConvertToInterpolationType.cpp
    mesh/ElementPrismCreate.cpp
    mesh/ImportFromGmsh.cpp
    mesh/MeshCompanion.cpp
    mesh/MeshGenerator.cpp
    )

set(MechanicsNodesSources
    nodes/NodeBase.cpp
    nodes/NodeDof.cpp
    nodes/NodeEnum.cpp
    )

set(MechanicsSectionSources
    sections/Section.cpp
    sections/SectionTruss.cpp
    sections/SectionVariableTruss.cpp
    sections/SectionPlane.cpp
    sections/SectionFibreMatrixBond.cpp
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
    timeIntegration/NewmarkDirect.cpp
    timeIntegration/NystroemBase.cpp
    timeIntegration/NystroemQinZhu.cpp
    timeIntegration/RungeKutta2.cpp
    timeIntegration/RungeKutta3.cpp
    timeIntegration/RungeKutta38.cpp
    timeIntegration/RungeKutta4.cpp
    timeIntegration/RungeKuttaCashKarp.cpp
    timeIntegration/RungeKuttaDormandPrince.cpp
    timeIntegration/RungeKuttaBase.cpp
    timeIntegration/TimeControl.cpp
    timeIntegration/TimeIntegrationBase.cpp
    timeIntegration/ImplicitExplicitBase.cpp
    timeIntegration/ImplEx.cpp
    timeIntegration/ImplExCallback.cpp
    timeIntegration/VelocityVerlet.cpp
    timeIntegration/postProcessing/PostProcessor.cpp
    timeIntegration/postProcessing/ResultBase.cpp
    timeIntegration/postProcessing/ResultElementIpData.cpp
    timeIntegration/postProcessing/ResultGroupNodeForce.cpp
    timeIntegration/postProcessing/ResultGroupNodeDof.cpp
    timeIntegration/postProcessing/ResultNodeAcceleration.cpp
    timeIntegration/postProcessing/ResultNodeDisp.cpp
    timeIntegration/postProcessing/ResultNodeDof.cpp
    timeIntegration/postProcessing/ResultTime.cpp
    )

set(MechanicsToolsSources
    tools/GlobalFractureEnergyIntegrator.cpp
    )

set(MechanicsSources
    ${MechanicsConstitutiveSources}
    ${MechanicsConstraintSources}
    ${MechanicsCrackSources}
    ${MechanicsDofSubMatrixStorageSources}
    ${MechanicsElementsSources}
    ${MechanicsGroupsSources}
    ${MechanicsIntegrationtypesSources}
    ${MechanicsInterpolationtypesSources}
    ${MechanicsIGASources}
    ${MechanicsLoadsSources}
    ${MechanicsMeshSources}
    ${MechanicsNodesSources}
    ${MechanicsSectionSources}
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
