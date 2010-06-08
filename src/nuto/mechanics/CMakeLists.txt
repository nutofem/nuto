SET(NuToMechanicsConstitutive_SOURCES
constitutive/ConstitutiveBase.cpp
constitutive/ConstitutiveStaticDataBase.cpp
constitutive/ConstitutiveTangentLocalBase.cpp
constitutive/ConstitutiveTangentLocal1x1.cpp
constitutive/ConstitutiveTangentLocal3x3.cpp
constitutive/ConstitutiveTangentLocal6x6.cpp
constitutive/mechanics/ConstitutiveEngineeringStressStrain.cpp
constitutive/mechanics/ConstitutiveMisesPlasticity.cpp
constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.cpp
constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity3D.cpp
constitutive/mechanics/ConstitutiveStaticDataMisesPlasticity3D.cpp
constitutive/mechanics/ConstitutiveStaticDataMisesPlasticityWithEnergy3D
constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.cpp
constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain3D.cpp
constitutive/mechanics/EngineeringStress1D.cpp
constitutive/mechanics/EngineeringStress2D.cpp
constitutive/mechanics/EngineeringStress3D.cpp
constitutive/mechanics/DeformationGradient1D.cpp
constitutive/mechanics/DeformationGradient2D.cpp
constitutive/mechanics/DeformationGradient3D.cpp
constitutive/mechanics/EngineeringStrain1D.cpp
constitutive/mechanics/EngineeringStrain2D.cpp
constitutive/mechanics/EngineeringStrain3D.cpp
constitutive/mechanics/GreenLagrangeStrain1D.cpp
constitutive/mechanics/GreenLagrangeStrain2D.cpp
constitutive/mechanics/GreenLagrangeStrain3D.cpp
constitutive/mechanics/InterfaceTractions1D.cpp
constitutive/mechanics/InterfaceTractions2D.cpp
constitutive/mechanics/InterfaceTractions3D.cpp
constitutive/mechanics/LinearElastic.cpp
constitutive/mechanics/NonlocalDamagePlasticity.cpp
constitutive/mechanics/RelativeInterfaceDisplacements1D.cpp
constitutive/mechanics/RelativeInterfaceDisplacements2D.cpp
constitutive/mechanics/RelativeInterfaceDisplacements3D.cpp
constitutive/mechanics/SecondPiolaKirchhoffStress1D.cpp
constitutive/mechanics/SecondPiolaKirchhoffStress2D.cpp
constitutive/mechanics/SecondPiolaKirchhoffStress3D.cpp
)

SET(NuToMechanicsConstraint_SOURCES
constraints/ConstraintBase.cpp
constraints/ConstraintEquation.cpp
constraints/ConstraintEquationTerm.cpp
constraints/ConstraintNode.cpp
constraints/ConstraintNodeDisplacements1D.cpp
constraints/ConstraintNodeDisplacements2D.cpp
constraints/ConstraintNodeDisplacements3D.cpp
constraints/ConstraintNodeGroup.cpp
constraints/ConstraintNodeGroupDisplacements1D.cpp
constraints/ConstraintNodeGroupDisplacements2D.cpp
constraints/ConstraintNodeGroupDisplacements3D.cpp
)

SET(NuToMechanicsElements_SOURCES
elements/Brick8N.cpp
elements/ElementBase.cpp
elements/ElementDataBase.cpp
elements/ElementDataConstitutiveBase.cpp
elements/ElementDataConstitutiveIp.cpp
elements/ElementDataConstitutiveIpNonlocal.cpp
elements/ElementDataIpBase.cpp
elements/ElementDataNonlocalBase.cpp
elements/ElementWithDataBase.cpp
elements/IpDataBase.cpp
elements/IpDataNonlocalBase.cpp
elements/IpDataStaticData.cpp
elements/IpDataStaticDataBase.cpp
elements/IpDataStaticDataNonlocal.cpp
elements/Plane.cpp
elements/Plane2D.cpp
elements/Plane2D4N.cpp
elements/Solid.cpp
elements/Tetrahedron10N.cpp
elements/Truss.cpp
elements/Truss1D.cpp
elements/Truss1D2N.cpp
elements/Truss1D3N.cpp
elements/Tetrahedron4N.cpp
elements/Voxel8N.cpp
)

SET(NuToMechanicsGroups_SOURCES
groups/Group.cpp
groups/GroupBase.cpp
)

SET(NuToMechanicsIntegrationtypes_SOURCES
integrationtypes/IntegrationType1D.cpp
integrationtypes/IntegrationType1D2NConstVariableIp.cpp
integrationtypes/IntegrationType1D2NGauss1Ip.cpp
integrationtypes/IntegrationType1D2NGauss2Ip.cpp
integrationtypes/IntegrationType2D.cpp
integrationtypes/IntegrationType2D4NGauss1Ip.cpp
integrationtypes/IntegrationType2D4NGauss4Ip.cpp
integrationtypes/IntegrationType3D.cpp
integrationtypes/IntegrationType3D4NGauss1Ip.cpp
integrationtypes/IntegrationType3D8NGauss1Ip.cpp
integrationtypes/IntegrationType3D8NGauss2x2x2Ip.cpp
integrationtypes/IntegrationTypeBase.cpp
)

SET(NuToMechanicsLoads_SOURCES
loads/LoadBase.cpp
loads/LoadNode.cpp
loads/LoadNodeGroup.cpp
loads/LoadNodeForces1D.cpp
loads/LoadNodeForces3D.cpp
loads/LoadNodeGroupForces1D.cpp
loads/LoadNodeGroupForces3D.cpp
)

SET(NuToMechanicsNodes_SOURCES
nodes/NodeBase.cpp
nodes/NodeCoordinates1D.cpp
nodes/NodeCoordinates2D.cpp
nodes/NodeCoordinates3D.cpp
nodes/NodeDisplacements1D.cpp
nodes/NodeDisplacements2D.cpp
nodes/NodeDisplacements3D.cpp
nodes/NodeGrid1D.cpp
nodes/NodeGrid2D.cpp
nodes/NodeGrid3D.cpp
nodes/NodeRotations2D.cpp
nodes/NodeRotations3D.cpp
nodes/NodeTemperature.cpp
nodes/NodeNonlocalDataBase.cpp
)

SET(NuToMechanicsSection_SOURCES
sections/SectionBase.cpp
sections/SectionTruss.cpp
sections/SectionPlane.cpp
)

SET(NuToMechanicsStructures_SOURCES
structures/StructureBase.cpp
structures/StructureBaseConstitutiveLaw.cpp
structures/StructureBaseConstraints.cpp
structures/StructureBaseElement.cpp
structures/StructureBaseGroup.cpp
structures/StructureBaseIntegrationTypes.cpp
structures/StructureBaseLoad.cpp
structures/StructureBaseNode.cpp
structures/StructureBaseSection.cpp
structures/unstructured/Structure.cpp
structures/unstructured/StructureElement.cpp
structures/unstructured/StructureNode.cpp
structures/grid/StructureGrid.cpp
structures/grid/StructureGridNode.cpp
structures/grid/StructureGridElement.cpp
)

SET(NuToMechanics_SOURCES
MechanicsException.cpp ${NuToMechanicsConstitutive_SOURCES} ${NuToMechanicsConstraint_SOURCES}
${NuToMechanicsElements_SOURCES} ${NuToMechanicsGroups_SOURCES} ${NuToMechanicsIntegrationtypes_SOURCES}
 ${NuToMechanicsLoads_SOURCES}
${NuToMechanicsNodes_SOURCES} ${NuToMechanicsSection_SOURCES} ${NuToMechanicsStructures_SOURCES}
)


SET(NuToMechanics_SWIG_MODULES
  ModulMechanics
)

# build static library
ADD_LIBRARY(NuToMechanics STATIC ${NuToMechanics_SOURCES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(NuToMechanics PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
# install static library
INSTALL(TARGETS NuToMechanics ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/mechanics)
FOREACH(NuToMechanics_SWIG_MODULE ${NuToMechanics_SWIG_MODULES})
  SET_SOURCE_FILES_PROPERTIES(${NuToMechanics_SWIG_MODULE}.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(${NuToMechanics_SWIG_MODULE}.i PROPERTIES SWIG_FLAGS "${NuTo_SWIG_FLAGS}")
  SWIG_ADD_MODULE(${NuToMechanics_SWIG_MODULE} python ${NuToMechanics_SWIG_MODULE}.i)
  SWIG_LINK_LIBRARIES(${NuToMechanics_SWIG_MODULE} NuToMechanics NuToBase NuToMath ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    SWIG_LINK_LIBRARIES(${NuToMechanics_SWIG_MODULE}  NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  # check for unresolved symbols
  IF( NOT CYGWIN AND NOT WIN32 )
    SET_TARGET_PROPERTIES(${SWIG_MODULE_${NuToMechanics_SWIG_MODULE}_REAL_NAME} PROPERTIES LINK_FLAGS -Wl,-z,defs)
  ENDIF( NOT CYGWIN AND NOT WIN32 )
  IF(MINGW)
    SET_TARGET_PROPERTIES(${SWIG_MODULE_${NuToMechanics_SWIG_MODULE}_REAL_NAME} PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  NUTO_INSTALL_SWIG_PYTHON_MODULE(${NuToMechanics_SWIG_MODULE} ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/mechanics)
ENDFOREACH(NuToMechanics_SWIG_MODULE)
