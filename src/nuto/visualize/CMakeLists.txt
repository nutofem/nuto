# $Id$
SET(NuToVisualize_SOURCES
  CellBase.cpp
  CellHexahedron.cpp
  CellLine.cpp
  CellPyramid.cpp
  CellQuad.cpp
  CellTetra.cpp
  CellTriangle.cpp
  Point.cpp
  VisualizeBase.cpp
  VisualizeComponentAcceleration.cpp
  VisualizeComponentAngularAcceleration.cpp
  VisualizeComponentAngularVelocity.cpp
  VisualizeComponentBase.cpp
  VisualizeComponentConstitutive.cpp
  VisualizeComponentDamage.cpp
  VisualizeComponentDisplacement.cpp
  VisualizeComponentCrack.cpp
  VisualizeComponentElement.cpp
  VisualizeComponentEngineeringPlasticStrain.cpp
  VisualizeComponentEngineeringStrain.cpp
  VisualizeComponentEngineeringStress.cpp
  VisualizeComponentHeatFlux.cpp
  VisualizeComponentLatticeStrain.cpp
  VisualizeComponentLatticeStress.cpp
  VisualizeComponentNonlocalWeight.cpp
  VisualizeComponentParticleRadius.cpp
  VisualizeComponentPrincipalEngineeringStress.cpp
  VisualizeComponentRotation.cpp
  VisualizeComponentSection.cpp
  VisualizeComponentTemperature.cpp
  VisualizeComponentVelocity.cpp
  VisualizeDataBase.cpp
  VisualizeDataField.cpp
  VisualizeDataScalar.cpp
  VisualizeDataTensor.cpp
  VisualizeDataVector.cpp
  VisualizeDataType.cpp
  VisualizeException.cpp
  VisualizeUnstructuredGrid.cpp
)

# build shared library
ADD_LIBRARY(NuToVisualize ${NuToVisualize_SOURCES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(NuToVisualize PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

# install library
INSTALL(TARGETS NuToVisualize
	ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
	LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})
