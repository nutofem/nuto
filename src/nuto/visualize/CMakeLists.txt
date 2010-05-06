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
  VisualizeComponentBase.cpp
  VisualizeComponentDisplacement.cpp
  VisualizeComponentNonlocalWeight.cpp
  VisualizeComponentEngineeringStrain.cpp
  VisualizeComponentEngineeringStress.cpp
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
ADD_LIBRARY(NuToVisualize STATIC ${NuToVisualize_SOURCES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(NuToVisualize PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

# install shared library
INSTALL(TARGETS NuToVisualize ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})
