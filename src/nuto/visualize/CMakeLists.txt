# $Id$
SET(NuToVisualize_SOURCES
  CellBase.cpp
  CellHexahedron.cpp
  CellVertex.cpp
  CellLine.cpp
  CellPyramid.cpp
  CellQuad.cpp
  CellTetra.cpp
  CellTriangle.cpp
  Point.cpp
  VisualizeBase.cpp
  VisualizeComponent.cpp
  VisualizeComponentBase.cpp
  VisualizeComponentNonlocalWeight.cpp
  VisualizeDataBase.cpp
  VisualizeDataField.cpp
  VisualizeDataScalar.cpp
  VisualizeDataTensor.cpp
  VisualizeDataVector.cpp
  VisualizeDataType.cpp
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
