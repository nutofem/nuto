set(NuToVisualize_SOURCES
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
    VisualizeDataBase.cpp
    VisualizeDataField.cpp
    VisualizeDataScalar.cpp
    VisualizeDataTensor.cpp
    VisualizeDataVector.cpp
    VisualizeDataType.cpp
    VisualizeUnstructuredGrid.cpp
)

# build shared library
add_library(NuToVisualize ${NuToVisualize_SOURCES})

# install library
install(TARGETS NuToVisualize
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
