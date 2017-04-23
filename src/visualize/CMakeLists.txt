set(NuToVisualizeSources
    VisualizeUnstructuredGrid.cpp
    CellBase.cpp
    Point.cpp
)

sources_to_objects("${NuToVisualizeSources}" NuToVisualizeObjects)

# build shared library
add_library(NuToVisualize ${NuToVisualizeObjects})

# install library
install(TARGETS NuToVisualize
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
