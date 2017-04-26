set(NuToVisualizeSources
    Cell.cpp
    ComponentName.cpp
    Point.cpp
    UnstructuredGrid.cpp
    XMLWriter.cpp 
)

sources_to_objects("${NuToVisualizeSources}" NuToVisualizeObjects)

# build shared library
add_library(NuToVisualize ${NuToVisualizeObjects})
target_link_libraries(NuToVisualize ${VTK_LIBRARIES})

# install library
install(TARGETS NuToVisualize
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
