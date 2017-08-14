set(VisualizeSources
    Cell.cpp
    ComponentName.cpp
    Point.cpp
    UnstructuredGrid.cpp
    XMLWriter.cpp 
    DataArray.cpp
)

sources_to_objects("${VisualizeSources}" VisualizeObjects)

add_library(Visualize ${VisualizeObjects})
set_target_properties(Visualize PROPERTIES
    OUTPUT_NAME NuToVisualize
    )
target_include_directories(Visualize PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/nuto>
    )

install(TARGETS Visualize EXPORT NuToTargets
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include/nuto
    )
