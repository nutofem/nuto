set(NuToVisualizeSources
    Cell.cpp
    ComponentName.cpp
    Point.cpp
    UnstructuredGrid.cpp
    XMLWriter.cpp 
    DataArray.cpp
)

sources_to_objects("${NuToVisualizeSources}" NuToVisualizeObjects)

add_library(NuToVisualize ${NuToVisualizeObjects})

install(TARGETS NuToVisualize
    LIBRARY DESTINATION lib
    )
