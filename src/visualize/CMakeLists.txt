set(VisualizeSources
    Cell.cpp
    ComponentName.cpp
    Point.cpp
    UnstructuredGrid.cpp
    XMLWriter.cpp
    DataArray.cpp
    QuadAverageHandler.cpp
    TensorProductVoronoiHandler.cpp
)

create_nuto_module(Visualize "${VisualizeSources}")
target_link_libraries(Visualize Math)
