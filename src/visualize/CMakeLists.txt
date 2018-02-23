set(VisualizeSources
    Cell.cpp
    Point.cpp
    UnstructuredGrid.cpp
    XMLWriter.cpp
    DataArray.cpp
    AverageHandler.cpp
    Visualizer.cpp
    PointHandler.cpp
    PostProcess.cpp
    VoronoiHandler.cpp
    VoronoiGeometries.cpp
)

create_nuto_module(Visualize "${VisualizeSources}")
target_link_libraries(Visualize PRIVATE Math Boost::filesystem)
