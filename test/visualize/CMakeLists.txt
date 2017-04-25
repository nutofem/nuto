
add_unit_test(UnstructuredGrid
    visualize/Cell.cpp
    visualize/Point.cpp
    visualize/XMLWriter.cpp
    )

target_link_libraries(UnstructuredGrid ${VTK_LIBRARIES})
