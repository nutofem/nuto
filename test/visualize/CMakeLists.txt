
add_unit_test(VisualizeUnstructuredGrid
    visualize/CellBase.cpp
    visualize/Point.cpp
    )

target_link_libraries(VisualizeUnstructuredGrid ${VTK_LIBRARIES})
