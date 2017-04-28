
add_unit_test(UnstructuredGrid
    visualize/Cell.cpp
    visualize/Point.cpp
    visualize/XMLWriter.cpp
    visualize/base64.cpp
    )

if (ZLIB_FOUND)
    target_link_libraries(UnstructuredGrid ${ZLIB_LIBRARY})
endif()
