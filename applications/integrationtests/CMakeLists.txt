add_subdirectory(geometryConcrete)

add_subdirectory(metamodel)

add_subdirectory(optimize)

add_subdirectory(mechanics)

if(ENABLE_VISUALIZE)
    add_subdirectory(visualize)
endif()


add_executable(SparseMatrixTests EXCLUDE_FROM_ALL SparseMatrixTests.cpp)
target_link_libraries(SparseMatrixTests Math Base)
add_test(integration::SparseMatrixTests.cpp SparseMatrixTests)
append_to_tests(SparseMatrixTests)

add_custom_target(integrationtests
    COMMAND ctest -R integration --output-on-failure
    DEPENDS ${all_integration_tests}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
