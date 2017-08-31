if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(geometryConcrete)
endif()

if(ENABLE_METAMODEL)
    add_subdirectory(metamodel)
endif()

if(ENABLE_OPTIMIZE)
    add_subdirectory(optimize)
endif()

if(ENABLE_MECHANICS)
    add_subdirectory(mechanics)
endif()

if(ENABLE_VISUALIZE)
    add_subdirectory(visualize)
endif()


add_executable(SparseMatrixTests SparseMatrixTests.cpp)
target_link_libraries(SparseMatrixTests Math Base)
add_test(integration::SparseMatrixTests.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/SparseMatrixTests)

add_custom_target(integrationtests COMMAND ctest -R integration)
