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

add_executable(SparseMatrixTests SparseMatrixTests.cpp)
target_link_libraries(SparseMatrixTests NuToMath NuToBase)
add_test(integration::SparseMatrixTests.cpp ${CMAKE_CURRENT_BINARY_DIR}/SparseMatrixTests)
