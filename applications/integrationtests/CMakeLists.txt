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
add_test(c::mechanics::SparseMatrixTests ${CMAKE_CURRENT_BINARY_DIR}/SparseMatrixTests)

add_test(p::mechanics::SparseMatrixCSR ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/SparseMatrixCSR.py
    ${CMAKE_SYSTEM_NAME}
${CMAKE_SYSTEM_PROCESSOR}
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST p::mechanics::SparseMatrixCSR PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
