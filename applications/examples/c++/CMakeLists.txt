if(ENABLE_MPI)
    add_subdirectory(feti)
endif()


# copy meshes
foreach(mesh
        HeatedPlate
        CSDAMesh
        PGDPlateMeso2D
        )
    create_symlink(${mesh}.msh)
endforeach()

# math/base tests
foreach(test
        MortonOrder
        TestOperators
        )

    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(${test} Math Base ${LAPACK_LIBRARIES})
    add_test(examples::cpp::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    append_to_examples(${test})
endforeach()

foreach(test
    CSDAPrism3D
    ExplicitTimeIntegration
    HeatedPlate
    InsulatedPlate
    PGDPlateMeso2D
    SparseDirectSolverMUMPS
    Temperature1D
    Truss1D2N
    )

    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base Visualize
        Boost::unit_test_framework ${LAPACK_LIBRARIES})
    add_test(examples::cpp::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    append_to_examples(${test})
endforeach()

if(ENABLE_OPENMP AND OPENMP_FOUND)
    add_executable(OpenMP OpenMP.cpp)
    target_link_libraries(OpenMP Math)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/openMPplot.gnu
        ${CMAKE_CURRENT_BINARY_DIR}/openMPplot.gnu COPYONLY)
    append_to_examples(OpenMP)
endif()

if(ENABLE_MKL AND MKLSolver_FOUND)
    add_executable(SparseDirectSolverMKLPardiso
        SparseDirectSolverMKLPardiso.cpp)
    target_link_libraries(SparseDirectSolverMKLPardiso Math Base
        ${LAPACK_LIBRARIES} ${MKLSolver_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT})
    add_test(examples::cpp::SparseDirectSolverMKLPardiso
        ${CMAKE_CURRENT_BUILD_DIR}/SparseDirectSolverMKLPardiso)
    append_to_examples(SparseDirectSolverMKLPardiso)
endif()

add_executable(GeometryConcrete2D GeometryConcrete2D.cpp)
target_link_libraries(GeometryConcrete2D Mechanics Math Base Visualize
    GeometryConcrete)
add_test(examples::cpp::GeometryConcrete2D
    ${CMAKE_CURRENT_BUILD_DIR}/GeometryConcrete2D)
append_to_examples(GeometryConcrete2D)
