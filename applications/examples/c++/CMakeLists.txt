# TODO: these examples have been commented out at some point in time, we should
# try to resurrect some of them:
# - ConstraintsLagrange1D
# - NonlocalDamagePlasticity
# - NonlocalDamagePlasticityBeam
# - NonlocalDamagePlasticityModel
# - ImportGmsh

# copy meshes
foreach(mesh
        HeatedPlate
        CSDAMesh
        PGDPlateMeso2D
        )
    execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${mesh}.msh" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}.msh")
endforeach()




# math/base tests
foreach(test
        MortonOrder
        TestOperators
        )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(examples::cpp::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

if (ENABLE_MECHANICS)
foreach(test
    CSDAPrism3D
    ExplicitTimeIntegration
    HeatedPlate
    InsulatedPlate
    PGDPlateMeso2D
    Plane2D4NFreeHangingNode
    SparseDirectSolverMUMPS
    Temperature1D
    Truss1D2N
    VisualizePlane2D4N
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(examples::cpp::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
endif()

if(ENABLE_OPENMP AND OPENMP_FOUND)
    add_executable(OpenMP OpenMP.cpp)
    target_link_libraries(OpenMP NuToMath)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/openMPplot.gnu ${CMAKE_CURRENT_BINARY_DIR}/openMPplot.gnu COPYONLY)
endif()

if(ENABLE_MKL AND MKLSolver_FOUND)
    add_executable(SparseDirectSolverMKLPardiso SparseDirectSolverMKLPardiso.cpp)
    target_link_libraries(SparseDirectSolverMKLPardiso NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${MKLSolver_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    add_test(examples::cpp::SparseDirectSolverMKLPardiso ${CMAKE_CURRENT_BUILD_DIR}/SparseDirectSolverMKLPardiso)
endif()

if(ENABLE_GEOMETRYCONCRETE AND ENABLE_MECHANICS)
    add_executable(GeometryConcrete2D GeometryConcrete2D.cpp)
    target_link_libraries(GeometryConcrete2D NuToMechanics NuToMath NuToBase NuToGeometryConcrete ${Boost_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(GeometryConcrete2D NuToVisualize)
    endif()
    add_test(examples::cpp::GeometryConcrete2D ${CMAKE_CURRENT_BUILD_DIR}/GeometryConcrete2D)
endif()
