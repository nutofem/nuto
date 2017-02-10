add_subdirectory(uniaxial)

if(ENABLE_PYTHON)
    foreach(test
        Brick8N
        DeleteElemNode
        GlobalDofs
        HeatConduction
        Quad2DPatchTest
        Tetrahedron10N
        )

        add_test(integration::mechanics::${test} ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py
            ${CMAKE_SYSTEM_NAME}
	    ${CMAKE_SYSTEM_PROCESSOR}
	    ${CMAKE_CURRENT_SOURCE_DIR}
	    ${CMAKE_CURRENT_BINARY_DIR})
        set_property(TEST integration::mechanics::${test} PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    endforeach()
endif()

# copy meshes
foreach(mesh
    2D
    PlateWithHole
    InterfaceElements
    ConstraintsNodeToElementMatrixMesh2d
    ConstraintsNodeToElementFiberMesh2d
    ConstraintsNodeToElementMatrixMesh3d
    ConstraintsNodeToElementFiberMesh3d
    TwoElements
    )
    execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/meshes/${mesh}.msh" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}.msh")
endforeach()

# generate tests
foreach(test
    BlockMatrices
    CoefficientChecks
    MoistureTransport
    ConstraintsNodeToElement
    CSDAInterface
    GradientDamage
    IGA
    ImplEx
    InterfaceElements
    IntegrationPointVoronoiCells
    InterpolationTypes
    LobattoIntegration
    MisesPlasticity
    MultipleConstitutiveLaws
    NewmarkPlane2D4N
    ParallelAssembly
    RungeKutta
    Shrinkage
    SpringDamperCombination
    SurfaceLoad
    ThermoMechanics2D
    TrussIn2D
    TrussIn3D
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(integration::mechanics::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()

#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND ENABLE_MPI)
#
#    # create symlinks to the necessary mesh files
#    foreach(task
#            _compare
#            0
#            1
#            )
#        execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/meshes/feti_rectangle_tri.msh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_rectangle_tri.msh${task}")
#    endforeach()
#
#    add_executable(Feti Feti.cpp)
#    target_link_libraries(Feti NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
#    if(ENABLE_VISUALIZE)
#        target_link_libraries(Feti NuToVisualize)
#    endif()
#    if(MUMPS_FOUND)
#        target_link_libraries(Feti ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#    endif()
#    set(MPIEXEC_NUMPROC_FLAG -np 2)
#    set(MPIEXEC mpirun)
#    add_test(NAME integration::mechanics::Feti COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${CMAKE_CURRENT_BINARY_DIR}/Feti)
#endif()


if(ENABLE_MPI)
    add_executable(MPI MPI.cpp)
    target_link_libraries(MPI ${MPI_LIBRARIES})
    set(MPIEXEC_NUMPROC_FLAG -np 8)
    set(MPIEXEC mpirun)
    set(OVERSUBSCRIBE --oversubscribe)
    add_test(NAME integration::mechanics::MPI COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${OVERSUBSCRIBE} ${CMAKE_CURRENT_BINARY_DIR}/MPI)
endif()