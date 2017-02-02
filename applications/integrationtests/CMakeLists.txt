if(ENABLE_PYTHON)
    foreach(test
        Brick8N
        CreateNodesAndGroupOperations
        GlobalDofs
        Sections
        Tetrahedron10N
        DeleteElemNode
        ElementUniaxialTest
        Quad2DPatchTest
        HeatConduction
        )

        add_test(p::mechanics::${test} ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py
            ${CMAKE_SYSTEM_NAME}
	    ${CMAKE_SYSTEM_PROCESSOR}
	    ${CMAKE_CURRENT_SOURCE_DIR}
	    ${CMAKE_CURRENT_BINARY_DIR})
        set_property(TEST p::mechanics::${test} PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    endforeach()
endif()


if(ENABLE_C++_TESTS)
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
        execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/meshFiles/${mesh}.msh" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}.msh")
    endforeach()

    # generate tests
    foreach(test
        BlockMatrices
        CoefficientChecks
        MoistureTransport
        ConstraintsNodeToElement
        CSDAInterface
        ElementUniaxialBrick3D
        ElementUniaxialPrism3D
        ElementUniaxialQuad2D
        ElementUniaxialTetrahedron3D
        ElementUniaxialTriangle2D
        ElementUniaxialTruss1D
        EngineeringInvariants
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
        add_test(c::mechanics::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    endforeach()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND ENABLE_MPI)

        # create symlinks to the necessary mesh files
        foreach(task
                _compare
                0
                1
                )
            execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/meshFiles/feti_rectangle_tri.msh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_rectangle_tri.msh${task}")
        endforeach()

        set(test Feti)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES} ${MPI_LIBRARIES})
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        if(MUMPS_FOUND)
            target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
        endif()
        set(MPIEXEC_NUMPROC_FLAG -np 2)
        set(MPIEXEC mpirun)
        add_test(NAME c::mechanics::${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    endif()


endif()
