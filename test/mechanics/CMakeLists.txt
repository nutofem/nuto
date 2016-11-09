# TODO: these test have been commented at some point in time, we should try
# to resurrect some of them:
# Python:
#   NonlocalDamagePlasticity
#   ConstraintsLagrange1D
# C++:
#   ParallelAssembly
#   Brick8NcycleJump
#   CycleJumpNonlocal
#   DogBoneDirectNewmark
#   NonlocalDamagePlasticityModel
#   StructureEvaluate
add_subdirectory(constitutive)

#test files for mechanics module
if(ENABLE_PYTHON_TESTS)
    foreach(test
        Brick8N
        ConstitutiveLaws
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
        set_property(TEST p::mechanics::${test} PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
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
        execute_process( COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${mesh}.msh" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}.msh")
    endforeach()

    # generate tests
    foreach(test
        GradientDamage
        NewmarkPlane2D4N
        RungeKutta
        ImplEx
        EngineeringInvariants
        SurfaceLoad
        CoefficientChecks
        ElementUniaxialTruss1D
        ElementUniaxialTriangle2D
        ElementUniaxialQuad2D
        ElementUniaxialBrick3D
        ElementUniaxialTetrahedron3D
        ThermoMechanics2D
        LobattoIntegration
        IPData
        IGA
        InterpolationTypes
        TrussIn2D
        TrussIn3D
        InterfaceElements
        BlockMatrices
        ConstraintsNodeToElement
        IntegrationPointVoronoiCells
        ConstitutiveLawMoistureTransport
        MultipleConstitutiveLaws
        Shrinkage
        SpringDamperCombination
        CSDAInterface
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

    if(NOT DEFINED ENV{TRAVIS})
        add_executable(MisesPlasticity MisesPlasticity.cpp)
        target_link_libraries(MisesPlasticity NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES})
        if(ENABLE_VISUALIZE)
          target_link_libraries(MisesPlasticity NuToVisualize)
        endif()
        if(MUMPS_FOUND)
            target_link_libraries(MisesPlasticity ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
        endif()
        add_test(c::mechanics::MisesPlasticity ${CMAKE_CURRENT_BINARY_DIR}/MisesPlasticity)
    endif()


    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

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





if(ENABLE_UNIT_TESTS)
    # TESTS WITH NuToMechancis
    foreach(test ContinuumElement PhaseField RestartFiles)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        add_test(c::unit::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    endforeach()

    # TESTS WITHOUT NuToMechanics
    foreach(test GlobalFractureEnergyIntegration)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMath NuToBase)
        add_test(c::unit::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    endforeach()


endif()
