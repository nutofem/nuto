add_subdirectory(uniaxial)

if(ENABLE_MPI)
    add_subdirectory(feti)
endif()

if(ENABLE_PYTHON)
    foreach(test
        Brick8N
        DeleteElemNode
        ElectricWave
        HeatConduction
        MassMatrixIntegration
        NonlinearLoad
        Quad2DPatchTest
        Tetrahedron10N
        TimeIntegrationResultForce
        )
        add_test(integration::mechanics::${test}.py ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py)
        set_property(TEST integration::mechanics::${test}.py PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    endforeach()
endif()

create_symlink("meshes")

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
    MisesPlasticity
    MultipleConstitutiveLaws
    NewmarkPlane2D4N
    ParallelAssembly
    PiezoelectricLaw
    PlateWithHole
    RestartFiles
    Shrinkage
    SpringDamperCombination
    StructureNodeTest
    SurfaceLoad
    ThermoElasticity1D
    TrussIn2D
    TrussIn3D
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(integration::mechanics::${test}.cpp ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()



if(ENABLE_MPI)
    add_executable(MPI MPI.cpp)
    target_link_libraries(MPI ${MPI_LIBRARIES})
    set(MPIEXEC_NUMPROC_FLAG -np 4)
    set(MPIEXEC mpirun)
    #set(OVERSUBSCRIBE --oversubscribe)
    add_test(NAME integration::mechanics::MPI COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${CMAKE_CURRENT_BINARY_DIR}/MPI)
endif()
