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
        add_test(integration::mechanics::${test}.py
            ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py)
        set_property(TEST integration::mechanics::${test}.py
            PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    endforeach()
endif()

create_symlink("meshes")

# generate tests
add_integrationtest(AdditiveInput)
add_integrationtest(AdditiveOutput)
add_integrationtest(BlockMatrices)
add_integrationtest(CoefficientChecks)
add_integrationtest(MoistureTransport)
add_integrationtest(ConstraintsNodeToElement)
add_integrationtest(ConstraintNodeToElement2D)
add_integrationtest(ConstraintNodeToElement3D)
add_integrationtest(ContinuumElement)
add_integrationtest(CreepUniaxial)
add_integrationtest(CSDAInterface)
add_integrationtest(GradientDamage)
add_integrationtest(IGA)
add_integrationtest(ImplEx)
add_integrationtest(InhomogeneousMaterial)
add_integrationtest(InterfaceElements)
add_integrationtest(IntegrationPointVoronoiCells)
add_integrationtest(InterpolationTypes)
add_integrationtest(MeshCompanion)
add_integrationtest(MisesPlasticity)
add_integrationtest(MultipleConstitutiveLaws)
add_integrationtest(NewmarkPlane2D4N)
add_integrationtest(ParallelAssembly)
add_integrationtest(PiezoelectricLaw)
add_integrationtest(PlateWithHole)
add_integrationtest(RestartFiles)
add_integrationtest(Shrinkage)
add_integrationtest(SpringDamperCombination)
add_integrationtest(StructureNodeTest)
add_integrationtest(SurfaceLoad)
add_integrationtest(ThermoElasticity1D)
add_integrationtest(TrussIn2D)
add_integrationtest(TrussIn3D)



if(ENABLE_MPI)
    add_executable(MPI EXCLUDE_FROM_ALL MPI.cpp)
    target_link_libraries(MPI ${MPI_LIBRARIES})
    set(MPIEXEC_NUMPROC_FLAG -np 4)
    set(MPIEXEC mpirun)
    #set(OVERSUBSCRIBE --oversubscribe)
    add_test(NAME integration::mechanics::MPI
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
        ${CMAKE_CURRENT_BINARY_DIR}/MPI)
    append_to_tests(MPI)
endif()
