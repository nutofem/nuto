add_subdirectory(paraviewScripts)

set(MESH_DIR
    "${CMAKE_SOURCE_DIR}/applications/examples/c++/feti/2d/phase_field/mesh")

foreach(task RANGE 0 1)
    foreach(mesh single_edge_notched_tension_test_2_subdomains_quad.mesh)
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
            "${MESH_DIR}/${mesh}${task}"
            "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}")
    endforeach()
endforeach()

foreach(file
        single_edge_notched_tension_test
        )
    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file}
        Mechanics Math Base FetiLib ${Boost_LIBRARIES} ${LAPACK_LIBRARIES}
        ${ANN_LIBRARIES} ${MPI_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${file} Visualize)
    endif()
    set(MPIEXEC_NUMPROC_FLAG -np 2)
    set(MPIEXEC mpirun)
    add_test(NAME examples::cpp::feti::${file}
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
        ${CMAKE_CURRENT_BINARY_DIR}/${file})
endforeach()

