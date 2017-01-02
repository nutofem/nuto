add_custom_target(benchmark COMMAND ctest -V -R benchmark)

# generate tests
foreach(test
    BuildGradient
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase -lprofiler)
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(benchmark::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
