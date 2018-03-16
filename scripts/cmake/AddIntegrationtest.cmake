function(append_to_tests TestName)
    set(all_integration_tests "${all_integration_tests};${TestName}"
        CACHE INTERNAL "The names of all the integration tests")
endfunction()

function(add_integrationtest test)
    if (ARGN) # additional arguments to this function are appended to the test
              # command and represent additional args passed to the execuable
        foreach (arg ${ARGN})
            set(AdditionalArgs "${AdditionalArgs};${arg}")
        endforeach()
    endif()

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base Visualize 
        Boost::unit_test_framework ${LAPACK_LIBRARIES})
    target_include_directories(${test} PRIVATE ${CMAKE_SOURCE_DIR}/test/tools)
    get_filename_component(module ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    add_test(integration::${module}::${test}.cpp ${test} ${AdditionalArgs})
    append_to_tests(${test})
endfunction()
