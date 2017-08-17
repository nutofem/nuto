foreach(test
    Dof
    Mesh
    Cache
    Assembly
    )

    add_executable(todo${test} todo${test}.cpp)
    target_link_libraries(todo${test} Base ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    #add_test(todo::${test} todo${test} --log_level=message)
endforeach()
