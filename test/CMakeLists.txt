set(all_unit_tests "" CACHE INTERNAL "The names of all the unit tests")

add_subdirectory(base)
add_subdirectory(math)
add_subdirectory(mechanics)
add_subdirectory(visualize)

add_custom_target(unit
    COMMAND ctest -R unit --output-on-failure
    DEPENDS ${all_unit_tests}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
