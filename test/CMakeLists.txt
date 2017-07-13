set(all_unit_tests "" CACHE INTERNAL "The names of all the unit tests")

add_subdirectory(base)
add_subdirectory(math)
if(ENABLE_MECHANICS)
    add_subdirectory(mechanics)
endif()

add_subdirectory(trilinos)

if(ENABLE_VISUALIZE)
    add_subdirectory(visualize)
endif()

add_custom_target(unit
    COMMAND ctest -R unit
    DEPENDS ${all_unit_tests}
    )
