unset(all_integration_tests CACHE)

add_subdirectory(geometryConcrete)
add_subdirectory(math)
add_subdirectory(mechanics)
add_subdirectory(metamodel)
add_subdirectory(optimize)
add_subdirectory(visualize)

add_custom_target(integrationtests
    COMMAND ctest -R integration --output-on-failure
    DEPENDS ${all_integration_tests}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
