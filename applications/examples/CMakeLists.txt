unset(all_examples CACHE)

add_subdirectory(c++)
add_subdirectory(python)

add_custom_target(examples
    COMMAND ctest -R example --output-on-failure
    DEPENDS ${all_examples}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
