if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(FATAL_ERROR "Coverage only available when compiling with GCC.")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
add_custom_target(coverage
    COMMAND ${CMAKE_MAKE_PROGRAM} -j8
    COMMAND ctest
    COMMAND lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file coverage.info
    COMMAND lcov -r coverage.info '/usr/include/*' -o coverage.info
    COMMAND genhtml coverage.info --output-directory coverageOutput
    COMMAND xdg-open coverageOutput/index.html)
