if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(FATAL_ERROR "Coverage only available when compiling with GCC.")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
