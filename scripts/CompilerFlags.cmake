# use newer c++ standard, enable additional warnings, support for sse4
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -msse4 -pedantic -std=c++14")
# before adding warnings that complain about the swig code anyway, we save the
# flags as they are
set(PYTHON_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Clang specific options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # prevent warnings from boost and eigen
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --system-header-prefix=boost/ --system-header-prefix=eigen3/")
    # prevent "suggest braces around initialization of subobject" warning
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
endif()

# GCC specific options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "5")
        # warns if override is missing
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
    endif()
endif()
