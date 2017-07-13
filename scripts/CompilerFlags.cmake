# use newer c++ standard, enable additional warnings, support for sse4
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -msse4 -pedantic -std=c++14 -fPIC")

# find openmp
if(ENABLE_OPENMP)
    message(STATUS "Checking for OpenMP support ...")
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
else()
    set(OPENMP_FOUND FALSE)
    find_package(Threads REQUIRED)
    #the package finds the THREADS_LIBRARY but does not add the pthread flag
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

# before adding warnings that complain about the swig code anyway, we save the
# flags as they are
set(PYTHON_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Clang specific options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # prevent warnings from boost and eigen
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --system-header-prefix=boost/ --system-header-prefix=eigen3/")
    # prevent "suggest braces around initialization of subobject" warning
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces -Wno-zero-length-array")
endif()

# GCC specific options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "5")
        # warns if override is missing
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
    endif()
endif()
