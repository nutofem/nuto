cmake_minimum_required(VERSION 3.0)

project(NuTo)

option(BUILD_SHARED_LIBS "build NuTo libraries as shared libraries" TRUE)
option(ENABLE_MKL "enables support for mkl blas, lapack, solvers" FALSE)
option(ENABLE_MUMPS "enables support for MUMPS solver" TRUE)
option(ENABLE_CUSTOM_EXAMPLES "enable private examples in separate directory (applications/custom)" TRUE)
option(ENABLE_OPENMP "enables OpenMP" TRUE)
option(ENABLE_PARDISO "enables support for PARDISO solver" FALSE)
option(ENABLE_PYTHON "create python wrapper for nuto" TRUE)
option(ENABLE_MPI "enables message passing interface" TRUE)


# set module path for custom cmake scripts
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake ${CMAKE_MODULE_PATH})
include(NuToMacros)

set_compiler_flags()
get_git_head_revision(git_refspec commit_nr)
check_for_dependencies()

add_subdirectory(external)
include_directories(external/)

if(ENABLE_PYTHON)
    add_subdirectory(python/nuto)
endif()

# build nuto
add_subdirectory(src)

#test Environment
include(CTest)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(applications/integrationtests EXCLUDE_FROM_ALL)
add_subdirectory(applications/examples EXCLUDE_FROM_ALL)
if (benchmark_FOUND)
    add_subdirectory(applications/benchmakr EXCLUDE_FROM_ALL)
endif()

# Custom examples are applications you don't want to add to the repository, but
# that should still be picked up by CMake. If this directory is empty, nothing
# happens.
if(ENABLE_CUSTOM_EXAMPLES)
    add_subdirectory(applications/custom)
endif()

install(DIRECTORY src/
    DESTINATION include/nuto
    FILES_MATCHING PATTERN "*.h"
    )

install(DIRECTORY scripts/cmake/
    DESTINATION lib/cmake/NuTo
    FILES_MATCHING PATTERN "*.cmake")

install(EXPORT NuToTargets
    FILE NuToConfig.cmake
    NAMESPACE NuTo::
    DESTINATION lib/cmake/NuTo
    )

# This makes qtcreator show header-only classes...
file(GLOB_RECURSE HeaderFiles "*.h")
add_custom_target(headers SOURCES ${HeaderFiles})
