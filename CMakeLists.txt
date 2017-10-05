cmake_minimum_required(VERSION 3.5)

project(NuTo)

option(BUILD_SHARED_LIBS "build NuTo libraries as shared libraries" TRUE)
option(ENABLE_COVERAGE "Set compiler flag for coverage analysis" FALSE)
option(ENABLE_EXAMPLES "enables examples" TRUE)
option(ENABLE_MKL "enables support for mkl blas, lapack, solvers" FALSE)
option(ENABLE_MUMPS "enables support for MUMPS solver" TRUE)
option(ENABLE_CUSTOM_EXAMPLES
    "enable private examples in separate directory (applications/custom)" FALSE)
option(ENABLE_OPENMP "enables OpenMP" FALSE)
option(ENABLE_PARDISO "enables support for PARDISO solver" FALSE)
option(ENABLE_PYTHON "create python wrapper for nuto" TRUE)
option(ENABLE_VISUALIZE "enables the visualization package" TRUE)
option(ENABLE_MPI "enables message passing interface" FALSE)

# set module path for custom cmake scripts
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts ${CMAKE_MODULE_PATH})
include(NuToMacros)
include(CompilerFlags)
include(GetGitRevisionDescription)
get_git_head_revision(git_refspec commit_nr)
configure_file(${CMAKE_SOURCE_DIR}/scripts/CTestCustom.cmake
    ${CMAKE_BINARY_DIR}/CTestCustom.cmake COPYONLY)

include_directories(SYSTEM external/)

include(CheckForDependencies)

if(ENABLE_VISUALIZE)
    add_definitions("-DENABLE_VISUALIZE")
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_VISUALIZE")
endif()

if(ENABLE_COVERAGE)
    include(Coverage)
endif()

if(ENABLE_PYTHON)
    add_subdirectory(python/nuto)
endif()

# build nuto
add_subdirectory(src)

#test Environment
include(CTest)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test/tools)
add_subdirectory(test)
add_subdirectory(applications/integrationtests)

if(ENABLE_EXAMPLES)
    add_subdirectory(applications/examples)
endif()

add_subdirectory(applications/benchmark)

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

install(FILES scripts/NuToMacros.cmake
    DESTINATION lib/cmake/NuTo
    )

install(EXPORT NuToTargets
    FILE NuToConfig.cmake
    NAMESPACE NuTo::
    DESTINATION lib/cmake/NuTo
    )

# This makes qtcreator show header-only classes...
file(GLOB_RECURSE HeaderFiles "*.h")
add_custom_target(headers SOURCES ${HeaderFiles})
