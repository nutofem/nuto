cmake_minimum_required(VERSION 3.5)

project(NuTo)

option(BUILD_SHARED_LIBS "build NuTo libraries as shared libraries" TRUE)

# set module path for custom cmake scripts
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake ${CMAKE_MODULE_PATH})
include(NuToMacros)

set_compiler_flags()
get_git_head_revision(git_refspec commit_nr)
check_for_dependencies()

add_subdirectory(external)

# build nuto
add_subdirectory(nuto)

# test Environment
include(CTest)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(applications/integrationtests EXCLUDE_FROM_ALL)
add_subdirectory(applications/examples EXCLUDE_FROM_ALL)
if(benchmark_FOUND)
    add_subdirectory(applications/benchmark EXCLUDE_FROM_ALL)
endif()

# Copy headers
install(DIRECTORY nuto/
    DESTINATION include/nuto
    FILES_MATCHING PATTERN "*.h"
    )

# Copy cmake macros
install(DIRECTORY scripts/cmake/
    DESTINATION lib/cmake/NuTo
    FILES_MATCHING PATTERN "*.cmake"
    )

install(EXPORT NuToTargets
    FILE NuToTargets.cmake
    NAMESPACE NuTo::
    DESTINATION lib/cmake/NuTo
    )
