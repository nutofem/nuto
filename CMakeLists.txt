project(NuTo)

cmake_minimum_required(VERSION 3.0)

option(BUILD_SHARED_LIBS "build NuTo libraries as shared libraries" TRUE)
option(ENABLE_ARPACK
    "enables ARPACK libraries (for the solution of eigenvalue problems)" FALSE)
option(ENABLE_BENCHMARKS "enables Benchmarks" FALSE)
option(ENABLE_TESTS "create makefiles for c++ tests" TRUE)
option(ENABLE_COVERAGE "Set compiler flag for coverage analysis" FALSE)
option(ENABLE_EXAMPLES "enables examples" TRUE)
option(ENABLE_GEOMETRYCONCRETE "geometry for concrete mesoscale modeling" FALSE)
option(ENABLE_MECHANICS "enables the mechanics package" TRUE)
option(ENABLE_METAMODEL "enables the metamodel package" FALSE)
option(ENABLE_MKL "enables support for mkl blas, lapack, solvers" FALSE)
option(ENABLE_MUMPS "enables support for MUMPS solver" TRUE)
option(ENABLE_CUSTOM_EXAMPLES
    "enable private examples in separate directory (applications/custom)" FALSE)
option(ENABLE_OPENMP "enables OpenMP" FALSE)
option(ENABLE_OPTIMIZE "enables the optimize package" FALSE)
option(ENABLE_PARDISO "enables support for PARDISO solver" FALSE)
option(ENABLE_PYTHON "create python wrapper for nuto" TRUE)
option(ENABLE_STATIC_EXTERNALS
    "enables linking with static external libraries" FALSE)
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

if(ENABLE_MPI)
    message(STATUS "Find MPI")

    find_package(MPI REQUIRED)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})

    # These variables are set by FindMPI
    message(STATUS "MPI_CXX_FOUND:              ${MPI_CXX_FOUND}")
    message(STATUS "MPI_CXX_COMPILER:           ${MPI_CXX_COMPILER}")
    message(STATUS "MPI_CXX_COMPILER_FLAGS:     ${MPI_CXX_COMPILER_FLAGS}")
    message(STATUS "MPI_CXX_INCLUDE_PATH:       ${MPI_CXX_INCLUDE_PATH}")
    message(STATUS "MPI_CXX_LINK_FLAGS:         ${MPI_CXX_LINK_FLAGS}")
    message(STATUS "MPI_CXX_LIBRARIES:          ${MPI_CXX_LIBRARIES}")

    # Additionally, FindMPI sets the following variables for running
    # MPI programs from the command line:
    message(STATUS "MPIEXEC:                    ${MPIEXEC}")
    message(STATUS "MPIEXEC_NUMPROC_FLAG:       ${MPIEXEC_NUMPROC_FLAG}")
    message(STATUS "MPIEXEC_PREFLAG:            ${MPIEXEC_PREFLAG}")
    message(STATUS "MPIEXEC_POSTFLAG:           ${MPIEXEC_POSTFLAG}")
endif()

if(ENABLE_METAMODEL)
    set(ENABLE_OPTIMIZE TRUE)
endif()

if(ENABLE_MECHANICS)
    add_definitions("-DENABLE_MECHANICS")
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_MECHANICS")
endif()

if(ENABLE_OPTIMIZE)
    add_definitions("-DENABLE_OPTIMIZE")
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_OPTIMIZE")
endif()

if(ENABLE_VISUALIZE)
    add_definitions("-DENABLE_VISUALIZE")
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_VISUALIZE")
endif()

# search path for external libraries
# check if NUTO_EXTERNAL_LIBRARIES_ROOT is set
if(NOT NUTO_EXTERNAL_LIBRARIES_ROOT
        AND DEFINED $ENV{NUTO_EXTERNAL_LIBRARIES_ROOT})
    set(NUTO_EXTERNAL_LIBRARIES_ROOT $ENV{NUTO_EXTERNAL_LIBRARIES_ROOT})
endif()
# convert path to unix style path and set search path
if(NUTO_EXTERNAL_LIBRARIES_ROOT)
    file(TO_CMAKE_PATH ${NUTO_EXTERNAL_LIBRARIES_ROOT}
        NUTO_EXTERNAL_LIBRARIES_ROOT)
endif()

if(ENABLE_COVERAGE)
    include(Coverage)
endif()

# search for MKL
if(ENABLE_MKL)
    message(STATUS "Checking for MKL LAPACK and MKL SOLVER ...")
    find_package(MKLSolver)
    if(MKLSolver_FOUND)
        message(STATUS "MKLSolver_LIBRARIES = ${MKLSolver_LIBRARIES}")
        message(STATUS "MKLSolver_INCLUDE_DIR = ${MKLSolver_INCLUDE_DIR}")
        message(STATUS "MKLSolver_DEFINITIONS = ${MKLSolver_DEFINITIONS}")
        add_definitions(${MKLSolver_DEFINITIONS})
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};${MKLSolver_DEFINITIONS}")
        include_directories(${MKLSolver_INCLUDE_DIR})
    endif()
endif()

include(CheckBoost)

# find OpenBLAS
if(ENABLE_MUMPS OR ENABLE_PARDISO)
    if(NOT OpenBLAS_ROOT AND NOT DEFINED $ENV{OpenBLAS_ROOT})
        if(NUTO_EXTERNAL_LIBRARIES_ROOT)
            set(OpenBLAS_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
        endif()
    endif()
    set(OpenBLAS_FIND_SHARED_LIBRARY TRUE)
    find_package(OPENBLAS REQUIRED)
    if(NOT OpenBLAS_FOUND)
        message(STATUS "OpenBLAS not found")
    endif()
endif()

# find mumps solver
if(ENABLE_MUMPS)
    if(NOT MUMPS_ROOT AND NOT DEFINED $ENV{MUMPS_ROOT})
        if(NUTO_EXTERNAL_LIBRARIES_ROOT)
            set(MUMPS_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
        endif()
    endif()
    if(ENABLE_STATIC_EXTERNALS)
        set(MUMPS_FIND_STATIC_LIBRARY TRUE)
        set(METIS_FIND_STATIC_LIBRARY TRUE)
    endif()
    if(NOT METIS_ROOT AND NOT DEFINED $ENV{METIS_ROOT})
        if(NUTO_EXTERNAL_LIBRARIES_ROOT)
            set(METIS_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
        endif()
    endif()
    set(MUMPS_DEBUG TRUE)
    find_package(MUMPS REQUIRED)
    if(MUMPS_FOUND)
        add_definitions("-DHAVE_MUMPS")
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DHAVE_MUMPS")
        include_directories(${MUMPS_INCLUDE_DIR})
    endif()
else()
    set(MUMPS_FOUND FALSE)
    message(STATUS "MUMPS not switched on")
endif()

if(ENABLE_PYTHON)
    add_subdirectory(python/nuto)
    add_subdirectory(external/pybind11)
    add_subdirectory(pybind)
endif()

# find pardiso solver
if(ENABLE_PARDISO)
    if(NOT PARDISO_ROOT AND NOT DEFINED $ENV{PARDISO_ROOT})
        if(NUTO_EXTERNAL_LIBRARIES_ROOT)
            set(PARDISO_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
        endif()
    endif()
    set(PARDISO_FIND_SHARED_LIBRARY TRUE)
    set(PARDISO_DEBUG TRUE)
    find_package(PARDISO REQUIRED)
    if(PARDISO_FOUND)
        add_definitions("-DHAVE_PARDISO")
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DHAVE_PARDISO")
    endif()
else()
    set(PARDISO_FOUND FALSE)
    message(STATUS "PARDISO not switched on")
endif()

# find ANN library for nearest neighbor search
if(NOT ANN_ROOT AND NOT DEFINED $ENV{ANN_ROOT})
    if(NUTO_EXTERNAL_LIBRARIES_ROOT)
        set(ANN_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
    endif()
endif()
if(ENABLE_STATIC_EXTERNALS)
    set(ANN_FIND_STATIC_LIBRARY TRUE)
endif()
set(ANN_DEBUG TRUE)
find_package(ANN REQUIRED)
if(ANN_FOUND)
    include_directories(${ANN_INCLUDE_DIR})
endif()

# find Eigen header files (Linear Algebra)
if(NOT EIGEN_ROOT AND NOT DEFINED $ENV{EIGEN_ROOT})
    if(NUTO_EXTERNAL_LIBRARIES_ROOT)
        set(EIGEN_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
    endif()
endif()
find_package(EIGEN REQUIRED)
if(EIGEN_FOUND)
    message(STATUS "EIGEN_VERSION_NUMBER = ${EIGEN_VERSION_NUMBER}")
    # the math directory has to be included for the matrixbase addons
    if(EIGEN_VERSION_NUMBER VERSION_LESS "3.2.0")
        message(FATAL_ERROR "Insufficient eigen version (at least 3.2.0)")
    endif()

    include_directories(${EIGEN_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Fatal error: Eigen (version >= 3.0) required.\n")
endif()

# find ARPACK library for eigenvalue analysis
if(ENABLE_ARPACK)
    if(NOT ARPACK_ROOT AND NOT DEFINED $ENV{ARPACK_ROOT})
        if(NUTO_EXTERNAL_LIBRARIES_ROOT)
            set(ARPACK_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
        endif()
    endif()
    if(ENABLE_STATIC_EXTERNALS)
        set(ARPACK_FIND_STATIC_LIBRARY TRUE)
    endif()
    set(ARPACK_DEBUG TRUE)
    find_package(ARPACK REQUIRED)
    if(ARPACK_FOUND)
        add_definitions("-DHAVE_ARPACK")
    endif()
endif()

# find lapack library
if(ENABLE_MKL)
    find_package(LAPACK REQUIRED)
endif()

# build nuto
add_subdirectory(src)

#test Environment
include(CTest)

if(ENABLE_TESTS)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test/tools)
    add_subdirectory(test)
    add_subdirectory(applications/integrationtests)
endif()

if(ENABLE_EXAMPLES)
    add_subdirectory(applications/examples)
endif()

if(ENABLE_BENCHMARKS)
    add_subdirectory(applications/benchmark)
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
