project(NuTo CXX C Fortran) # C and Fortran are required to find LAPACK
cmake_minimum_required(VERSION 3.0)

option(BUILD_SHARED_LIBS "build NuTo libraries as shared libraries" TRUE)
option(ENABLE_ARPACK "enables ARPACK libraries (for the solution of eigenvalue problems)" FALSE)
option(ENABLE_C++_TESTS "create makefiles for c++ tests" TRUE)
option(ENABLE_CGAL "enables CGAL libraries (for cracks only)" FALSE)
option(ENABLE_COVERAGE "Set compiler flag for coverage analysis" FALSE)
option(ENABLE_DEBUG "enables debug information" TRUE)
option(ENABLE_DEBUG_SERIALIZATION "prints additional debug information during serialization" FALSE)
option(ENABLE_EXAMPLES "enables examples" FALSE)
option(ENABLE_GEOMETRYCONCRETE "geometry for concrete mesoscale modeling" FALSE)
option(ENABLE_GUI "enable building the GUI application" FALSE)
option(ENABLE_MECHANICS "enables the mechanics package" TRUE)
option(ENABLE_METAMODEL "enables the metamodel package" FALSE)
option(ENABLE_MKL "enables support for mkl blas, lapack, solvers" FALSE)
option(ENABLE_MUMPS "enables support for MUMPS solver" TRUE)
option(ENABLE_MY_EXAMPLES "enable building private examples in a separate directory (nuto/myNutoExamples)" FALSE)
option(ENABLE_NUMPY "enables numpy" TRUE)
option(ENABLE_OPENACC "enables OpenACC" FALSE)
option(ENABLE_OPENMP "enables OpenMP" FALSE)
option(ENABLE_OPTIMIZE "enables the optimize package" FALSE)
option(ENABLE_PARDISO "enables support for PARDISO solver" FALSE)
option(ENABLE_PYTHON_TESTS "create makefiles for python tests" TRUE)
option(ENABLE_SHOW_TIME "shows for each python command the execution time" TRUE)
option(ENABLE_STATIC_EXTERNALS "enables linking with static external libraries" FALSE)
option(ENABLE_SERIALIZATION "enables serialization" FALSE)
option(ENABLE_UNIT_TESTS "use Boost Test library for unit tests" TRUE)
option(ENABLE_VERBOSELEVEL "enables status messages during the execution by setting the verboselevel" FALSE)
option(ENABLE_VISUALIZE "enables the visualization package" FALSE)
option(ENABLE_MPI "enables message passing interface" FALSE)

# set module path for custom cmake scripts
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeScripts ${CMAKE_MODULE_PATH})
include(NuToMacros)

add_subdirectory(external)
include_directories(external/jsoncpp/dist/)

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    message("MPI_INCLUDE_PATH:  " ${MPI_INCLUDE_PATH})
    message("MPI_COMPILE_FLAGS: " ${MPI_COMPILE_FLAGS})
    message("MPI_LINK_FLAGS:    " ${MPI_LINK_FLAGS})
    message("MPI_LIBRARIES:     " ${MPI_LIBRARIES})
endif()

if(ENABLE_METAMODEL)
    set(ENABLE_OPTIMIZE TRUE)
endif()

if(ENABLE_MECHANICS)
    add_definitions("-DENABLE_MECHANICS")
    set(NuTo_DOXYGEN_PREDEFINED "ENABLE_MECHANICS" ${NuTo_DOXYGEN_PREDEFINED})
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_MECHANICS")
endif()

if(ENABLE_OPTIMIZE)
    add_definitions("-DENABLE_OPTIMIZE")
    set(NuTo_DOXYGEN_PREDEFINED "ENABLE_OPTIMIZE" ${NuTo_DOXYGEN_PREDEFINED})
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_OPTIMIZE")
endif()

if(ENABLE_VISUALIZE)
    add_definitions("-DENABLE_VISUALIZE")
    set(NuTo_DOXYGEN_PREDEFINED "ENABLE_VISUALIZE" ${NuTo_DOXYGEN_PREDEFINED})
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_VISUALIZE")
endif()

if(ENABLE_SERIALIZATION)
    add_definitions("-DENABLE_SERIALIZATION")
    set(NuTo_DOXYGEN_PREDEFINED "ENABLE_SERIALIZATION" ${NuTo_DOXYGEN_PREDEFINED})
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_SERIALIZATION")
endif()

if(ENABLE_VERBOSELEVEL)
    add_definitions("-DENABLE_VERBOSELEVEL")
    set(NuTo_DOXYGEN_PREDEFINED "ENABLE_VERBOSELEVEL" ${NuTo_DOXYGEN_PREDEFINED})
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_VERBOSELEVEL")
endif()

if(ENABLE_SHOW_TIME)
    add_definitions("-DSHOW_TIME")
    set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DSHOW_TIME")
endif()

# search path for external libraries
# check if NUTO_EXTERNAL_LIBRARIES_ROOT is set
if(NOT NUTO_EXTERNAL_LIBRARIES_ROOT AND DEFINED $ENV{NUTO_EXTERNAL_LIBRARIES_ROOT})
    set(NUTO_EXTERNAL_LIBRARIES_ROOT $ENV{NUTO_EXTERNAL_LIBRARIES_ROOT})
endif()
# convert path to unix style path and set search path
if(NUTO_EXTERNAL_LIBRARIES_ROOT)
    file(TO_CMAKE_PATH ${NUTO_EXTERNAL_LIBRARIES_ROOT} NUTO_EXTERNAL_LIBRARIES_ROOT)
endif()

# set build type
if(ENABLE_DEBUG)
    add_definitions("-DDEBUG")
    set(CMAKE_BUILD_TYPE "Debug")
else()
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release")
    endif()
endif()
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

if(ENABLE_DEBUG_SERIALIZATION)
    add_definitions("-DDEBUG_SERIALIZATION")
endif()

if(ENABLE_OPENACC)
    add_definitions("-D_OPENACC")
endif()

# find openmp
if(ENABLE_OPENMP)
    message(STATUS "Checking for OpenMP support ...")
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        endif()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
else()
    set(OPENMP_FOUND FALSE)
    find_package(Threads REQUIRED)
    #the package finds the THREADS_LIBRARY but does not add the pthread flag
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

# use newer c++ standard, enable additional warnings, support for sse4
# set(CMAKE_CXX_STANDARD 14) # use this if cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -msse4 -pedantic -std=c++14")
# prevent warnings from boost and eigen
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --system-header-prefix=boost/ --system-header-prefix=eigen3/ -Wno-inconsistent-missing-override -Wno-missing-braces")
endif()

if(ENABLE_COVERAGE)
    include(Coverage)
endif()

if(ENABLE_NUMPY)
    find_package(Numpy REQUIRED)
    if(NUMPY_FOUND)
        add_definitions("-DENABLE_NUMPY")
        include_directories(${NUMPY_INCLUDE_PATH})
        set(NuTo_DOXYGEN_PREDEFINED "ENABLE_NUMPY" ${NuTo_DOXYGEN_PREDEFINED})
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DENABLE_NUMPY")
        message(STATUS "Numpy found")
    else()
        message(STATUS "Numpy not found")
    endif()
else()
    message(STATUS "Numpy not activated")
endif()

# Python required for interface files c++ and python
message(STATUS "Checking for Python...")
find_package(PythonInterp 2 REQUIRED)
find_package(PythonLibs 2 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# get python version
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;print(sys.version[:3]);"
        OUTPUT_VARIABLE PYTHON_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")

# set relative paths for install
set(LIB_SUFFIX "" CACHE STRING "Define suffix of directory name (32/64)")
set(NUTO_LIBRARIES_INSTALL_PATH "lib${LIB_SUFFIX}")
set(NUTO_INCLUDE_INSTALL_PATH "include")
set(NUTO_PYTHON_MODULES_INSTALL_PATH "lib${LIB_SUFFIX}/python${PYTHON_VERSION}/site-packages")

# SWIG required for interface files c++ and python
message(STATUS "Checking for SWIG...")
find_program(SWIG_EXECUTABLE NAMES swig3.0 swig swig2.0)
message(STATUS "MYSWIG_EXECUTABLE = ${SWIG_EXECUTABLE1}")
find_package(SWIG REQUIRED)
message(STATUS "SWIG found: ${SWIG_EXECUTABLE}")
message(STATUS "SWIG_VERSION = ${SWIG_VERSION}")
if(SWIG_VERSION VERSION_LESS "3.0")
  message(FATAL_ERROR "Insufficient swig version (at least 3.0)")
endif()
include(${SWIG_USE_FILE})
set(NuTo_SWIG_FLAGS "-Wall;-fvirtual;${NuTo_SWIG_FLAGS}")

# search for MKL
if(ENABLE_MKL)
    message(STATUS "Checking for MKL LAPACK and MKL SOLVER ...")
    find_package(MKLSolver)
    if(MKLSolver_FOUND)
        message(STATUS "MKLSolver_LIBRARIES = ${MKLSolver_LIBRARIES}")
        message(STATUS "MKLSolver_INCLUDE_DIR = ${MKLSolver_INCLUDE_DIR}")
        message(STATUS "MKLSolver_DEFINITIONS = ${MKLSolver_DEFINITIONS}")
        add_definitions(${MKLSolver_DEFINITIONS})
        set(NuTo_DOXYGEN_PREDEFINED "${MKLSolver_DOXYGEN_DEFINITIONS}" ${NuTo_DOXYGEN_PREDEFINED})
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};${MKLSolver_DEFINITIONS}")
        include_directories(${MKLSolver_INCLUDE_DIR})
    endif()
endif()

# boost
message(STATUS "Checking for Boost...")
set(Boost_ADDITIONAL_VERSIONS "1.61;1.60;1.59;1.57;1.56;1.55;1.41.0;1.41;1.42;1.43;1.44;1.45;1.46;1.46.1;1.47;1.50;1.51;1.54")
if(ENABLE_STATIC_EXTERNALS)
    set(Boost_USE_STATIC_LIBS TRUE)
endif()
set(Boost_DETAILED_FAILURE_MSG TRUE)
if(NOT BOOST_ROOT AND NOT DEFINED $ENV{BOOST_ROOT})
    if(NUTO_EXTERNAL_LIBRARIES_ROOT)
        set(BOOST_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
    endif()
endif()

# collect the required boost components and versions
set(NuToBoostComponents system filesystem)
set(NuToBoostRequiredVersion 1.44.0)

if(ENABLE_MPI)
    set(NuToBoostComponents ${NuToBoostComponents} mpi)
endif()

if(ENABLE_SERIALIZATION)
    set(NuToBoostRequiredVersion 1.58.0)
    set(NuToBoostComponents ${NuToBoostComponents} serialization)
endif()

if(ENABLE_UNIT_TESTS)
    set(NuToBoostComponents ${NuToBoostComponents} unit_test_framework)
endif()

message(STATUS "Requesting BOOST components ${NuToBoostComponents} - at least version ${NuToBoostRequiredVersion} required.")

find_package(Boost ${NuToBoostRequiredVersion} COMPONENTS ${NuToBoostComponents} REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
message(STATUS "Boost version = ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_SERIALIZATION_FOUND = ${Boost_SERIALIZATION_FOUND}")
message(STATUS "Boost_SERIALIZATION_LIBRARY = ${Boost_SERIALIZATION_LIBRARY}")
message(STATUS "BOOST_INCLUDE_DIR = ${BOOST_INCLUDE_DIR}")
message(STATUS "BOOST_INCLUDE_DIRS = ${BOOST_INCLUDE_DIRS}")

# find OpenBLAS
if(ENABLE_MUMPS OR ENABLE_PARDISO)
    if(NOT OpenBLAS_ROOT AND NOT DEFINED $ENV{OpenBLAS_ROOT})
        if(NUTO_EXTERNAL_LIBRARIES_ROOT)
            set(OpenBLAS_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
        endif()
    endif()
    set(OpenBLAS_FIND_SHARED_LIBRARY TRUE)
    set(OpenBLAS_DEBUG TRUE)
    find_package(OpenBLAS REQUIRED)
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
        set(NuTo_DOXYGEN_PREDEFINED "HAVE_MUMPS" ${NuTo_DOXYGEN_PREDEFINED})
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DHAVE_MUMPS")
        include_directories(${MUMPS_INCLUDE_DIR})
        message(STATUS "MUMPS_INCLUDE_DIR = ${MUMPS_INCLUDE_DIR}")
    endif()
else()
    set(MUMPS_FOUND FALSE)
    message(STATUS "MUMPS not switched on")
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
        set(NuTo_DOXYGEN_PREDEFINED "HAVE_PARDISO" ${NuTo_DOXYGEN_PREDEFINED})
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
find_package(Eigen REQUIRED)
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


# find CGAL installation
if(ENABLE_CGAL)
    # initialize variables
    message(STATUS "Checking for CGAL headers ...")
    # check if CGAL_ROOT is set
    if(NOT CGAL_ROOT AND NOT $ENV{CGAL_ROOT} STREQUAL "")
        set(CGAL_ROOT $ENV{CGAL_ROOT})
    endif()
    if(NOT CGAL_ROOT)
        set(CGAL_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
    endif()
    find_package(CGAL)
    if(CGAL_FOUND)
        message(STATUS "CGAL_LIBRARIES = ${CGAL_LIBRARIES}")
        message(STATUS "CGAL_INCLUDE_DIR = ${CGAL_INCLUDE_DIR}")
        include_directories(${CGAL_INCLUDE_DIR})
        add_definitions("-DHAVE_CGAL")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
    else()
        message(WARNING "CGAL not found")
    endif()
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
    find_package(Arpack REQUIRED)
    if(ARPACK_FOUND)
        add_definitions("-DHAVE_ARPACK")
        include_directories(${ARPACK_INCLUDE_DIR})
    endif()
endif()

# find lapack library
if(ENABLE_MKL)
    find_package(LAPACK REQUIRED)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# build nuto
add_subdirectory(src)

# Documentation
add_subdirectory(doc)

#test Environment
include(CTest)

add_subdirectory(test/math)
if(ENABLE_MECHANICS)
    add_subdirectory(test/mechanics)
endif()

if(ENABLE_METAMODEL)
    add_subdirectory(test/metamodel)
endif()

if(ENABLE_OPTIMIZE)
    add_subdirectory(test/optimize)
endif()

if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(test/geometryConcrete)
endif()

if(ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

# myExamples should be modified locally, but are not supposed to be added to the repository
# if it is empty, there is just an empty include
# it is used for adding working examples locally, which should not be exposed to the public
if(ENABLE_MY_EXAMPLES)
    add_subdirectory(myNutoExamples)
endif()

# This makes qtcreator show header-only classes...
file(GLOB_RECURSE HeaderFiles "*.h")
add_custom_target(headers SOURCES ${HeaderFiles})
