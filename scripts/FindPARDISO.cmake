# Find the PARDISO solver package
#
# variables used by this module (can be also defined as environment variables):
#   PARDISO_ROOT - preferred installation prefix for searching for mumps
#
# variables defined by this module
#   PARDISO_FOUND - defines whether mumps was found or not
#   PARDISO_LIBRARIES - mumps libraries

# set variables
set(PARDISO_LIBRARIES)

message(STATUS "Checking for PARDISO Solver package ...")

# check if PARDISO_ROOT is set
if(NOT PARDISO_ROOT AND NOT $ENV{PARDISO_ROOT} STREQUAL "")
    set(PARDISO_ROOT $ENV{PARDISO_ROOT})
endif()

# convert path to unix style path and set search path
if(PARDISO_ROOT)
    file(TO_CMAKE_PATH ${PARDISO_ROOT} PARDISO_ROOT)
    set(_pardiso_LIBRARIES_SEARCH_DIRS ${PARDISO_ROOT}/lib ${PARDISO_ROOT}
        ${_pardiso_LIBRARIES_SEARCH_DIRS})
endif()

# check for pardiso libraries
find_library(_pardiso_LIB_PARDISO
    NAMES pardiso
    HINTS ${_pardiso_LIBRARIES_SEARCH_DIRS})

if(_pardiso_LIB_PARDISO)
    set(PARDISO_LIBRARIES ${_pardiso_LIB_PARDISO} ${BLAS_LIBRARIES})
endif()

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARDISO DEFAULT_MSG PARDISO_LIBRARIES)

message(STATUS "PARDISO_FOUND=${PARDISO_FOUND}")
message(STATUS "PARDISO_LIBRARIES=${PARDISO_LIBRARIES}")
