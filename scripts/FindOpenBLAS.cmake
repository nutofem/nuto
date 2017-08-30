# find the OpenBLAS lib
#
# variables used by this module (can be also defined as environment variables):
#   OpenBLAS_ROOT - preferred installation prefix for searching for Openblas
#   OpenBLAS_DEBUG - print debug messages
#
# variables defined by this module
#   OpenBLAS_FOUND - defines whether Openblas was found or not
#   OpenBLAS_LIBRARIES - Openblas libraries

# set variables
set(OpenBLAS_LIBRARIES)

message(STATUS "Checking for BLAS lib ...")

# check if BLAS_ROOT is set
if(NOT OpenBLAS_ROOT AND NOT $ENV{OpenBLAS_ROOT} STREQUAL "")
    set(OpenBLAS_ROOT $ENV{OpenBLAS_ROOT})
endif()

# convert path to unix style path and set search path
if(OpenBLAS_ROOT)
    file(TO_CMAKE_PATH ${OpenBLAS_ROOT} OpenBLAS_ROOT)
    set(_Openblas_LIBRARIES_SEARCH_DIRS ${OpenBLAS_ROOT}/lib
        ${OpenBLAS_ROOT}/openblas/lib ${OpenBLAS_ROOT}
        ${_Openblas_LIBRARIES_SEARCH_DIRS})
endif()

# search for BLAS libraries
if(UNIX AND OpenBLAS_FIND_SHARED_LIBRARY)
    set(OpenBLAS_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES
        ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
endif()

# check for BLAS libraries
find_library(_Openblas_LIB_OpenBLAS
    NAMES openblas
    HINTS ${_Openblas_LIBRARIES_SEARCH_DIRS})

if(_Openblas_LIB_OpenBLAS)
    set(OpenBLAS_LIBRARIES ${_Openblas_LIB_OpenBLAS})
    set(OpenBLAS_FOUND TRUE)
endif()

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenBLAS DEFAULT_MSG OpenBLAS_LIBRARIES)


if(OpenBLAS_DEBUG)
    message(STATUS "OpenBLAS_FOUND=${OpenBLAS_FOUND}")
    message(STATUS "OpenBLAS_LIBRARIES=${OpenBLAS_LIBRARIES}")
endif()
