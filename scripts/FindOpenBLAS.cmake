# $Id$

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
SET(OpenBLAS_LIBRARIES)

MESSAGE(STATUS "Checking for BLAS lib ...")

# check if BLAS_ROOT is set
IF(NOT OpenBLAS_ROOT AND NOT $ENV{OpenBLAS_ROOT} STREQUAL "")
  SET(OpenBLAS_ROOT $ENV{OpenBLAS_ROOT})
ENDIF(NOT OpenBLAS_ROOT AND NOT $ENV{OpenBLAS_ROOT} STREQUAL "")

# convert path to unix style path and set search path
IF(OpenBLAS_ROOT)
  FILE(TO_CMAKE_PATH ${OpenBLAS_ROOT} OpenBLAS_ROOT)
  SET(_Openblas_LIBRARIES_SEARCH_DIRS ${OpenBLAS_ROOT}/lib ${OpenBLAS_ROOT}/openblas/lib ${OpenBLAS_ROOT} ${_Openblas_LIBRARIES_SEARCH_DIRS})
ENDIF(OpenBLAS_ROOT)

# search for BLAS libraries
IF(UNIX AND OpenBLAS_FIND_SHARED_LIBRARY)
  SET(OpenBLAS_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
ENDIF(UNIX AND OpenBLAS_FIND_SHARED_LIBRARY)

# check for BLAS libraries
FIND_LIBRARY(_Openblas_LIB_OpenBLAS NAMES openblas HINTS ${_Openblas_LIBRARIES_SEARCH_DIRS})

IF(_Openblas_LIB_OpenBLAS)
  SET(OpenBLAS_LIBRARIES ${_Openblas_LIB_OpenBLAS})
  SET(OpenBLAS_FOUND TRUE)
ENDIF(_Openblas_LIB_OpenBLAS)

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenBLAS DEFAULT_MSG OpenBLAS_LIBRARIES)


IF(OpenBLAS_DEBUG)
  MESSAGE(STATUS "OpenBLAS_FOUND=${OpenBLAS_FOUND}")
  MESSAGE(STATUS "OpenBLAS_LIBRARIES=${OpenBLAS_LIBRARIES}")
ENDIF(OpenBLAS_DEBUG)
