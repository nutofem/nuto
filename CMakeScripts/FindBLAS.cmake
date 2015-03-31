# $Id$

# find the BLAS lib
#
# variables used by this module (can be also defined as environment variables):
#   BLAS_ROOT - preferred installation prefix for searching for blas
#   BLAS_DEBUG - print debug messages
#
# variables defined by this module
#   BLAS_FOUND - defines whether blas was found or not
#   BLAS_LIBRARIES - blas libraries

# set variables
SET(BLAS_LIBRARIES)

MESSAGE(STATUS "Checking for BLAS lib ...")

# check if BLAS_ROOT is set
IF(NOT BLAS_ROOT AND NOT $ENV{BLAS_ROOT} STREQUAL "")
  SET(BLAS_ROOT $ENV{BLAS_ROOT})
ENDIF(NOT BLAS_ROOT AND NOT $ENV{BLAS_ROOT} STREQUAL "")

# convert path to unix style path and set search path
IF(BLAS_ROOT)
  FILE(TO_CMAKE_PATH ${BLAS_ROOT} BLAS_ROOT)
  SET(_blas_LIBRARIES_SEARCH_DIRS ${BLAS_ROOT}/lib ${BLAS_ROOT} ${_blas_LIBRARIES_SEARCH_DIRS})
ENDIF(BLAS_ROOT)

# search for BLAS libraries
IF(UNIX AND BLAS_FIND_SHARED_LIBRARY)
  SET(BLAS_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
ENDIF(UNIX AND BLAS_FIND_SHARED_LIBRARY)

# check for BLAS libraries
FIND_LIBRARY(_blas_LIB_BLAS NAMES openblas HINTS ${_blas_LIBRARIES_SEARCH_DIRS})

IF(_blas_LIB_BLAS)
  SET(BLAS_LIBRARIES ${_blas_LIB_BLAS})
ENDIF(_blas_LIB_BLAS)

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLAS DEFAULT_MSG BLAS_LIBRARIES)


IF(BLAS_DEBUG)
  MESSAGE(STATUS "BLAS_FOUND=${BLAS_FOUND}")
  MESSAGE(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
ENDIF(BLAS_DEBUG)
