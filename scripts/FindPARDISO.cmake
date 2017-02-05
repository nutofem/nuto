# $Id$

# find the PARDISO solver package
#
# variables used by this module (can be also defined as environment variables):
#   PARDISO_ROOT - preferred installation prefix for searching for mumps
#   PARDISO_DEBUG - print debug messages
#
# variables defined by this module
#   PARDISO_FOUND - defines whether mumps was found or not
#   PARDISO_LIBRARIES - mumps libraries

# set variables
SET(PARDISO_LIBRARIES)

MESSAGE(STATUS "Checking for PARDISO Solver package ...")

# check if PARDISO_ROOT is set
IF(NOT PARDISO_ROOT AND NOT $ENV{PARDISO_ROOT} STREQUAL "")
  SET(PARDISO_ROOT $ENV{PARDISO_ROOT})
ENDIF(NOT PARDISO_ROOT AND NOT $ENV{PARDISO_ROOT} STREQUAL "")

# convert path to unix style path and set search path
IF(PARDISO_ROOT)
  FILE(TO_CMAKE_PATH ${PARDISO_ROOT} PARDISO_ROOT)
  SET(_pardiso_LIBRARIES_SEARCH_DIRS ${PARDISO_ROOT}/lib ${PARDISO_ROOT} ${_pardiso_LIBRARIES_SEARCH_DIRS})
ENDIF(PARDISO_ROOT)

# search for pardiso libraries
IF(UNIX AND PARDISO_FIND_SHARED_LIBRARY)
  SET(PARDISO_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
ENDIF(UNIX AND PARDISO_FIND_SHARED_LIBRARY)
# check for pardiso libraries
FIND_LIBRARY(_pardiso_LIB_PARDISO NAMES pardiso HINTS ${_pardiso_LIBRARIES_SEARCH_DIRS})

IF(_pardiso_LIB_PARDISO)
  SET(PARDISO_LIBRARIES ${_pardiso_LIB_PARDISO} ${OpenBLAS_LIBRARIES})
ENDIF(_pardiso_LIB_PARDISO)

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARDISO DEFAULT_MSG PARDISO_LIBRARIES)

IF(PARDISO_DEBUG)
  MESSAGE(STATUS "PARDISO_FOUND=${PARDISO_FOUND}")
  MESSAGE(STATUS "PARDISO_LIBRARIES=${PARDISO_LIBRARIES}")
ENDIF(PARDISO_DEBUG)
