# $Id$

# find the mumps solver package, sequential version (try to find the mumps header and the mumps libraries)
#
# Stefan Eckardt, Institute of Structural Mechanics, Bauhaus-University Weimar, August 2009
#
# variables used by this module (can be also defined as environment variables):
#   MUMPS_ROOT - preferred installation prefix for searching for mumps
#   MUMPS_FIND_STATIC_LIBRARY - searches for static libraries (UNIX only)
#   MUMPS_DEBUG - print debug messages
#
# variables defined by this module
#   MUMPS_FOUND - defines whether mumps was found or not
#   MUMPS_INCLUDE_DIR - mumps include directory
#   MUMPS_LIBRARIES - mumps libraries

# set variables
SET(MUMPS_LIBRARIES)
SET(MUMPS_INCLUDE_DIR)

MESSAGE(STATUS "Checking for MUMPS Solver package ...")

# check if MUMPS_ROOT is set
IF(NOT MUMPS_ROOT AND NOT $ENV{MUMPS_ROOT} STREQUAL "")
  SET(MUMPS_ROOT $ENV{MUMPS_ROOT})
ENDIF(NOT MUMPS_ROOT AND NOT $ENV{MUMPS_ROOT} STREQUAL "")

IF(NOT METIS_FOUND)
  FIND_PACKAGE(Metis REQUIRED)
ENDIF(NOT METIS_FOUND)

# convert path to unix style path and set search path
IF(MUMPS_ROOT)
  FILE(TO_CMAKE_PATH ${MUMPS_ROOT} MUMPS_ROOT)
  SET(_mumps_INCLUDE_SEARCH_DIRS ${MUMPS_ROOT}/include ${MUMPS_ROOT} ${_mumps_INCLUDE_SEARCH_DIRS})
  SET(_mumps_LIBRARIES_SEARCH_DIRS ${MUMPS_ROOT}/lib ${MUMPS_ROOT} ${METIS_ROOT}/lib ${_mumps_LIBRARIES_SEARCH_DIRS})
ENDIF(MUMPS_ROOT)

# search for header dmumps_c.h
FIND_PATH(MUMPS_INCLUDE_DIR NAMES dmumps_c.h HINTS ${_mumps_INCLUDE_SEARCH_DIRS})

# search for mumps libraries
IF(UNIX AND MUMPS_FIND_STATIC_LIBRARY)
  SET(MUMPS_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
ENDIF(UNIX AND MUMPS_FIND_STATIC_LIBRARY)
# check for mumps libraries
FIND_LIBRARY(_mumps_LIB_DMUMPS NAMES dmumps HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
FIND_LIBRARY(_mumps_LIB_MUMPS_COMMON NAMES mumps_common HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
FIND_LIBRARY(_mumps_LIB_MPISEQ NAMES mpiseq HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
IF(_mumps_LIB_DMUMPS AND _mumps_LIB_MUMPS_COMMON AND NOT _mumps_LIB_MPISEQ)
  MESSAGE(WARNING "Only sequential version of MUMPS package is supported")
ENDIF(_mumps_LIB_DMUMPS AND _mumps_LIB_MUMPS_COMMON AND NOT _mumps_LIB_MPISEQ)
IF(_mumps_LIB_DMUMPS AND _mumps_LIB_MUMPS_COMMON AND _mumps_LIB_MPISEQ)
  SET(MUMPS_LIBRARIES ${_mumps_LIB_DMUMPS} ${_mumps_LIB_MUMPS_COMMON} ${_mumps_LIB_MPISEQ})
  # search for ordering packages
  FIND_LIBRARY(_mumps_LIB_PORD NAMES pord HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
  IF(_mumps_LIB_PORD)
    SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${_mumps_LIB_PORD})
  ENDIF(_mumps_LIB_PORD)
  FIND_LIBRARY(_mumps_LIB_METIS NAMES metis HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
  IF(_mumps_LIB_METIS)
    SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${_mumps_LIB_METIS})
  ENDIF(_mumps_LIB_METIS)
  SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${OpenBLAS_LIBRARIES})
  
  # add fortran library for gnu and clang compilers
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} "-lgfortran")
  endif()
ENDIF(_mumps_LIB_DMUMPS AND _mumps_LIB_MUMPS_COMMON AND _mumps_LIB_MPISEQ)
IF(UNIX AND MUMPS_FIND_STATIC_LIBRARY)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ${MUMPS_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
ENDIF(UNIX AND MUMPS_FIND_STATIC_LIBRARY)

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_LIBRARIES MUMPS_INCLUDE_DIR)

IF(MUMPS_DEBUG)
  MESSAGE(STATUS "MUMPS_FOUND=${MUMPS_FOUND}")
  MESSAGE(STATUS "MUMPS_LIBRARIES=${MUMPS_LIBRARIES}")
  MESSAGE(STATUS "MUMPS_INCLUDE_DIR=${MUMPS_INCLUDE_DIR}")
ENDIF(MUMPS_DEBUG)
  
  
