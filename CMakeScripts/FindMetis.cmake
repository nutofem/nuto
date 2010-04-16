# $Id$

# find the metis package (try to find the metis header and the metis library)
#
# Stefan Eckardt, Institute of Structural Mechanics, Bauhaus-University Weimar, August 2009
#
# variables used by this module (can be also defined as environment variables):
#   METIS_ROOT - preferred installation prefix for searching for metis
#
# variables defined by this module
#   METIS_FOUND - defines whether metis was found or not
#   METIS_INCLUDE_DIR - metis include directory
#   METIS_LIBRARIES - metis libraries


# initialize variables
MESSAGE(STATUS "Checking for Metis library ...")
# check if METIS_ROOT is set
IF(NOT METIS_ROOT AND NOT $ENV{METIS_ROOT} STREQUAL "")
  SET(METIS_ROOT $ENV{METIS_ROOT})
ENDIF(NOT METIS_ROOT AND NOT $ENV{METIS_ROOT} STREQUAL "")

# convert path to unix style path and set search path
IF(METIS_ROOT)
  FILE(TO_CMAKE_PATH ${METIS_ROOT} METIS_ROOT)
  SET(_metis_INCLUDE_SEARCH_DIRS ${METIS_ROOT}/include ${METIS_ROOT} ${_metis_INCLUDE_SEARCH_DIRS})
  SET(_metis_LIBRARIES_SEARCH_DIRS ${METIS_ROOT}/lib ${METIS_ROOT} ${_metis_LIBRARIES_SEARCH_DIRS})
ENDIF(METIS_ROOT)

# search for header metis.h
FIND_PATH(METIS_INCLUDE_DIR PATH_SUFFIXES metis NAMES metis.h HINTS ${_metis_INCLUDE_SEARCH_DIRS})

# search for metis library
FIND_LIBRARY(METIS_LIBRARIES NAMES metis HINTS ${_metis_LIBRARIES_SEARCH_DIRS})

# handle the QUIETLY and REQUIRED arguments
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Metis DEFAULT_MSG METIS_LIBRARIES METIS_INCLUDE_DIR)