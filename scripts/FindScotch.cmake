
# find the scotch package for metis (try to find the scotch header and the scotch libraries)
#
# Stefan Eckardt, Institute of Structural Mechanics, Bauhaus-University Weimar, August 2009
#
# variables used by this module (can be also defined as environment variables):
#   SCOTCH_ROOT - preferred installation prefix for searching for scotch
#
# variables defined by this module
#   SCOTCH_FOUND - defines whether scotch was found or not
#   SCOTCH_INCLUDE_DIR - scotch include directory
#   SCOTCH_LIBRARIES - scotch libraries

# set variables
SET(Scotch_LIBRARIES)
SET(Scotch_INCLUDE_DIR)

MESSAGE(STATUS "Checking for Scotch library ...")

# check if SCOTCH_ROOT is set
IF(NOT SCOTCH_ROOT AND NOT $ENV{SCOTCH_ROOT} STREQUAL "")
  SET(SCOTCH_ROOT $ENV{SCOTCH_ROOT})
ENDIF(NOT SCOTCH_ROOT AND NOT $ENV{SCOTCH_ROOT} STREQUAL "")

# convert path to unix style path and set search path
IF(SCOTCH_ROOT)
  FILE(TO_CMAKE_PATH ${SCOTCH_ROOT} SCOTCH_ROOT)
  SET(_scotch_INCLUDE_SEARCH_DIRS ${SCOTCH_ROOT}/include ${SCOTCH_ROOT} ${_scotch_INCLUDE_SEARCH_DIRS})
  SET(_scotch_LIBRARIES_SEARCH_DIRS ${SCOTCH_ROOT}/lib ${SCOTCH_ROOT} ${_scotch_LIBRARIES_SEARCH_DIRS})
ENDIF(SCOTCH_ROOT)

# search for header scotch.h
FIND_PATH(SCOTCH_INCLUDE_DIR NAMES scotch.h HINTS ${_scotch_INCLUDE_SEARCH_DIRS})

# search for scotch libraries
FIND_LIBRARY(_scotch_LIB_ESMUMPS NAMES esmumps HINTS ${_scotch_LIBRARIES_SEARCH_DIRS})
FIND_LIBRARY(_scotch_LIB_SCOTCH NAMES scotch HINTS ${_scotch_LIBRARIES_SEARCH_DIRS})
FIND_LIBRARY(_scotch_LIB_SCOTCHERR NAMES scotcherr HINTS ${_scotch_LIBRARIES_SEARCH_DIRS})
IF(_scotch_LIB_ESMUMPS AND _scotch_LIB_SCOTCH AND _scotch_LIB_SCOTCHERR)
  SET(SCOTCH_LIBRARIES ${_scotch_LIB_ESMUMPS} ${_scotch_LIB_SCOTCH} ${_scotch_LIB_SCOTCHERR})
ENDIF(_scotch_LIB_ESMUMPS AND _scotch_LIB_SCOTCH AND _scotch_LIB_SCOTCHERR)

# handle the QUIETLY and REQUIRED arguments
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Scotch DEFAULT_MSG SCOTCH_LIBRARIES SCOTCH_INCLUDE_DIR)
