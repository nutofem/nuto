# $Id$

# find the CGAL package (try to find the CGAL header and the CGAL library)
#
# Daniel Arnold, Institute of Structural Mechanics, Bauhaus-Universitaet Weimar, April 2012
#
# variables used by this module (can be also defined as environment variables):
#   CGAL_ROOT - preferred installation prefix for searching for CGAL package
#
# variables defined by this module
#   CGAL_FOUND - defines whether CGAL was found or not
#   CGAL_INCLUDE_DIR - CGAL include directory
#   CGAL_LIBRARIES - CGAL libraries


# convert path to unix style path and set search path
IF(CGAL_ROOT)
  FILE(TO_CMAKE_PATH ${CGAL_ROOT} CGAL_ROOT)
  SET(_CGAL_INCLUDE_SEARCH_DIRS ${CGAL_ROOT}/include ${CGAL_ROOT} ${_CGAL_INCLUDE_SEARCH_DIRS})
  SET(_CGAL_LIBRARIES_SEARCH_DIRS ${CGAL_ROOT}/lib ${CGAL_ROOT} ${_CGAL_LIBRARIES_SEARCH_DIRS})
ENDIF(CGAL_ROOT)

# search for header metis.h
FIND_PATH(CGAL_INCLUDE_DIR NAMES "CGAL/basic.h" HINTS ${_CGAL_INCLUDE_SEARCH_DIRS})
IF(CGAL_INCLUDE_DIR)
  set( CGAL_FOUND TRUE )
ELSE(CGAL_INCLUDE_DIR)
  set( CGAL_FOUND FALSE )
ENDIF(CGAL_INCLUDE_DIR)

# search for CGAL library
FIND_LIBRARY(CGAL_LIBRARIES NAMES "CGAL" HINTS ${_CGAL_LIBRARIES_SEARCH_DIRS})
IF(CGAL_LIBRARIES AND CGAL_FOUND)
    set( CGAL_FOUND TRUE )
ELSE(CGAL_LIBRARIES AND CGAL_FOUND)
  set( CGAL_FOUND FALSE )
ENDIF(CGAL_LIBRARIES AND CGAL_FOUND)

# handle the QUIETLY and REQUIRED arguments
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CGAL DEFAULT_MSG CGAL_INCLUDE_DIR)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CGAL DEFAULT_MSG CGAL_LIBRARIES)
