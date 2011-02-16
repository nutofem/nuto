# $Id$

# Eigen header files (try to find the header file Eigen/Core)
#
# Jörg F. Unger, Institute of Structural Mechanics, Bauhaus-University Weimar, October 2009
#
# variables used by this module (can be also defined as environment variables):
#   EIGEN_ROOT - preferred installation prefix for searching for eigen package
#
# variables defined by this module
#   EIGEN_FOUND - defines whether metis was found or not
#   EIGEN_INCLUDE_DIR - mersenne include directory


# initialize variables
MESSAGE(STATUS "Checking for Eigen2 headers ...")
# check if EIGEN_ROOT is set
IF(NOT EIGEN_ROOT AND NOT $ENV{EIGEN_ROOT} STREQUAL "")
  SET(EIGEN_ROOT $ENV{EIGEN_ROOT})
ENDIF(NOT EIGEN_ROOT AND NOT $ENV{EIGEN_ROOT} STREQUAL "")

# convert path to unix style path and set search path
IF(EIGEN_ROOT)
  FILE(TO_CMAKE_PATH ${EIGEN_ROOT} EIGEN_ROOT)
  SET(_eigen_INCLUDE_SEARCH_DIRS ${EIGEN_ROOT}/include ${EIGEN_ROOT} ${_eigen_INCLUDE_SEARCH_DIRS})
ENDIF(EIGEN_ROOT)

# search for header dSFMT.h
FIND_PATH(EIGEN_INCLUDE_DIR NAMES eigen2/Eigen/Core HINTS ${_eigen_INCLUDE_SEARCH_DIRS})

# handle the QUIETLY and REQUIRED arguments
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR)
