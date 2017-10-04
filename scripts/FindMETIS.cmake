# Find the metis package (try to find the metis header and the metis library)
#
# variables used by this module (can be also defined as environment variables):
#   METIS_ROOT - preferred installation prefix for searching for metis
#
# variables defined by this module
#   METIS_FOUND - defines whether metis was found or not
#   METIS_INCLUDE_DIR - metis include directory
#   METIS_LIBRARIES - metis libraries

# initialize variables
message(STATUS "Checking for Metis library ...")

# check if METIS_ROOT is set
if(NOT METIS_ROOT AND NOT $ENV{METIS_ROOT} STREQUAL "")
    set(METIS_ROOT $ENV{METIS_ROOT})
endif()

# convert path to unix style path and set search path
if(METIS_ROOT)
    file(TO_CMAKE_PATH ${METIS_ROOT} METIS_ROOT)
    set(_metis_INCLUDE_SEARCH_DIRS ${METIS_ROOT}/include ${METIS_ROOT}
        ${_metis_INCLUDE_SEARCH_DIRS})
    set(_metis_LIBRARIES_SEARCH_DIRS ${METIS_ROOT}/lib ${METIS_ROOT}
        ${_metis_LIBRARIES_SEARCH_DIRS})
endif()

# search for header metis.h
find_path(METIS_INCLUDE_DIR PATH_SUFFIXES metis
    NAMES metis.h
    HINTS ${_metis_INCLUDE_SEARCH_DIRS})

# search for metis library
find_library(METIS_LIBRARIES NAMES metis HINTS ${_metis_LIBRARIES_SEARCH_DIRS})

message(STATUS "Metis include: ${METIS_INCLUDE_DIR}")
message(STATUS "Metis lib: ${METIS_LIBRARIES}")

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG METIS_LIBRARIES
    METIS_INCLUDE_DIR)
