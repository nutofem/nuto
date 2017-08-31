# find the Arpack library
#
# variables used by this module (can be also defined as environment variables):
#   ARPACK_ROOT - preferred installation prefix for searching for ARPACK
#   ARPACK_FIND_STATIC_LIBRARY - searches for static libraries (UNIX only)
#   ARPACK_DEBUG - print debug messages
#
# variables defined by this module
#   ARPACK_FOUND - defines whether metis was found or not
#   ARPACK_INCLUDE_DIR - ARPACK include directory
#   ARPACK_LIBRARIES   - ARPACK libraries

# initialize variables
message(STATUS "Checking for ARPACK Library ...")

# check if ARPACK_ROOT is set
if(NOT ARPACK_ROOT AND NOT $ENV{ARPACK_ROOT} STREQUAL "")
    set(ARPACK_ROOT $ENV{ARPACK_ROOT})
endif()

# convert path to unix style path and set search path
if(ARPACK_ROOT)
    file(TO_CMAKE_PATH ${ARPACK_ROOT} ARPACK_ROOT)
    set(_ARPACK_LIBRARIES_SEARCH_DIRS ${ARPACK_ROOT}/lib ${ARPACK_ROOT}
        ${_ARPACK_LIBRARIES_SEARCH_DIRS})
endif()

# search for ARPACK library
if(UNIX AND ARPACK_FIND_STATIC_LIBRARY)
    set(ARPACK_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

find_library(ARPACK_LIBRARIES NAMES arpack
    HINTS ${_ARPACK_LIBRARIES_SEARCH_DIRS})

if(UNIX AND ARPACK_FIND_STATIC_LIBRARY)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${ARPACK_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG ARPACK_LIBRARIES)

if(ARPACK_DEBUG)
    message(STATUS "ARPACK_FOUND=${ARPACK_FOUND}")
    message(STATUS "ARPACK_LIBRARIES=${ARPACK_LIBRARIES}")
endif()

