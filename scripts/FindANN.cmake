# Find the ANN library
#
# variables used by this module (can be also defined as environment variables):
#   ANN_ROOT - preferred installation prefix for searching for ANN
#   ANN_FIND_STATIC_LIBRARY - searches for static libraries (UNIX only)
#   ANN_DEBUG - print debug messages
#
# variables defined by this module
#   ANN_FOUND - defines whether metis was found or not
#   ANN_INCLUDE_DIR - ANN include directory
#   ANN_LIBRARIES   - ANN libraries

# initialize variables
message(STATUS "Checking for ANN Library ...")

# check if ANN_ROOT is set
if(NOT ANN_ROOT AND NOT $ENV{ANN_ROOT} STREQUAL "")
    set(ANN_ROOT $ENV{ANN_ROOT})
endif()

# convert path to unix style path and set search path
if(ANN_ROOT)
    file(TO_CMAKE_PATH ${ANN_ROOT} ANN_ROOT)
    set(_ANN_INCLUDE_SEARCH_DIRS ${ANN_ROOT}/include ${ANN_ROOT}
        ${_ANN_INCLUDE_SEARCH_DIRS})
    set(_ANN_LIBRARIES_SEARCH_DIRS ${ANN_ROOT}/lib ${ANN_ROOT}
        ${_ANN_LIBRARIES_SEARCH_DIRS})
endif()

# search for header ANN.h
find_path(ANN_INCLUDE_DIR NAMES ANN/ANN.h HINTS ${_ANN_INCLUDE_SEARCH_DIRS})

# search for ann library
if(UNIX AND ANN_FIND_STATIC_LIBRARY)
    set(ANN_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

find_library(ANN_LIBRARIES NAMES ANN ann HINTS ${_ANN_LIBRARIES_SEARCH_DIRS})

if(UNIX AND ANN_FIND_STATIC_LIBRARY)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${ANN_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ANN DEFAULT_MSG ANN_LIBRARIES ANN_INCLUDE_DIR)

if(ANN_DEBUG)
    message(STATUS "ANN_FOUND=${ANN_FOUND}")
    message(STATUS "ANN_INCLUDE_DIR=${ANN_INCLUDE_DIR}")
    message(STATUS "ANN_LIBRARIES=${ANN_LIBRARIES}")
endif()
