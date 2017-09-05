# Find the mumps solver package (sequential version)
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
set(MUMPS_LIBRARIES)
set(MUMPS_INCLUDE_DIR)

message(STATUS "Checking for MUMPS Solver package ...")

# check if MUMPS_ROOT is set
if(NOT MUMPS_ROOT AND NOT $ENV{MUMPS_ROOT} STREQUAL "")
    set(MUMPS_ROOT $ENV{MUMPS_ROOT})
endif()

if(NOT METIS_FOUND)
    find_package(METIS REQUIRED)
endif()

# convert path to unix style path and set search path
if(MUMPS_ROOT)
    file(TO_CMAKE_PATH ${MUMPS_ROOT} MUMPS_ROOT)
    set(_mumps_INCLUDE_SEARCH_DIRS ${MUMPS_ROOT}/include ${MUMPS_ROOT}
        ${_mumps_INCLUDE_SEARCH_DIRS})
    set(_mumps_LIBRARIES_SEARCH_DIRS ${MUMPS_ROOT}/lib ${MUMPS_ROOT}
        ${METIS_ROOT}/lib ${_mumps_LIBRARIES_SEARCH_DIRS})
endif()

# search for header dmumps_c.h
find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h
          HINTS ${_mumps_INCLUDE_SEARCH_DIRS})

if(UNIX AND MUMPS_FIND_STATIC_LIBRARY)
    set(MUMPS_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    find_library(_mumps_LIB_DMUMPS NAMES dmumps
        HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
    find_library(_mumps_LIB_MUMPS_COMMON
        NAMES mumps_common
        HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
    find_library(_mumps_LIB_MPISEQ
        NAMES mpiseq mpiseq_seq-4.10.0
        HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
else()
    find_library(_mumps_LIB_DMUMPS NAMES dmumps_seq-4.10.0)
    find_library(_mumps_LIB_MUMPS_COMMON NAMES mumps_common_seq-4.10.0)
    find_library(_mumps_LIB_MPISEQ NAMES mpiseq mpiseq_seq-4.10.0)
endif()

if(_mumps_LIB_DMUMPS AND _mumps_LIB_MUMPS_COMMON AND NOT _mumps_LIB_MPISEQ)
    message(WARNING "Only sequential version of MUMPS package is supported")
endif()

if(_mumps_LIB_DMUMPS AND _mumps_LIB_MUMPS_COMMON AND _mumps_LIB_MPISEQ)
    set(MUMPS_LIBRARIES ${_mumps_LIB_DMUMPS} ${_mumps_LIB_MUMPS_COMMON}
        ${_mumps_LIB_MPISEQ})
    # search for ordering packages
    find_library(_mumps_LIB_PORD
        NAMES pord
        HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
    if(_mumps_LIB_PORD)
        set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${_mumps_LIB_PORD})
    endif()

    find_library(_mumps_LIB_METIS
        NAMES metis
        HINTS ${_mumps_LIBRARIES_SEARCH_DIRS})
    if(_mumps_LIB_METIS)
        set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${_mumps_LIB_METIS})
    endif()

    set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${OpenBLAS_LIBRARIES})
endif()

if(UNIX AND MUMPS_FIND_STATIC_LIBRARY)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${MUMPS_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_LIBRARIES
    MUMPS_INCLUDE_DIR)

if(MUMPS_DEBUG)
    message(STATUS "MUMPS_FOUND=${MUMPS_FOUND}")
    message(STATUS "MUMPS_LIBRARIES=${MUMPS_LIBRARIES}")
    message(STATUS "MUMPS_INCLUDE_DIR=${MUMPS_INCLUDE_DIR}")
endif()
