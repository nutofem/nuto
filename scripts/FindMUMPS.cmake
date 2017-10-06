# Find the mumps solver package (sequential version)
# Creates IMPORTED target Mumps::Mumps
message(STATUS "Checking for MUMPS Solver package ...")

find_package(METIS REQUIRED)
find_package(BLAS REQUIRED)

find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h
    HINTS /usr/include/mumps-seq-shared/)

find_library(_mumps_LIB_DMUMPS NAMES dmumps_seq dmumps)
find_library(_mumps_LIB_MUMPS_COMMON NAMES mumps_common_seq mumps_common)
find_library(_mumps_LIB_PORD NAMES pord)

set(MUMPS_LIBRARIES ${_mumps_LIB_DMUMPS} ${_mumps_LIB_MUMPS_COMMON}
    ${_mumps_LIB_PORD})

set(link_libs
    ${_mumps_LIB_MUMPS_COMMON}
    ${_mumps_LIB_PORD}
    ${METIS_LIBRARIES}
    ${BLAS_LIBRARIES})

add_library(Mumps::Mumps SHARED IMPORTED)
set_target_properties(Mumps::Mumps PROPERTIES
    IMPORTED_LOCATION "${_mumps_LIB_DMUMPS}"
    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${link_libs}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_LIBRARIES
    MUMPS_INCLUDE_DIR)
