# Creates IMPORTED target SuiteSparse::SuiteSparse
message(STATUS "Checking for SuiteSparse Solver package ...")

find_path(SUITESPARSE_INCLUDE_DIR NAMES umfpack.h HINTS /usr/include/suitesparse)

find_library(SUITESPARSE_LIBRARIES NAMES umfpack)

add_library(SuiteSparse::SuiteSparse SHARED IMPORTED)
set_target_properties(SuiteSparse::SuiteSparse PROPERTIES
    IMPORTED_LOCATION "${SUITESPARSE_LIBRARIES}"
    INTERFACE_COMPILE_DEFINITIONS "HAVE_SUITESPARSE"
    INTERFACE_INCLUDE_DIRECTORIES "${SUITESPARSE_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUITESPARSE DEFAULT_MSG SUITESPARSE_LIBRARIES SUITESPARSE_INCLUDE_DIR)
