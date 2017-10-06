# Find the metis package.
message(STATUS "Checking for Metis library ...")

find_library(METIS_LIBRARIES NAMES metis)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG METIS_LIBRARIES)
