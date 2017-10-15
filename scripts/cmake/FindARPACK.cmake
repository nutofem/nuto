# Find the Arpack library. Creates Arpack::Arpack target.
# ARPACK_FOUND tells you if Arpack was found.
message(STATUS "Checking for ARPACK Library ...")

find_library(ARPACK_LIBRARIES NAMES arpack)

add_library(Arpack::Arpack SHARED IMPORTED)
set_target_properties(Arpack::Arpack PROPERTIES
    IMPORTED_LOCATION "${ARPACK_LIBRARIES}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG ARPACK_LIBRARIES)

