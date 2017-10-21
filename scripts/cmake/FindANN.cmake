# Find the ANN library
# Creates IMPORTED target Ann::Ann
message(STATUS "Checking for ANN Library ...")

find_path(ANN_INCLUDE_DIR NAMES ANN/ANN.h)
find_library(ANN_LIBRARIES NAMES ANN ann)

add_library(Ann::Ann SHARED IMPORTED)
set_target_properties(Ann::Ann PROPERTIES
    IMPORTED_LOCATION "${ANN_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${ANN_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ANN DEFAULT_MSG ANN_LIBRARIES ANN_INCLUDE_DIR)
