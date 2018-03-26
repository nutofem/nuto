set(Components filesystem unit_test_framework)
set(RequiredVersion 1.49.0)

message(STATUS "Looking for Boost components ${Components} of at least version ${RequiredVersion}.")
find_package(Boost ${RequiredVersion} COMPONENTS ${Components} REQUIRED)
