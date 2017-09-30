message(STATUS "Checking for Boost...")

set(Boost_DETAILED_FAILURE_MSG TRUE)
if(NOT BOOST_ROOT AND NOT DEFINED $ENV{BOOST_ROOT})
    if(NUTO_EXTERNAL_LIBRARIES_ROOT)
        set(BOOST_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
    endif()
endif()

# collect the required boost components and versions
set(NuToBoostComponents system filesystem)
set(NuToBoostRequiredVersion 1.49.0)

if(ENABLE_MPI)
    set(NuToBoostComponents ${NuToBoostComponents} mpi)
endif()

set(NuToBoostComponents ${NuToBoostComponents} unit_test_framework)

message(STATUS "Requesting BOOST components ${NuToBoostComponents}")
message(STATUS "At least version ${NuToBoostRequiredVersion} required.")

find_package(Boost ${NuToBoostRequiredVersion}
    COMPONENTS ${NuToBoostComponents} REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
set(boostVersion
    "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
message(STATUS "Boost version = ${boostVersion}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
message(STATUS "BOOST_INCLUDE_DIR = ${BOOST_INCLUDE_DIR}")
message(STATUS "BOOST_INCLUDE_DIRS = ${BOOST_INCLUDE_DIRS}")
