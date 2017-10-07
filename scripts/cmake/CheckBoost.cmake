set(Components filesystem unit_test_framework)
set(RequiredVersion 1.49.0)

if(ENABLE_MPI)
    set(Components ${Components} mpi)
endif()

message(STATUS "Requesting BOOST components ${Components}")
message(STATUS "At least version ${RequiredVersion} required.")

find_package(Boost ${RequiredVersion} COMPONENTS ${Components} REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
set(boostVersion
    "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
message(STATUS "Boost version = ${boostVersion}")
