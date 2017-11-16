function(setup_google_benchmark)
    # Store nuto cxx flags and set library specific flags
    set(ORIGINAL_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-std=c++14 -fPIC -pthread")

    add_subdirectory(external/benchmark EXCLUDE_FROM_ALL)

    # Restore nuto cxx flags
    set(CMAKE_CXX_FLAGS "${ORIGINAL_CXX_FLAGS}")
endfunction()
