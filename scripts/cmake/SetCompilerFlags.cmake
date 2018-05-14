macro(set_compiler_flags)
    # use newer c++ standard, enable additional warnings, support for sse4
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wall -Wextra -msse4 -pedantic -std=c++14 -fPIC")

    # before adding warnings that complain about the swig code anyway, we save
    # the flags as they are
    set(PYTHON_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

    # openmp specific
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

    # Clang specific options
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # prevent warnings from boost and eigen
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --system-header-prefix=boost/")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --system-header-prefix=eigen3/")
        # warn about missing/inconsistent documentation
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdocumentation")
        # prevent "suggest braces around initialization of subobject" warning
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -Wno-missing-braces -Wno-zero-length-array")
    endif()

    # GCC specific options
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "5")
            # warns if override is missing
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
        endif()
    endif()
endmacro()
