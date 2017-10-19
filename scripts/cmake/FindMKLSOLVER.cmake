# Find MKL Solver
# This module finds an installed mkl solver libraries
#
# This module sets the following variables:
#  MKLSolver_FOUND - set to true if the mkl solver libraries are found
#  MKLSolver_LIBRARIES - uncached list of mkl libraries (using full path name)
#  MKLSolver_INCLUDE_DIR - mkl include directory: uncached
#  MKLSolver_DEFINITIONS - definitions which solver is found

# set variables
set(MKLSolver_FOUND FALSE)
set(MKLSolver_LIBRARIES)
set(MKLSolver_INCLUDE_DIR)
set(MKLSolver_DEFINITIONS)
set(MKLSolver_DOXYGEN_DEFINITIONS)

# check for mkl_solvers
find_path(MKLSolver_INCLUDE_DIR mkl_solver.h)
if(MKLSolver_INCLUDE_DIR)
    find_package(Threads REQUIRED)

    # find dss solver
    find_path(MKLSolver_DSS_INCLUDE_DIR mkl_dss.h)
    if(MKLSolver_DSS_INCLUDE_DIR)
        set(MKLSolver_DEFINITIONS "-DHAVE_MKL_DSS" ${MKLSolver_DEFINITIONS})
    endif()

    # find pardiso solver
    find_path(MKLSolver_PARDISO_INCLUDE_DIR mkl_pardiso.h)
    if(MKLSolver_PARDISO_INCLUDE_DIR)
        set(MKLSolver_DEFINITIONS "-DHAVE_MKL_PARDISO" ${MKLSolver_DEFINITIONS})
    endif()

    # find mkl libraries
    if(MKLSolver_DSS_INCLUDE_DIR OR MKLSolver_PARDISO_INCLUDE_DIR)
        set(MKLSolver_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES
            ${CMAKE_FIND_LIBRARY_SUFFIXES})
        # static or shared libraries
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        find_library(MKLSolver_LIB_MKL_INTEL_LIBRARY
            NAMES mkl_intel_lp64
            PATHS ENV LD_LIBRARY_PATH)
        find_library(MKLSolver_LIB_MKL_LIBRARY
            NAMES mkl_core
            PATHS ENV LD_LIBRARY_PATH)
        find_library(MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY
            NAMES mkl_sequential
            PATHS ENV LD_LIBRARY_PATH)
        set(CMAKE_FIND_LIBRARY_SUFFIXES
            ${MKLSolver_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

        if(MKLSolver_LIB_MKL_INTEL_LIBRARY
                AND MKLSolver_LIB_MKL_LIBRARY
                AND MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY)
            set(MKLSolver_FOUND TRUE)
            set(BLAS_FOUND TRUE)
            set(BLAS_LIBRARIES "-Wl,--start-group "
                "${MKLSolver_LIB_MKL_INTEL_LIBRARY} "
                "${MKLSolver_LIB_MKL_LIBRARY} "
                "${MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY}  -Wl,--end-group"
                ${CMAKE_THREAD_LIBS_INIT})
            set(MKLSolver_LIBRARIES "-Wl,--start-group "
                "${MKLSolver_LIB_MKL_INTEL_LIBRARY} "
                "${MKLSolver_LIB_MKL_LIBRARY} "
                "${MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY}  -Wl,--end-group"
                ${CMAKE_THREAD_LIBS_INIT})
        endif()
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKLSOLVER DEFAULT_MSG MKLSolver_LIBRARIES
    MKLSolver_INCLUDE_DIR)
