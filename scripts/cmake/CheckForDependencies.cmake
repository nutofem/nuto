macro(check_for_dependencies)
    if(ENABLE_MPI)
        message(STATUS "Find MPI")

        find_package(MPI REQUIRED)
        include_directories(SYSTEM ${MPI_INCLUDE_PATH})

        # These variables are set by FindMPI
        message(STATUS "MPI_CXX_FOUND:              ${MPI_CXX_FOUND}")
        message(STATUS "MPI_CXX_COMPILER:           ${MPI_CXX_COMPILER}")
        message(STATUS "MPI_CXX_COMPILER_FLAGS:     ${MPI_CXX_COMPILER_FLAGS}")
        message(STATUS "MPI_CXX_INCLUDE_PATH:       ${MPI_CXX_INCLUDE_PATH}")
        message(STATUS "MPI_CXX_LINK_FLAGS:         ${MPI_CXX_LINK_FLAGS}")
        message(STATUS "MPI_CXX_LIBRARIES:          ${MPI_CXX_LIBRARIES}")

        # Additionally, FindMPI sets the following variables for running
        # MPI programs from the command line:
        message(STATUS "MPIEXEC:                    ${MPIEXEC}")
        message(STATUS "MPIEXEC_NUMPROC_FLAG:       ${MPIEXEC_NUMPROC_FLAG}")
        message(STATUS "MPIEXEC_PREFLAG:            ${MPIEXEC_PREFLAG}")
        message(STATUS "MPIEXEC_POSTFLAG:           ${MPIEXEC_POSTFLAG}")
    endif()

    if(ENABLE_MKL)
        message(STATUS "Checking for MKL LAPACK and MKL SOLVER ...")
        find_package(MKLSolver)
        if(MKLSolver_FOUND)
            message(STATUS "MKLSolver_LIBRARIES = ${MKLSolver_LIBRARIES}")
            message(STATUS "MKLSolver_INCLUDE_DIR = ${MKLSolver_INCLUDE_DIR}")
            message(STATUS "MKLSolver_DEFINITIONS = ${MKLSolver_DEFINITIONS}")
            add_definitions(${MKLSolver_DEFINITIONS})
            set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};${MKLSolver_DEFINITIONS}")
            include_directories(${MKLSolver_INCLUDE_DIR})
        endif()
        find_package(LAPACK REQUIRED)
    endif()

    include(CheckBoost)

    find_package(MUMPS QUIET)
    if(MUMPS_FOUND)
        message(STATUS "MUMPS solver found.")
    else()
        warning("MUMPS solver not found - will not be available.")
    endif()

    find_package(PARDISO QUIET)
    if(PARDISO_FOUND)
        message(STATUS "Pardsio solver found.")
    else()
        warning("Pardiso solver not found - will not be available.")
    endif()

    # find Eigen header files (Linear Algebra)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    message(STATUS "Eigen version = ${EIGEN3_VERSION_STRING}")

    # find ARPACK library for eigenvalue analysis
    find_package(ARPACK QUIET)
    if(ARPACK_FOUND)
        message(STATUS "Arpack solver found.")
    else()
        warning("Arpack not found. EigenSolverArpack will not be available.")
    endif()

    find_package(benchmark QUIET)
    if(benchmark_FOUND)
        message(STATUS "Found google benchmark")
    else()
        warning("Google benchmark NOT found - Benchmarks won't be build!")
    endif()

    find_package(SUITESPARSE QUIET)
    if(SUITESPARSE_FOUND)
        message(STATUS "Found SuiteSparse solver library")
    else()
        warning("SuiteSparse NOT found - certain solvers won't be available")
    endif()
endmacro()
