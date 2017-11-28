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
    endif()

    include(CheckBoost)

    # find mumps solver
    if(ENABLE_MUMPS)
        find_package(MUMPS REQUIRED)
        add_definitions("-DHAVE_MUMPS")
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DHAVE_MUMPS")
    endif()

    # find pardiso solver
    if(ENABLE_PARDISO)
        find_package(PARDISO REQUIRED)
        add_definitions("-DHAVE_PARDISO")
        set(NuTo_SWIG_FLAGS "${NuTo_SWIG_FLAGS};-DHAVE_PARDISO")
    endif()

    find_package(ANN REQUIRED)

    # find Eigen header files (Linear Algebra)
    find_package(EIGEN 3.2 REQUIRED)
    message(STATUS "EIGEN_VERSION_NUMBER = ${EIGEN_VERSION_NUMBER}")
    include_directories(${EIGEN_INCLUDE_DIR})

    # find ARPACK library for eigenvalue analysis
    find_package(ARPACK)
    if(ARPACK_FOUND)
        add_definitions("-DHAVE_ARPACK")
    else()
        warning("Arpack not found. EigenSolverArpack will not be available.")
    endif()

    # find lapack library
    if(ENABLE_MKL)
        find_package(LAPACK REQUIRED)
    endif()

    find_package(benchmark QUIET)
    if(benchmark_FOUND)
        message(STATUS "Found google benchmark")
    else()
        warning("Google benchmark NOT found - Benchmarks won't be build!")
    endif()
endmacro()
