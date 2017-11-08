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

    # TRILINOS_PATH to the path to your Trilinos install.
# You do _not_ need to edit this line.
FIND_PACKAGE(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})

# If FIND_PACKAGE successfully found your Trilinos install, it will
# set the Boolean flag Trilinos_FOUND.  The following IF statement
# fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
# found, it prints out the values of some Trilinos configuration
# details.  You may find them useful for building your application
# that uses Trilinos.
IF(Trilinos_FOUND)
   MESSAGE("\nFound Trilinos!  Here are the details: ")
   MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
   MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
   MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
   MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
   MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
   MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
   MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
   MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
   MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
   MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
   MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS =
     ${Trilinos_Fortran_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
   MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
   MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
   MESSAGE("End of Trilinos details\n")

   include_directories(${Trilinos_INCLUDE_DIRS})
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos!")
ENDIF()


endmacro()
