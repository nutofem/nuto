# - Find MKL Solver 
# This module finds an installed mkl solver libraries
#
# This module sets the following variables:
#  MKLSolver_FOUND - set to true if the mkl solver libraries are found
#  MKLSolver_LIBRARIES - uncached list of mkl libraries (using full path name)
#  MKLSolver_INCLUDE_DIR - mkl include directory: uncached
#  MKLSolver_DEFINITIONS - definitions which solver is found
#  MKLSolver_DOXYGEN_DEFINITIONS - definitions for doxygen which solver is found
#  MKLSolver_STATIC  if set on this determines what kind of linkage we do (static)
##########

# set variables
SET(MKLSolver_FOUND FALSE)
SET(MKLSolver_LIBRARIES)
SET(MKLSolver_INCLUDE_DIR)
SET(MKLSolver_DEFINITIONS)
SET(MKLSolver_DOXYGEN_DEFINITIONS)

# check for mkl_solvers
FIND_PATH(MKLSolver_INCLUDE_DIR mkl_solver.h)
IF(MKLSolver_INCLUDE_DIR)
  FIND_PACKAGE(Threads REQUIRED)

  # find dss solver
  FIND_PATH(MKLSolver_DSS_INCLUDE_DIR mkl_dss.h)
  IF(MKLSolver_DSS_INCLUDE_DIR)
    #MESSAGE(STATUS "MKL DSS-Solver found.")
    SET(MKLSolver_DEFINITIONS "-DHAVE_MKL_DSS" ${MKLSolver_DEFINITIONS})
    SET(MKLSolver_DOXYGEN_DEFINITIONS "HAVE_MKL_DSS" ${MKLSolver_DOXYGEN_DEFINITIONS})
  ENDIF(MKLSolver_DSS_INCLUDE_DIR)

  # find pardiso solver
  FIND_PATH(MKLSolver_PARDISO_INCLUDE_DIR mkl_pardiso.h)
  IF(MKLSolver_PARDISO_INCLUDE_DIR)
    #MESSAGE(STATUS "MKL PARDISO-Solver found.")
    SET(MKLSolver_DEFINITIONS "-DHAVE_MKL_PARDISO" ${MKLSolver_DEFINITIONS})
    SET(MKLSolver_DOXYGEN_DEFINITIONS "HAVE_MKL_PARDISO" ${MKLSolver_DOXYGEN_DEFINITIONS})
  ENDIF(MKLSolver_PARDISO_INCLUDE_DIR)
  IF(MKLSolver_DSS_INCLUDE_DIR OR MKLSolver_PARDISO_INCLUDE_DIR)
    # find mkl libraries
    IF(UNIX)
      SET(MKLSolver_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
      # static or shared libraries
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
      FIND_LIBRARY(MKLSolver_LIB_MKL_INTEL_LIBRARY NAMES mkl_intel_lp64 PATHS ENV LD_LIBRARY_PATH)
      MESSAGE(STATUS "MKLSolver_LIB_MKL_INTEL_LIBRARY = ${MKLSolver_LIB_MKL_INTEL_LIBRARY}.")
      FIND_LIBRARY(MKLSolver_LIB_MKL_LIBRARY NAMES mkl_core PATHS ENV LD_LIBRARY_PATH)
      MESSAGE(STATUS "MKLSolver_LIB_MKL_LIBRARY = ${MKLSolver_LIB_MKL_LIBRARY}.")
      FIND_LIBRARY(MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY NAMES mkl_sequential PATHS ENV LD_LIBRARY_PATH)
      MESSAGE(STATUS "MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY = ${MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY}.")
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ${MKLSolver_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
      IF(MKLSolver_LIB_MKL_INTEL_LIBRARY AND MKLSolver_LIB_MKL_LIBRARY AND MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY)
        SET(MKLSolver_FOUND TRUE)
        SET(BLAS_FOUND TRUE)
        SET(BLAS_LIBRARIES "-Wl,--start-group ${MKLSolver_LIB_MKL_INTEL_LIBRARY} ${MKLSolver_LIB_MKL_LIBRARY} ${MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY}  -Wl,--end-group" ${CMAKE_THREAD_LIBS_INIT})
        SET(MKLSolver_LIBRARIES "-Wl,--start-group ${MKLSolver_LIB_MKL_INTEL_LIBRARY} ${MKLSolver_LIB_MKL_LIBRARY} ${MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY}  -Wl,--end-group" ${CMAKE_THREAD_LIBS_INIT})
      ENDIF(MKLSolver_LIB_MKL_INTEL_LIBRARY AND MKLSolver_LIB_MKL_LIBRARY AND MKLSolver_LIB_MKL_SEQUENTIAL_LIBRARY)
    ENDIF(UNIX)
  ENDIF(MKLSolver_DSS_INCLUDE_DIR OR MKLSolver_PARDISO_INCLUDE_DIR)
ENDIF(MKLSolver_INCLUDE_DIR)
