# $Id$
SET(NuToMath_SOURCES
  FullMatrix.cpp
  Matrix.cpp
  MathException.cpp
  NuToMath.cpp
  SparseMatrix.cpp
  SparseMatrixCSR.cpp
  SparseMatrixCSRGeneral.cpp
  SparseMatrixCSRSymmetric.cpp
  SparseMatrixCSRVector2.cpp
  SparseMatrixCSRVector2General.cpp
  SparseMatrixCSRVector2Symmetric.cpp
  SparseDirectSolverMKLDSS.cpp
  SparseDirectSolverMKLPardiso.cpp
  SparseDirectSolverMUMPS.cpp
)

SET(NuToMath_SWIG_MODULES
  NuToMath
  ModulMatrix
  ModulFullMatrix
  ModulFullVector
  ModulOperator
  ModulSparseMatrix
  ModulSparseSolver
)

# build static library
ADD_LIBRARY(NuToMath ${NuToMath_SOURCES})
TARGET_LINK_LIBRARIES(NuToMath NuToBase ${LAPACK_LIBRARIES} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MKLSolver_LIBRARIES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
 SET_TARGET_PROPERTIES(NuToMath PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

IF(ENABLE_NUMPY)
    configure_file(${NUMPY_I_PATH}/numpy.i ${CMAKE_CURRENT_BINARY_DIR}/numpy.i COPYONLY)
ENDIF(ENABLE_NUMPY)

SET (SEARCH_REGEX "NO_IMPORT_ARRAY")
# replace text in all files
message("Replacing text")

 FILE(READ ${CMAKE_CURRENT_BINARY_DIR}/numpy.i FILE_CONTENT)
 STRING(REGEX REPLACE "${SEARCH_REGEX}" "IMPORT_ARRAY" 
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
 FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/numpy.i "${MODIFIED_FILE_CONTENT}")

                       ################

include_directories(/home/Desktop/project/eigen)
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(ModulFullMatrix.i PROPERTIES CPLUSPLUS ON)


#### exist in NuToMacros.cmake ln98
#SWIG_ADD_MODULE(ModulFullMatrix python ModulFullMatrix.i)
#SWIG_LINK_LIBRARIES(ModulFullMatrix ${PYTHON_LIBRARIES})





#####################################

#TARGET_LINK_LIBRARIES(NuToMath NuToBase ${Boost_LIBRARIES})
#SET_TARGET_PROPERTIES(NuToMath PROPERTIES LINK_FLAGS -Wl,-z,defs)

# install library
INSTALL(TARGETS NuToMath
	ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
	LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/math)
FOREACH(NuToMath_SWIG_MODULE ${NuToMath_SWIG_MODULES})
  NUTO_SWIG_MODULE(${NuToMath_SWIG_MODULE} ${NuToMath_SWIG_MODULE}.i nuto/math NuToMath)
ENDFOREACH(NuToMath_SWIG_MODULE)
