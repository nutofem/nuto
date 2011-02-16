# $Id$
SET(NuToMath_SOURCES
  FullMatrix.cpp
  Matrix.cpp
  NuToMath.cpp
  MathException.cpp
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
  ModulOperator
  ModulSparseMatrix
  ModulSparseSolver
)

# build static library
ADD_LIBRARY(NuToMath STATIC ${NuToMath_SOURCES})
TARGET_LINK_LIBRARIES(NuToMath NuToBase ${LAPACK_LIBRARIES} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MKLSolver_LIBRARIES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
 SET_TARGET_PROPERTIES(NuToMath PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

#TARGET_LINK_LIBRARIES(NuToMath NuToBase ${Boost_LIBRARIES})
#SET_TARGET_PROPERTIES(NuToMath PROPERTIES LINK_FLAGS -Wl,-z,defs)

# install shared library
INSTALL(TARGETS NuToMath ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/math)
FOREACH(NuToMath_SWIG_MODULE ${NuToMath_SWIG_MODULES})
  NUTO_SWIG_MODULE(${NuToMath_SWIG_MODULE} ${NuToMath_SWIG_MODULE}.i nuto/math NuToMath)
ENDFOREACH(NuToMath_SWIG_MODULE)
