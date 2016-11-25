set(NuToMath_SOURCES
    CubicSplineInterpolation.cpp
    EigenSolverArpack.cpp
    FullMatrix.cpp
    Interpolation.cpp
    LinearInterpolation.cpp
    Matrix.cpp
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
    SparseDirectSolverPardiso.cpp
    SparseDirectSolverMUMPS.cpp
)

set(NuToMath_SWIG_MODULES
    NuToMath
    ModulMatrix
    ModulFullMatrix
    ModulFullVector
    ModulOperator
    ModulSparseMatrix
    ModulSparseSolver
)

set(SWIG_MODULE_ModulFullMatrix_EXTRA_DEPS
    FullMatrix_Def.h
)

set(SWIG_MODULE_ModulFullVectorx_EXTRA_DEPS
    FullVector_Def.h
)

set(SWIG_MODULE_ModulSparseMatrix_EXTRA_DEPS
    SparseMatrix.h
    SparseMatrixCSR.h
    SparseMatrixCSRGeneral_Def.h
    SparseMatrixCSRSymmetric_Def.h
    SparseMatrixCSRVector2General_Def.h
    SparseMatrixCSRVector2Symmetric_Def.h
)

# build static library
add_library(NuToMath ${NuToMath_SOURCES})
target_link_libraries(NuToMath NuToBase ${LAPACK_LIBRARIES} ${MUMPS_LIBRARIES} ${PARDISO_LIBRARIES} ${ARPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MKLSolver_LIBRARIES})
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set_target_properties(NuToMath PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

if(ENABLE_NUMPY)
    configure_file(${NUMPY_I_PATH}/numpy.i ${CMAKE_CURRENT_BINARY_DIR}/numpy.i COPYONLY)
endif()

set(SEARCH_REGEX "NO_IMPORT_ARRAY")

if(ENABLE_NUMPY)
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/numpy.i FILE_CONTENT)
endif()

string(REGEX REPLACE "${SEARCH_REGEX}" "IMPORT_ARRAY"
    MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/numpy.i "${MODIFIED_FILE_CONTENT}")

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(ModulFullMatrix.i PROPERTIES CPLUSPLUS ON)

# install library
install(TARGETS NuToMath
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/math)
foreach(NuToMath_SWIG_MODULE ${NuToMath_SWIG_MODULES})
    nuto_swig_module(${NuToMath_SWIG_MODULE} ${NuToMath_SWIG_MODULE}.i nuto/math NuToMath)
endforeach()
