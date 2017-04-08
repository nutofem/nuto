set(NuToMath_SWIG_MODULES
    NuToMath
    ModuleTestInterface
    ModuleSparseMatrix
    ModuleSparseSolver
)

set(SWIG_MODULE_ModuleSparseMatrix_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrix.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSR.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRGeneral_Def.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRSymmetric_Def.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRVector2General_Def.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRVector2Symmetric_Def.h
)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/src)

# build python modules
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
foreach(NuToMath_SWIG_MODULE ${NuToMath_SWIG_MODULES})
    nuto_swig_module(${NuToMath_SWIG_MODULE} ${NuToMath_SWIG_MODULE}.i NuToMath)
endforeach()
