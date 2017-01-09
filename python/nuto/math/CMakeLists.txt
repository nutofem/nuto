set(NuToMath_SWIG_MODULES
    NuToMath
    ModulMatrix
    ModulOperator
    ModulSparseMatrix
    ModulSparseSolver
)

set(SWIG_MODULE_ModulSparseMatrix_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrix.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSR.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRGeneral_Def.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRSymmetric_Def.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRVector2General_Def.h
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRVector2Symmetric_Def.h
)

configure_file(${NUMPY_I_PATH}/numpy.i ${CMAKE_CURRENT_BINARY_DIR}/numpy.i COPYONLY)

set(SEARCH_REGEX "NO_IMPORT_ARRAY")

file(READ ${CMAKE_CURRENT_BINARY_DIR}/numpy.i FILE_CONTENT)

string(REGEX REPLACE "${SEARCH_REGEX}" "IMPORT_ARRAY"
    MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/numpy.i "${MODIFIED_FILE_CONTENT}")

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/src)

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(ModulFullMatrix.i PROPERTIES CPLUSPLUS ON)

# build python modules
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/math)
foreach(NuToMath_SWIG_MODULE ${NuToMath_SWIG_MODULES})
    nuto_swig_module(${NuToMath_SWIG_MODULE} ${NuToMath_SWIG_MODULE}.i nuto/math NuToMath)
endforeach()
