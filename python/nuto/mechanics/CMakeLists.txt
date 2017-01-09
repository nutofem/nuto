set(NuToMechanics_SWIG_MODULES
    ModulMechanics
    )

set(SWIG_MODULE_ModulMechanics_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/unstructured/Structure.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/elements/ElementEnum.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockStorageBase.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/DofStatus.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockScalar.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockFullMatrix.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockFullVector.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/timeIntegration/TimeIntegrationBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/timeIntegration/NewmarkBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/timeIntegration/NewmarkDirect.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureOutputBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureOutputBlockMatrix.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureOutputBlockVector.h

    ${CMAKE_SOURCE_DIR}/external/eigen.i
    )

# build python modules
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/mechanics)
if(ENABLE_VISUALIZE)
    nuto_swig_module(ModulMechanics ModulMechanics.i nuto/mechanics NuToMechanics NuToBase NuToMath NuToVisualize)
else()
    nuto_swig_module(ModulMechanics ModulMechanics.i nuto/mechanics NuToMechanics NuToBase NuToMath)
endif()
