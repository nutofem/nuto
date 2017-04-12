set(SWIG_MODULE_ModuleMechanics_EXTRA_DEPS
    ${CMAKE_SOURCE_DIR}/src/mechanics/constraints/Constraints.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/constraints/ConstraintCompanion.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/DofStatus.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockScalar.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockFullMatrix.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockFullVector.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/dofSubMatrixStorage/BlockStorageBase.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/elements/ElementEnum.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/nodes/NodeBase.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/timeIntegration/TimeIntegrationBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/timeIntegration/NewmarkBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/timeIntegration/NewmarkDirect.h

    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/Assembler.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureOutputBase.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureOutputBlockMatrix.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/StructureOutputBlockVector.h
    ${CMAKE_SOURCE_DIR}/src/mechanics/structures/unstructured/Structure.h

    ${CMAKE_SOURCE_DIR}/external/eigen.i
    ${CMAKE_CURRENT_SOURCE_DIR}/Sections.i
    )

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
if(ENABLE_VISUALIZE)
    nuto_swig_module(ModuleMechanics ModuleMechanics.i NuToMechanics NuToBase NuToMath NuToVisualize)
else()
    nuto_swig_module(ModuleMechanics ModuleMechanics.i NuToMechanics NuToBase NuToMath)
endif()
