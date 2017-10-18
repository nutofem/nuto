set(mechanicsDir ${CMAKE_SOURCE_DIR}/src/mechanics)
set(SWIG_MODULE_ModuleMechanics_EXTRA_DEPS
    ${mechanicsDir}/constraints/Constraints.h
    ${mechanicsDir}/constraints/ConstraintCompanion.h
    ${mechanicsDir}/constraints/Equation.h

    ${mechanicsDir}/dofSubMatrixStorage/DofStatus.h
    ${mechanicsDir}/dofSubMatrixStorage/BlockScalar.h
    ${mechanicsDir}/dofSubMatrixStorage/BlockFullMatrix.h
    ${mechanicsDir}/dofSubMatrixStorage/BlockFullVector.h
    ${mechanicsDir}/dofSubMatrixStorage/BlockSparseMatrix.h
    ${mechanicsDir}/dofSubMatrixStorage/BlockStorageBase.h

    ${mechanicsDir}/elements/ElementEnum.h

    ${mechanicsDir}/PDEs/HygroHelpers.h
    ${mechanicsDir}/PDEs/PorousMedium.h

    ${mechanicsDir}/nodes/NodeBase.h

    ${mechanicsDir}/timeIntegration/TimeIntegrationBase.h
    ${mechanicsDir}/timeIntegration/NewmarkDirect.h
    ${mechanicsDir}/timeIntegration/TimeControl.h
    ${mechanicsDir}/timeIntegration/postProcessing/PostProcessor.h

    ${mechanicsDir}/structures/Assembler.h
    ${mechanicsDir}/structures/StructureBase.h
    ${mechanicsDir}/structures/StructureOutputBase.h
    ${mechanicsDir}/structures/StructureOutputBlockMatrix.h
    ${mechanicsDir}/structures/StructureOutputBlockVector.h
    ${mechanicsDir}/structures/unstructured/Structure.h

    ${CMAKE_SOURCE_DIR}/external/eigen.i
    ${CMAKE_CURRENT_SOURCE_DIR}/Sections.i
    )

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
if(ENABLE_VISUALIZE)
    nuto_swig_module(ModuleMechanics ModuleMechanics.i
        Mechanics Base Math Visualize)
else()
    nuto_swig_module(ModuleMechanics ModuleMechanics.i Mechanics Base Math)
endif()
