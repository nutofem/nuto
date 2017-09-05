set(solverSources
    base/Timer.cpp
    mechanics/dofSubMatrixStorage/BlockSparseMatrix.cpp
    mechanics/dofSubMatrixStorage/DofStatus.cpp
    mechanics/dofSubMatrixStorage/BlockFullVector.cpp
    mechanics/dofSubMatrixStorage/BlockStorageBase.cpp
    mechanics/dofSubMatrixStorage/BlockScalar.cpp
    mechanics/nodes/NodeEnum.cpp
    math/SparseMatrixCSRSymmetric.cpp
    math/SparseMatrixCSRGeneral.cpp
    math/SparseMatrixCSR.cpp
    )

add_unit_test(SolverEigen
    ${solverSources}
    )

add_unit_test(SolverMUMPS
    math/SparseDirectSolverMUMPS.cpp
    ${solverSources}
    )
target_link_libraries(SolverMUMPS ${MUMPS_LIBRARIES})

if(PARDISO_FOUND)
    add_unit_test(SolverPardiso
        math/SparseDirectSolverPardiso.cpp
        ${solverSources}
        )
    target_link_libraries(SolverPardiso ${PARDISO_LIBRARIES}
        ${OpenBLAS_LIBRARIES})
endif()

#add_unit_test(SparseSolver
#        mechanics/dofSubMatrixStorage/BlockSparseMatrix.cpp
#        mechanics/dofSubMatrixStorage/DofStatus.cpp
#        mechanics/dofSubMatrixStorage/BlockFullVector.cpp
#        mechanics/dofSubMatrixStorage/BlockStorageBase.cpp
#        mechanics/dofSubMatrixStorage/BlockScalar.cpp
#        mechanics/nodes/NodeEnum.cpp
#        math/SparseMatrixCSRVector2.cpp
#        math/SparseMatrixCSRVector2General.cpp
#        math/SparseMatrixCSRSymmetric.cpp
#        math/SparseMatrixCSRVector2Symmetric.cpp
#        math/SparseMatrixCSRGeneral.cpp
#        math/SparseMatrixCSR.cpp
#        math/SparseMatrix.cpp
#        math/SparseDirectSolverMUMPS.cpp
#        math/SparseDirectSolverPardiso.cpp
#        )
#if (PARDISO_FOUND)
#
#add_unit_test(SpatialContainer)
#target_link_libraries(SpatialContainer ${ANN_LIBRARIES})
#
#if(ENABLE_ARPACK)
#    add_unit_test(EigenSolverArpack)
#    target_link_libraries(EigenSolverArpack NuToMath)
#endif()
