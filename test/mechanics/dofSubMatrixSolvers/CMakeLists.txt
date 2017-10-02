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
