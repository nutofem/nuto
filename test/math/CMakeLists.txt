add_unit_test(LinearInterpolation
    math/Interpolation.cpp 
    )

add_unit_test(CubicSplineInterpolation
    math/Interpolation.cpp 
    )

add_unit_test(EigenCompanion)

add_unit_test(SparseMatrixCSRGeneral
    math/SparseMatrixCSR.cpp
    math/SparseMatrix.cpp
    )

add_unit_test(SparseMatrixCSRVector2General
        math/SparseMatrixCSRVector2.cpp
        math/SparseMatrixCSRGeneral.cpp
        math/SparseMatrixCSR.cpp
        math/SparseMatrix.cpp
    )


add_unit_test(SparseSolver
        mechanics/dofSubMatrixStorage/BlockSparseMatrix.cpp
        mechanics/dofSubMatrixStorage/DofStatus.cpp
        mechanics/dofSubMatrixStorage/BlockFullVector.cpp
        mechanics/dofSubMatrixStorage/BlockStorageBase.cpp
        mechanics/dofSubMatrixStorage/BlockScalar.cpp
        mechanics/nodes/NodeEnum.cpp
        math/SparseMatrixCSRVector2.cpp
        math/SparseMatrixCSRVector2General.cpp
        math/SparseMatrixCSRSymmetric.cpp
        math/SparseMatrixCSRVector2Symmetric.cpp
        math/SparseMatrixCSRGeneral.cpp
        math/SparseMatrixCSR.cpp
        math/SparseMatrix.cpp
        math/SparseDirectSolverMUMPS.cpp
        )
target_link_libraries(SparseSolver ${MUMPS_LIBRARIES})


add_unit_test(SpatialContainer)
target_link_libraries(SpatialContainer ${ANN_LIBRARIES})

if(ENABLE_ARPACK)
    add_unit_test(EigenSolverArpack)
    target_link_libraries(EigenSolverArpack NuToMath)
endif()
