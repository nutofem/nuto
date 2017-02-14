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

add_unit_test(SpacialContainer)
target_link_libraries(SpacialContainer ${ANN_LIBRARIES})

if(ENABLE_ARPACK)
    add_unit_test(EigenSolverArpack)
    target_link_libraries(EigenSolverArpack NuToMath)
endif()
