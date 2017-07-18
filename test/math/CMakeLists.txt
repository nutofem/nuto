add_unit_test(Average)

add_unit_test(LinearInterpolation
    math/Interpolation.cpp
    )

add_unit_test(CubicSplineInterpolation
    math/Interpolation.cpp
    )

add_unit_test(EigenCompanion)

add_unit_test(SparseMatrixCSRGeneral
    math/SparseMatrixCSR.cpp
    )

add_unit_test(SparseMatrixCSRVector2General
        math/SparseMatrixCSRGeneral.cpp
        math/SparseMatrixCSR.cpp
    )


add_unit_test(Legendre)
add_unit_test(NaturalCoordinateMemoizer)
add_unit_test(NewtonRaphson
    math/SparseMatrixCSR.cpp
    math/SparseDirectSolverMUMPS.cpp
    base/Timer.cpp
    )
target_link_libraries(NewtonRaphson ${MUMPS_LIBRARIES})

add_unit_test(Gmres)


add_unit_test(SpatialContainer)
target_link_libraries(SpatialContainer ${ANN_LIBRARIES})

if(ENABLE_ARPACK)
    add_unit_test(EigenSolverArpack)
    target_link_libraries(EigenSolverArpack Math Base)
endif()
