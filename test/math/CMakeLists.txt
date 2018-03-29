add_unit_test(Average)
add_unit_test(LinearInterpolation math/Interpolation.cpp)
add_unit_test(CubicSplineInterpolation math/Interpolation.cpp)
add_unit_test(EigenCompanion)
add_unit_test(EigenIO)

add_unit_test(EigenSparseSolve)
if(SUITESPARSE_FOUND)
    target_link_libraries(EigenSparseSolve SuiteSparse::UmfPack SuiteSparse::Cholmod)
endif()
if(MUMPS_FOUND)
    target_link_libraries(EigenSparseSolve EigenMumpsSupport)
endif()
if(MKL_FOUND)
    target_link_libraries(EigenSparseSolve ${mkl_link_libs})
endif()


add_unit_test(Legendre)
add_unit_test(NaturalCoordinateMemoizer)
add_unit_test(NewtonRaphson
    math/EigenSparseSolve.cpp
    )

add_unit_test(Gmres)
add_unit_test(Shapes
    mechanics/elements/ElementShapeFunctions.cpp
    )
