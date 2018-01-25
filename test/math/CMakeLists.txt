add_unit_test(Average)
add_unit_test(LinearInterpolation math/Interpolation.cpp)
add_unit_test(CubicSplineInterpolation math/Interpolation.cpp)
add_unit_test(EigenCompanion)
add_unit_test(EigenIO)
add_unit_test(Legendre)
add_unit_test(NaturalCoordinateMemoizer)
add_unit_test(NewtonRaphson
    base/Timer.cpp
    )

add_unit_test(Gmres)
add_unit_test(SpatialContainer)
target_link_libraries(SpatialContainer Ann::Ann)
