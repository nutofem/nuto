set(MathSources
    CubicSplineInterpolation.cpp
    EigenIO.cpp
    Interpolation.cpp
    LinearInterpolation.cpp
    Legendre.cpp
)

create_nuto_module(Math "${MathSources}")
target_link_libraries(Math Base ${LAPACK_LIBRARIES})
