set(MathSources
    CubicSplineInterpolation.cpp
    EigenIO.cpp
    EigenSparseSolve.cpp
    Interpolation.cpp
    Legendre.cpp
    LinearInterpolation.cpp
    PolynomialLeastSquaresFitting.cpp
    Quadrature.cpp
)

create_nuto_module(Math "${MathSources}")
target_link_libraries(Math PRIVATE Base ${LAPACK_LIBRARIES})
target_link_libraries(Math PUBLIC Eigen3::Eigen)

if(SUITESPARSE_FOUND)
    target_link_libraries(Math PRIVATE SuiteSparse::UmfPack SuiteSparse::Cholmod)
endif()
if(MUMPS_FOUND)
    target_link_libraries(Math PRIVATE EigenMumpsSupport)
endif()
