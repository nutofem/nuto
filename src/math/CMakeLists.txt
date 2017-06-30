set(NuToMathSources
    CubicSplineInterpolation.cpp
    EigenCompanion.cpp
    EigenSolverArpack.cpp
    Interpolation.cpp
    LinearInterpolation.cpp
    Legendre.cpp
    SparseMatrixCSR.cpp
    SparseMatrixCSRGeneral.cpp
    SparseMatrixCSRSymmetric.cpp
    SparseMatrixCSRVector2General.cpp
    SparseMatrixCSRVector2Symmetric.cpp
    SparseDirectSolverMKLPardiso.cpp
    SparseDirectSolverPardiso.cpp
    SparseDirectSolverMUMPS.cpp
)

sources_to_objects("${NuToMathSources}" NuToMathObjects)

add_library(NuToMath ${NuToMathObjects})
target_link_libraries(NuToMath NuToBase ${LAPACK_LIBRARIES} ${MUMPS_LIBRARIES} ${PARDISO_LIBRARIES} ${ARPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MKLSolver_LIBRARIES})

install(TARGETS NuToMath
    LIBRARY DESTINATION lib
    )
