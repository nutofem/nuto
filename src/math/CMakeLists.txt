set(MathSources
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

sources_to_objects("${MathSources}" MathObjects)

add_library(Math ${MathObjects})
set_target_properties(Math PROPERTIES
    OUTPUT_NAME NuToMath
    )

target_link_libraries(Math Base ${LAPACK_LIBRARIES} ${MUMPS_LIBRARIES} ${PARDISO_LIBRARIES} ${ARPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MKLSolver_LIBRARIES})

install(TARGETS Math EXPORT NuToTargets
    LIBRARY DESTINATION lib
    )
