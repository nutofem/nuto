set(NuToMathSources
    CubicSplineInterpolation.cpp
    EigenCompanion.cpp
    EigenSolverArpack.cpp
    Interpolation.cpp
    LinearInterpolation.cpp
    NuToMath.cpp
    SparseMatrix.cpp
    SparseMatrixCSR.cpp
    SparseMatrixCSRGeneral.cpp
    SparseMatrixCSRSymmetric.cpp
    SparseMatrixCSRVector2.cpp
    SparseMatrixCSRVector2General.cpp
    SparseMatrixCSRVector2Symmetric.cpp
    SparseDirectSolverMKLDSS.cpp
    SparseDirectSolverMKLPardiso.cpp
    SparseDirectSolverPardiso.cpp
    SparseDirectSolverMUMPS.cpp
)

sources_to_objects("${NuToMathSources}" NuToMathObjects)

# build static library
add_library(NuToMath ${NuToMathObjects})
target_link_libraries(NuToMath NuToBase ${LAPACK_LIBRARIES} ${MUMPS_LIBRARIES} ${PARDISO_LIBRARIES} ${ARPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MKLSolver_LIBRARIES})

# install library
install(TARGETS NuToMath
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})
