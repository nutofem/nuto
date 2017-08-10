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
target_include_directories(Math PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/nuto>
    )

install(TARGETS Math EXPORT NuToTargets
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include/nuto
    )
