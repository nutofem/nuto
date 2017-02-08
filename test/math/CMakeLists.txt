add_unit_test(LinearInterpolation
    ${CMAKE_SOURCE_DIR}/src/math/Interpolation.cpp 
    )

add_unit_test(CubicSplineInterpolation
    ${CMAKE_SOURCE_DIR}/src/math/Interpolation.cpp 
    )

add_unit_test(EigenCompanion)

add_unit_test(SparseMatrixCSRGeneral
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSR.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrix.cpp
    )

add_unit_test(SparseMatrixCSRVector2General
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRVector2.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRGeneral.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSR.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrix.cpp
    )

add_unit_test(SpacialContainer)

if(ENABLE_ARPACK)
    add_executable(EigenSolverArpack EigenSolverArpack.cpp)
    target_link_libraries(EigenSolverArpack NuToMath NuToBase)
    add_test(unit::math::EigenSolverArpack ${CMAKE_CURRENT_BUILD_DIR}/EigenSolverArpack)
endif()
