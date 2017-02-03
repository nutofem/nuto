add_executable(LinearInterpolation LinearInterpolation.cpp)
target_link_libraries(LinearInterpolation NuToMath NuToBase)
add_test(c::unit::LinearInterpolation ${CMAKE_CURRENT_BUILD_DIR}/LinearInterpolation)

add_executable(CubicSplineInterpolation CubicSplineInterpolation.cpp)
target_link_libraries(CubicSplineInterpolation NuToMath NuToBase)
add_test(c::unit::CubicSplineInterpolation ${CMAKE_CURRENT_BUILD_DIR}/CubicSplineInterpolation)

if(ENABLE_ARPACK)
    add_executable(EigenSolverArpack EigenSolverArpack.cpp)
    target_link_libraries(EigenSolverArpack NuToMath NuToBase)
    add_test(c::unit::EigenSolverArpack ${CMAKE_CURRENT_BUILD_DIR}/EigenSolverArpack)
endif()

add_executable(EigenCompanion EigenCompanion.cpp ${CMAKE_SOURCE_DIR}/src/math/EigenCompanion.cpp)
target_link_libraries(EigenCompanion ${Boost_LIBRARIES})
add_test(c::unit::EigenCompanion ${CMAKE_CURRENT_BUILD_DIR}/EigenCompanion)

add_executable(SparseMatrixCSRGeneral SparseMatrixCSRGeneral.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRGeneral.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSR.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrix.cpp
    )
target_link_libraries(SparseMatrixCSRGeneral ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(c::unit::SparseMatrixCSRGeneral ${CMAKE_CURRENT_BUILD_DIR}/SparseMatrixCSRGeneral)

add_executable(SparseMatrixCSRVector2General SparseMatrixCSRVector2General.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRVector2General.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRVector2.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSRGeneral.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrixCSR.cpp
    ${CMAKE_SOURCE_DIR}/src/math/SparseMatrix.cpp
    )
target_link_libraries(SparseMatrixCSRVector2General ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(c::unit::SparseMatrixCSRVector2General ${CMAKE_CURRENT_BUILD_DIR}/SparseMatrixCSRVector2General)
