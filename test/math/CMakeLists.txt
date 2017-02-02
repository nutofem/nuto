add_executable(LinearInterpolation LinearInterpolation.cpp)
target_link_libraries(LinearInterpolation NuToMath NuToBase)
add_test(c::unit::LinearInterpolation ${CMAKE_CURRENT_BUILD_DIR}/LinearInterpolation)

add_executable(CubicSplineInterpolation CubicSplineInterpolation.cpp)
target_link_libraries(CubicSplineInterpolation NuToMath NuToBase)
add_test(c::unit::CubicSplineInterpolation ${CMAKE_CURRENT_BUILD_DIR}/CubicSplineInterpolation)

if(ENABLE_ARPACK)
    add_executable(Arpack Arpack.cpp)
    target_link_libraries(Arpack NuToMath NuToBase)
    add_test(c::unit::Arpack ${CMAKE_CURRENT_BUILD_DIR}/Arpack)
endif()

add_executable(EigenCompanion EigenCompanion.cpp ${CMAKE_SOURCE_DIR}/src/math/EigenCompanion.cpp)
target_link_libraries(EigenCompanion ${Boost_LIBRARIES})
add_test(c::unit::EigenCompanion ${CMAKE_CURRENT_BUILD_DIR}/EigenCompanion)
