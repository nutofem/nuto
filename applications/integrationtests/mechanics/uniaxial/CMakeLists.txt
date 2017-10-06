foreach(test
    ElementUniaxial3D
    ElementUniaxial2D
    ElementUniaxial1D
    )

    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(${test}
        Mechanics Math Base Visualize
        Boost::unit_test_framework ${LAPACK_LIBRARIES})
    target_include_directories(${test} PRIVATE ${CMAKE_SOURCE_DIR}/test/tools)
    add_test(integration::mechanics::uniaxial::${test}.cpp ${test})
    append_to_tests(${test})
endforeach()
