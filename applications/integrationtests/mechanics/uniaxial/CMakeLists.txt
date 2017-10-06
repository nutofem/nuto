foreach(test
    ElementUniaxial3D
    ElementUniaxial2D
    ElementUniaxial1D
    )

    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(${test}
        Mechanics Math Base Boost::unit_test_framework
        ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    if(MUMPS_FOUND)
        target_link_libraries(${test} Mumps::Mumps ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(integration::mechanics::uniaxial::${test}.cpp ${test})
    append_to_tests(${test})
endforeach()
