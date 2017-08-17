foreach(test
    ElementUniaxial3D
    ElementUniaxial2D
    ElementUniaxial1D
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(integration::mechanics::uniaxial::${test}.cpp ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
