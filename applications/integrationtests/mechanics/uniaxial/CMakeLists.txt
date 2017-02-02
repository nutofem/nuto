if(ENABLE_PYTHON)
    add_test(p::mechanics::ElementUniaxialTest ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/ElementUniaxialTest.py
        ${CMAKE_SYSTEM_NAME}
    ${CMAKE_SYSTEM_PROCESSOR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})
    set_property(TEST p::mechanics::ElementUniaxialTest PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endif()

foreach(test
    ElementUniaxialBrick3D
    ElementUniaxialPrism3D
    ElementUniaxialQuad2D
    ElementUniaxialTetrahedron3D
    ElementUniaxialTriangle2D
    ElementUniaxialTruss1D
    )

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    if(MUMPS_FOUND)
        target_link_libraries(${test} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif()
    add_test(c::mechanics::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
