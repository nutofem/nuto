foreach(test
        ConstraintNodeToElement2D
        ConstraintNodeToElement3D
        Constraints
        )
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    add_test(almost::mechanics::constraints::${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()

