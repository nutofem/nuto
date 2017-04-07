foreach(test
        ConstraintNodeToElement2D
        ConstraintNodeToElement3D
        )
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    add_test(almost::mechanics::constitutive::${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()

#add_unit_test(ConstraintLinear
    #mechanics/constraints/ConstraintLinear
    #mechanics/constraints/ConstraintNode
    #mechanics/constraints/ConstraintLinearNodeDisplacements2D
    #mechanics/constraints/ConstraintBase)
