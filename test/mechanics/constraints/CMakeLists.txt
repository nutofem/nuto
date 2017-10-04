foreach(test
        ConstraintNodeToElement2D
        ConstraintNodeToElement3D
        )
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base
        Boost::unit_test_framework)
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    add_test(almost::mechanics::constraints::${test}
        ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()

add_unit_test(ConstraintCompanion
    mechanics/groups/GroupBase.cpp
    mechanics/groups/Group.cpp)

add_unit_test(Constraints
    mechanics/nodes/NodeEnum.cpp
    mechanics/nodes/NodeBase.cpp
    math/SparseMatrixCSRGeneral.cpp
    math/SparseMatrixCSR.cpp)

