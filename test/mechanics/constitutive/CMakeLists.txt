if(ENABLE_UNIT_TESTS)
    foreach(test HeatConduction ConstitutiveIO AdditiveInput AdditiveOutput StaticData)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        add_test(c::unit::${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
    endforeach()

    add_executable(Leaf Leaf.cpp
        ${CMAKE_SOURCE_DIR}/src/nuto/mechanics/constitutive/staticData/DataMoistureTransport.cpp
        ${CMAKE_SOURCE_DIR}/src/nuto/math/Matrix.cpp ${CMAKE_SOURCE_DIR}/src/nuto/math/FullMatrix.cpp)
    target_link_libraries(Leaf ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(c::unit::Leaf ${CMAKE_CURRENT_BUILD_DIR}/Leaf)

    add_executable(Composite Composite.cpp
        ${CMAKE_SOURCE_DIR}/src/nuto/mechanics/constitutive/staticData/Composite.cpp)
    target_link_libraries(Composite ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(c::unit::Composite ${CMAKE_CURRENT_BUILD_DIR}/Composite)

endif()
