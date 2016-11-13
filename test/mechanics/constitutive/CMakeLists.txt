if(ENABLE_UNIT_TESTS)
    foreach(test HeatConduction ConstitutiveIO AdditiveInput AdditiveOutput GradientDamageFatigue LocalDamageModel)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        add_test(c::unit::${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
    endforeach()


    # unit tests without linking the full nutomechanics lib
    add_executable(DataContainer DataContainer.cpp)
    target_link_libraries(DataContainer ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase)
    add_test(c::unit::DataContainer ${CMAKE_CURRENT_BUILD_DIR}/DataContainer)

    add_executable(IPConstitutiveLaw IPConstitutiveLaw.cpp
            ${CMAKE_SOURCE_DIR}/src/nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp)

    target_link_libraries(IPConstitutiveLaw ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase)
    add_test(c::unit::IPConstitutiveLaw ${CMAKE_CURRENT_BUILD_DIR}/IPConstitutiveLaw)

endif()
