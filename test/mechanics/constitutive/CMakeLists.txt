if(ENABLE_UNIT_TESTS)
    foreach(test HeatConduction ConstitutiveIO AdditiveInput AdditiveOutput GradientDamageFatigue)
        add_executable(${test} ${test}.cpp)
        target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
        if(ENABLE_VISUALIZE)
            target_link_libraries(${test} NuToVisualize)
        endif()
        if(DEFINED ENV{JENKINS_HOME})
            add_test(c::unit::${test} sh -c "${CMAKE_CURRENT_BINARY_DIR}/${test} -f JUNIT > TEST-${test}.xml")
        else()
            add_test(c::unit::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
        endif()
    endforeach()


    # unit tests without linking the full nutomechanics lib
    add_executable(DataContainer DataContainer.cpp)
    target_link_libraries(DataContainer ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase NuToMath)
    if(DEFINED ENV{JENKINS_HOME})
        add_test(c::unit::DataContainer sh -c "${CMAKE_CURRENT_BINARY_DIR}/DataContainer -f JUNIT > TEST-DataContainer.xml")
    else()
        add_test(c::unit::DataContainer ${CMAKE_CURRENT_BINARY_DIR}/DataContainer)
    endif()

    add_executable(IPConstitutiveLaw IPConstitutiveLaw.cpp
            ${CMAKE_SOURCE_DIR}/src/nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp)
    target_link_libraries(IPConstitutiveLaw ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase NuToMath)
    if(DEFINED ENV{JENKINS_HOME})
        add_test(c::unit::IPConstitutiveLaw sh -c "${CMAKE_CURRENT_BINARY_DIR}/IPConstitutiveLaw -f JUNIT > TEST-IPConstitutiveLaw.xml")
    else()
        add_test(c::unit::IPConstitutiveLaw ${CMAKE_CURRENT_BUILD_DIR}/IPConstitutiveLaw)
    endif()

endif()
