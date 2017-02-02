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
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp)
target_link_libraries(IPConstitutiveLaw ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase)
add_test(c::unit::IPConstitutiveLaw ${CMAKE_CURRENT_BUILD_DIR}/IPConstitutiveLaw)

add_executable(EngineeringStrain EngineeringStrain.cpp
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/EngineeringStrain.cpp
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/ConstitutiveEnum.cpp
        )
target_link_libraries(EngineeringStrain ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase)
add_test(c::unit::EngineeringStrain ${CMAKE_CURRENT_BUILD_DIR}/EngineeringStrain)

add_executable(EquivalentStrain EquivalentStrain.cpp
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/EquivalentStrain.cpp
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/EngineeringStrain.cpp
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
        ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/ConstitutiveEnum.cpp
        )
target_link_libraries(EquivalentStrain ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} NuToBase)
add_test(c::unit::EquivalentStrain ${CMAKE_CURRENT_BUILD_DIR}/EquivalentStrain)

