foreach(test HeatConduction ConstitutiveIO AdditiveInput AdditiveOutput GradientDamageFatigue LocalDamageModel)
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} NuToMechanics NuToMath NuToBase)
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} NuToVisualize)
    endif()
    add_test(c::unit::${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()

add_unit_test(DataContainer)

add_unit_test(IPConstitutiveLaw
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp
    )

add_subdirectory(inputoutput)
