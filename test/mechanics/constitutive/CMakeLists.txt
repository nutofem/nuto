foreach(test AdditiveInput AdditiveOutput) 
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} Mechanics Math Base)
    if(ENABLE_VISUALIZE)
        target_link_libraries(${test} Visualize)
    endif()
    add_test(almost::mechanics::constitutive::${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()

add_unit_test(DataContainer
    base/serializeStream/SerializeStreamBase.cpp
    base/serializeStream/SerializeStreamIn.cpp
    base/serializeStream/SerializeStreamOut.cpp
    )


add_unit_test(IPConstitutiveLaw
    mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp
    )

add_unit_test(EngineeringStrainInvariants)
add_unit_test(EngineeringStrainPDE)
add_unit_test(EngineeringStressPDE)

add_subdirectory(damageLaws)
add_subdirectory(inputoutput)
add_subdirectory(laws)
