add_unit_test(EngineeringStrain
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/ConstitutiveEnum.cpp
    )

add_unit_test(EquivalentStrain
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/ConstitutiveEnum.cpp
    )
