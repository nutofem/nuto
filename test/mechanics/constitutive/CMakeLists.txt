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
