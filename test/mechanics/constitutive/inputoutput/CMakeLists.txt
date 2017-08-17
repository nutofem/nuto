add_unit_test(EngineeringStrain
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )

add_unit_test(EngineeringStress
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )

add_unit_test(EquivalentStrain
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )

add_unit_test(ConstitutiveIOMap
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )
    
