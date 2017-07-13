#add_unit_test(GradientDamageFatigueEngineeringStress
#    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
#    mechanics/constitutive/inputoutput/EquivalentStrain.cpp
#    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
#    mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp
#    mechanics/constitutive/laws/EngineeringStressHelper.cpp
#    mechanics/constitutive/laws/GradientDamageEngineeringStress.cpp
#    mechanics/constitutive/ConstitutiveBase.cpp
#    mechanics/constitutive/ConstitutiveEnum.cpp
#    mechanics/timeIntegration/ImplExCallback.cpp
#    )
#
add_unit_test(LocalDamageModel
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/inputoutput/EquivalentStrain.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp
    mechanics/constitutive/laws/EngineeringStressHelper.cpp
    mechanics/constitutive/ConstitutiveBase.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )

add_unit_test(HeatConduction
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp
    mechanics/constitutive/ConstitutiveBase.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )

add_unit_test(LinearElasticAnisotropic
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp
    mechanics/constitutive/ConstitutiveBase.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )

add_unit_test(PhaseField
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/laws/EngineeringStressHelper.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOMap.cpp
    mechanics/constitutive/ConstitutiveBase.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    )
