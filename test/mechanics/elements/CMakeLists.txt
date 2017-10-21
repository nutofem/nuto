add_unit_test(IPData
    mechanics/constitutive/ConstitutiveBase.cpp
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/laws/EngineeringStressHelper.cpp
    mechanics/constitutive/inputoutput/EquivalentStrain.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/laws/GradientDamageEngineeringStress.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp
    math/Legendre.cpp
    mechanics/integrationtypes/IntegrationTypeBase.cpp
    mechanics/integrationtypes/IntegrationTypeEnum.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    mechanics/timeIntegration/ImplExCallback.cpp
    )

add_unit_test(ElementIga)
add_unit_test(ElementFem mechanics/elements/ElementShapeFunctions.cpp)
add_unit_test(ElementCollection mechanics/elements/ElementShapeFunctions.cpp)
