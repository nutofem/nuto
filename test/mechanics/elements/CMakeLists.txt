add_executable(ContinuumElement ContinuumElement.cpp)
target_link_libraries(ContinuumElement Mechanics Math Base)
if(ENABLE_VISUALIZE)
    target_link_libraries(ContinuumElement Visualize)
endif()
add_test(almost::mechanics::elements::ContinuumElement ${CMAKE_CURRENT_BINARY_DIR}/ContinuumElement)

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
