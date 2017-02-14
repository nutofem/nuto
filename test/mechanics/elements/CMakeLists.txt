add_executable(ContinuumElement ContinuumElement.cpp)
target_link_libraries(ContinuumElement NuToMechanics NuToMath NuToBase)
if(ENABLE_VISUALIZE)
    target_link_libraries(ContinuumElement NuToVisualize)
endif()
add_test(unit::mechanics::elements::ContinuumElement ${CMAKE_CURRENT_BINARY_DIR}/ContinuumElement)

add_unit_test(IPData
    mechanics/constitutive/ConstitutiveBase.cpp
    mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    mechanics/constitutive/laws/EngineeringStressHelper.cpp
    mechanics/constitutive/inputoutput/EquivalentStrain.cpp
    mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    mechanics/constitutive/laws/GradientDamageEngineeringStress.cpp
    mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.cpp
    mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.cpp
    mechanics/integrationtypes/IntegrationTypeBase.cpp
    mechanics/integrationtypes/IntegrationType1D.cpp
    mechanics/integrationtypes/IntegrationTypeEnum.cpp
    mechanics/constitutive/ConstitutiveEnum.cpp
    mechanics/timeIntegration/ImplExCallback.cpp
    )
