add_executable(ContinuumElement ContinuumElement.cpp)
target_link_libraries(ContinuumElement NuToMechanics NuToMath NuToBase)
if(ENABLE_VISUALIZE)
    target_link_libraries(ContinuumElement NuToVisualize)
endif()
add_test(unit::mechanics::elements::ContinuumElement ${CMAKE_CURRENT_BINARY_DIR}/ContinuumElement)

add_unit_test(IPData
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/ConstitutiveBase.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/EngineeringStrain.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/laws/EngineeringStressHelper.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/EquivalentStrain.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/inputoutput/ConstitutiveIOBase.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/laws/GradientDamageEngineeringStress.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/integrationtypes/IntegrationTypeBase.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/integrationtypes/IntegrationType1D.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/integrationtypes/IntegrationTypeEnum.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/constitutive/ConstitutiveEnum.cpp
    ${CMAKE_SOURCE_DIR}/src/mechanics/timeIntegration/ImplExCallback.cpp
    )
