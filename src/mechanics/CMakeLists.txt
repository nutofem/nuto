set(MechanicsSources
    constraints/Constraints.cpp
    constraints/ConstraintCompanion.cpp
    dofs/DofNumbering.cpp
    elements/ElementShapeFunctions.cpp
    mesh/MeshFemDofConvert.cpp
    mesh/UnitMeshFem.cpp
    timeIntegration/TimeControl.cpp
    tools/GlobalFractureEnergyIntegrator.cpp
    mesh/MeshFem.cpp
    integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.cpp
    integrationtypes/IntegrationType0DBoundary.cpp
    integrationtypes/IntegrationType2D3NGauss13Ip.cpp
    integrationtypes/IntegrationType2D3NGauss16Ip.cpp
    integrationtypes/IntegrationType2D3NGauss1Ip.cpp
    integrationtypes/IntegrationType2D3NGauss3Ip.cpp
    integrationtypes/IntegrationType2D3NGauss4Ip.cpp
    integrationtypes/IntegrationType2D3NGauss6Ip.cpp
    integrationtypes/IntegrationType2D3NGauss12Ip.cpp
    integrationtypes/IntegrationType2D3NGauss12IpDetail.cpp
    integrationtypes/IntegrationType3D4NGauss1Ip.cpp
    integrationtypes/IntegrationType3D4NGauss4Ip.cpp
    integrationtypes/IntegrationType3D6NGauss1Ip.cpp
    integrationtypes/IntegrationType3D6NGauss2x3Ip.cpp
    integrationtypes/IntegrationTypeTensorProduct.cpp
    integrationtypes/IntegrationTypeBase.cpp
    integrationtypes/IntegrationTypeEnum.cpp
    )

create_nuto_module(Mechanics "${MechanicsSources}")
target_link_libraries(Mechanics Base Math Boost::filesystem Ann::Ann)
