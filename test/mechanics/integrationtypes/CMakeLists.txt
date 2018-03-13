add_unit_test(IntegrationTypeTensorProduct
    math/Quadrature.cpp
    math/Legendre.cpp
    )
    
add_unit_test(IntegrationCompanion
    mechanics/integrationtypes/IntegrationType3D6NGauss1Ip.cpp
    mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp
    mechanics/integrationtypes/IntegrationTypeTriangle.cpp
    mechanics/integrationtypes/IntegrationTypeTetrahedron.cpp

    math/Quadrature.cpp
    math/Legendre.cpp
    )

add_unit_test(IntegrationTypesTriangle
    mechanics/integrationtypes/IntegrationCompanion.cpp
    mechanics/integrationtypes/IntegrationType3D6NGauss1Ip.cpp
    mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp
    mechanics/integrationtypes/IntegrationTypeTriangle.cpp
    mechanics/integrationtypes/IntegrationTypeTetrahedron.cpp
    math/Quadrature.cpp
    math/Legendre.cpp
    )

add_unit_test(IntegrationTypesTetrahedra
    mechanics/integrationtypes/IntegrationCompanion.cpp
    mechanics/integrationtypes/IntegrationType3D6NGauss1Ip.cpp
    mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.cpp
    mechanics/integrationtypes/IntegrationTypeTensorProduct.cpp
    mechanics/integrationtypes/IntegrationTypeTriangle.cpp
    mechanics/integrationtypes/IntegrationTypeTetrahedron.cpp
    math/Quadrature.cpp
    math/Legendre.cpp
    )
