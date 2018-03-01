#include <typeinfo>
#include "BoostUnitTest.h"
#include "mechanics/integrationtypes/IntegrationCompanion.h"

#include "math/shapes/Triangle.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"

#include "math/shapes/Hexahedron.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(CreateTriangleIntegration)
{
    auto triangle1ips = CreateIntegrationType(Triangle(), 1);
    BOOST_CHECK(typeid(*triangle1ips) == typeid(IntegrationType2D3NGauss1Ip));

    auto triangle4ips = CreateIntegrationType(Triangle(), 3);
    BOOST_CHECK(typeid(*triangle4ips) == typeid(IntegrationType2D3NGauss4Ip));
}

BOOST_AUTO_TEST_CASE(CreateHexIntegration)
{
    auto brickOrder3 = CreateIntegrationType(Hexahedron(), 3);
    BOOST_CHECK(typeid(*brickOrder3) == typeid(IntegrationTypeTensorProduct<3>));
    BOOST_CHECK((*brickOrder3).GetNumIntegrationPoints() == 3 * 3 * 3);
}
