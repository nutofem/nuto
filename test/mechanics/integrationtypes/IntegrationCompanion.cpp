#include <typeinfo>
#include "BoostUnitTest.h"
#include "mechanics/integrationtypes/IntegrationCompanion.h"

#include "math/shapes/Triangle.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"

#include "math/shapes/Hexahedron.h"
#include "math/shapes/Quad.h"
#include "math/shapes/Line.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(CreateTriangleIntegration)
{
    auto triangle1ips = CreateGaussIntegrationType(Triangle(), 1);
    BOOST_CHECK(typeid(*triangle1ips) == typeid(IntegrationType2D3NGauss1Ip));

    auto triangle4ips = CreateGaussIntegrationType(Triangle(), 3);
    BOOST_CHECK(typeid(*triangle4ips) == typeid(IntegrationType2D3NGauss4Ip));
}

BOOST_AUTO_TEST_CASE(CreateLineIntegrationGauss)
{
    auto lineOrder3 = CreateGaussIntegrationType(Line(), 3);
    BOOST_CHECK(typeid(*lineOrder3) == typeid(IntegrationTypeTensorProduct<1>));
    BOOST_CHECK((*lineOrder3).GetNumIntegrationPoints() == 3);
}

BOOST_AUTO_TEST_CASE(CreateQuadIntegrationGauss)
{
    auto quadOrder3 = CreateGaussIntegrationType(Quad(), 3);
    BOOST_CHECK(typeid(*quadOrder3) == typeid(IntegrationTypeTensorProduct<2>));
    BOOST_CHECK((*quadOrder3).GetNumIntegrationPoints() == 3 * 3);
}

BOOST_AUTO_TEST_CASE(CreateHexIntegrationGauss)
{
    auto brickOrder3 = CreateGaussIntegrationType(Hexahedron(), 3);
    BOOST_CHECK(typeid(*brickOrder3) == typeid(IntegrationTypeTensorProduct<3>));
    BOOST_CHECK((*brickOrder3).GetNumIntegrationPoints() == 3 * 3 * 3);
}

BOOST_AUTO_TEST_CASE(CreateLineIntegrationLobatto)
{
    auto lineOrder3 = CreateLobattoIntegrationType(Line(), 3);
    BOOST_CHECK(typeid(*lineOrder3) == typeid(IntegrationTypeTensorProduct<1>));
    BOOST_CHECK((*lineOrder3).GetNumIntegrationPoints() == 3);
    BOOST_CHECK((*lineOrder3).GetNumIntegrationPoints() == 3);
}

BOOST_AUTO_TEST_CASE(CreateQuadIntegrationLobatto)
{
    auto quadOrder3 = CreateLobattoIntegrationType(Quad(), 3);
    BOOST_CHECK(typeid(*quadOrder3) == typeid(IntegrationTypeTensorProduct<2>));
    BOOST_CHECK((*quadOrder3).GetNumIntegrationPoints() == 3 * 3);
}

BOOST_AUTO_TEST_CASE(CreateTriangleIntegrationLobatto)
{
    BOOST_CHECK_THROW(auto triOrder3 = CreateLobattoIntegrationType(Triangle(), 3), Exception);
}

BOOST_AUTO_TEST_CASE(CreateHexIntegrationLobatto)
{
    auto brickOrder3 = CreateLobattoIntegrationType(Hexahedron(), 3);
    BOOST_CHECK(typeid(*brickOrder3) == typeid(IntegrationTypeTensorProduct<3>));
    BOOST_CHECK((*brickOrder3).GetNumIntegrationPoints() == 3 * 3 * 3);
}
