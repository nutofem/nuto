#include "BoostUnitTest.h"
#include "nuto/mechanics/integrationtypes/IntegrationCompanion.h"

#include "nuto/math/shapes/Triangle.h"

#include "nuto/math/shapes/Hexahedron.h"
#include "nuto/math/shapes/Quadrilateral.h"
#include "nuto/math/shapes/Line.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(CreateTriangleIntegration)
{
    auto triangle1ips = CreateGaussIntegrationType(Triangle(), 1);
    BOOST_CHECK_EQUAL(triangle1ips->GetNumIntegrationPoints(), 1);

    auto triangle4ips = CreateGaussIntegrationType(Triangle(), 3);
    BOOST_CHECK_EQUAL(triangle4ips->GetNumIntegrationPoints(), 4);
}

BOOST_AUTO_TEST_CASE(CreateLineIntegrationGauss)
{
    auto lineOrder3 = CreateGaussIntegrationType(Line(), 3);
    BOOST_CHECK_EQUAL(lineOrder3->GetNumIntegrationPoints(), 3);
}

BOOST_AUTO_TEST_CASE(CreateQuadIntegrationGauss)
{
    auto quadOrder3 = CreateGaussIntegrationType(Quadrilateral(), 3);
    BOOST_CHECK_EQUAL(quadOrder3->GetNumIntegrationPoints(), 3 * 3);
}

BOOST_AUTO_TEST_CASE(CreateHexIntegrationGauss)
{
    auto brickOrder3 = CreateGaussIntegrationType(Hexahedron(), 3);
    BOOST_CHECK_EQUAL(brickOrder3->GetNumIntegrationPoints(), 3 * 3 * 3);
}

BOOST_AUTO_TEST_CASE(CreateLineIntegrationLobatto)
{
    auto lineOrder3 = CreateLobattoIntegrationType(Line(), 3);
    BOOST_CHECK_EQUAL(lineOrder3->GetNumIntegrationPoints(), 3);
}

BOOST_AUTO_TEST_CASE(CreateQuadIntegrationLobatto)
{
    auto quadOrder3 = CreateLobattoIntegrationType(Quadrilateral(), 3);
    BOOST_CHECK_EQUAL(quadOrder3->GetNumIntegrationPoints(), 3 * 3);
}

BOOST_AUTO_TEST_CASE(CreateTriangleIntegrationLobatto)
{
    BOOST_CHECK_THROW(auto triOrder3 = CreateLobattoIntegrationType(Triangle(), 3), Exception);
}

BOOST_AUTO_TEST_CASE(CreateHexIntegrationLobatto)
{
    auto brickOrder3 = CreateLobattoIntegrationType(Hexahedron(), 3);
    BOOST_CHECK_EQUAL(brickOrder3->GetNumIntegrationPoints(), 3 * 3 * 3);
}
