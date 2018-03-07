#include "BoostUnitTest.h"
#include "mechanics/interpolation/InterpolationCompanion.h"

#include "math/shapes/Triangle.h"
#include "math/shapes/Hexahedron.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(CreateTriangleInterpolation)
{
    auto linearTriangle = CreateSerendipityInterpolation(Triangle(), 1);
    BOOST_CHECK_EQUAL(linearTriangle->GetNumNodes(), 3);

    auto quadraticTriangle = CreateSerendipityInterpolation(Triangle(), 2);
    BOOST_CHECK_EQUAL(quadraticTriangle->GetNumNodes(), 6);
}

BOOST_AUTO_TEST_CASE(CreateHexInterpolation)
{
    auto linearBrick = CreateSerendipityInterpolation(Hexahedron(), 1);
    BOOST_CHECK_EQUAL(linearBrick->GetNumNodes(), 8);
}
