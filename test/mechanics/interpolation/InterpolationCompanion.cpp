#include <typeinfo>
#include "BoostUnitTest.h"
#include "mechanics/interpolation/InterpolationCompanion.h"

#include "math/shapes/Triangle.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/interpolation/InterpolationTriangleQuadratic.h"

#include "math/shapes/Hexahedron.h"
#include "mechanics/interpolation/InterpolationBrickLinear.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(CreateTriangleInterpolation)
{
    auto linearTriangle = CreateInterpolation(Triangle(), 1);
    BOOST_CHECK(typeid(*linearTriangle) == typeid(InterpolationTriangleLinear));

    auto quadraticTriangle = CreateInterpolation(Triangle(), 2);
    BOOST_CHECK(typeid(*quadraticTriangle) == typeid(InterpolationTriangleQuadratic));
}

BOOST_AUTO_TEST_CASE(CreateHexInterpolation)
{
    auto linearBrick = CreateInterpolation(Hexahedron(), 1);
    BOOST_CHECK(typeid(*linearBrick) == typeid(InterpolationBrickLinear));
}
