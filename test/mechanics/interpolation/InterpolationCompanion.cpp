#include "BoostUnitTest.h"
#include "nuto/mechanics/interpolation/InterpolationCompanion.h"

#include "nuto/math/shapes/Line.h"
#include "nuto/math/shapes/Triangle.h"
#include "nuto/math/shapes/Quadrilateral.h"
#include "nuto/math/shapes/Hexahedron.h"
#include "nuto/math/shapes/Tetrahedron.h"
#include "nuto/math/shapes/Prism.h"
#include "nuto/math/shapes/Pyramid.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(CreateLineInterpolation)
{
    auto line = CreateLagrangeInterpolation(Line(), 1);
    BOOST_CHECK_EQUAL(line->GetNumNodes(), 2);
}

BOOST_AUTO_TEST_CASE(CreateTriangleInterpolation)
{
    auto linearTriangle = CreateLagrangeInterpolation(Triangle(), 1);
    BOOST_CHECK_EQUAL(linearTriangle->GetNumNodes(), 3);

    auto quadraticTriangle = CreateLagrangeInterpolation(Triangle(), 2);
    BOOST_CHECK_EQUAL(quadraticTriangle->GetNumNodes(), 6);
}

BOOST_AUTO_TEST_CASE(CreateQuadInterpolation)
{
    auto linearQuad = CreateLagrangeInterpolation(Quadrilateral(), 1);
    BOOST_CHECK_EQUAL(linearQuad->GetNumNodes(), 4);

    auto quadraticQuad = CreateLagrangeInterpolation(Quadrilateral(), 2);
    BOOST_CHECK_EQUAL(quadraticQuad->GetNumNodes(), 8);
}

BOOST_AUTO_TEST_CASE(CreateTetInterpolation)
{
    auto linearTet = CreateLagrangeInterpolation(Tetrahedron(), 1);
    BOOST_CHECK_EQUAL(linearTet->GetNumNodes(), 4);

    auto quadraticTet = CreateLagrangeInterpolation(Tetrahedron(), 2);
    BOOST_CHECK_EQUAL(quadraticTet->GetNumNodes(), 10);
}

BOOST_AUTO_TEST_CASE(CreateHexInterpolation)
{
    auto linearBrick = CreateLagrangeInterpolation(Hexahedron(), 1);
    BOOST_CHECK_EQUAL(linearBrick->GetNumNodes(), 8);
}

BOOST_AUTO_TEST_CASE(CreatePrismInterpolation)
{
    auto linearPrism = CreateLagrangeInterpolation(Prism(), 1);
    BOOST_CHECK_EQUAL(linearPrism->GetNumNodes(), 6);

    auto quadraticPrism = CreateLagrangeInterpolation(Prism(), 2);
    BOOST_CHECK_EQUAL(quadraticPrism->GetNumNodes(), 18);
}

BOOST_AUTO_TEST_CASE(CreatePyramidInterpolation)
{
    auto linearPyramid = CreateLagrangeInterpolation(Pyramid(), 1);
    BOOST_CHECK_EQUAL(linearPyramid->GetNumNodes(), 5);
}
