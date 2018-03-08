#include "BoostUnitTest.h"
#include "math/shapes/Line.h"
#include "math/shapes/Triangle.h"
#include "math/shapes/Quad.h"
#include "math/shapes/Tetrahedron.h"
#include "math/shapes/Hexahedron.h"
#include "math/shapes/Prism.h"
#include "math/shapes/Pyramid.h"

BOOST_AUTO_TEST_CASE(LineIsInside)
{
    NuTo::Line shape;

    Eigen::VectorXd pIn = Eigen::VectorXd::Constant(1, 0.1);
    Eigen::VectorXd pOut = Eigen::VectorXd::Constant(1, -5.1);
    Eigen::VectorXd pBoundary = Eigen::VectorXd::Constant(1, -1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(!shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(TriangleIsInside)
{
    NuTo::Triangle shape;

    Eigen::Vector2d pIn(0.1, 0.1);
    Eigen::Vector2d pOut(1., 1.);
    Eigen::Vector2d pBoundary(0., 0.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(!shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(QuadIsInside)
{
    NuTo::Quad shape;

    Eigen::Vector2d pIn(0.1, 0.1);
    Eigen::Vector2d pOut(1.2, 1.);
    Eigen::Vector2d pBoundary(-1., -1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(!shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(TetrahedronIsInside)
{
    NuTo::Tetrahedron shape;

    Eigen::Vector3d pIn(0.1, 0.1, 0.1);
    Eigen::Vector3d pOut(0.1, 0.1, -0.1);
    Eigen::Vector3d pBoundary(0.5, 0.5, 0.5);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(!shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(HexahedronIsInside)
{
    NuTo::Hexahedron shape;

    Eigen::Vector3d pIn(0.1, 0.1, 0.1);
    Eigen::Vector3d pOut(-1.1, 0.1, -0.1);
    Eigen::Vector3d pBoundary(-1., 1., 1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(!shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(PrismIsInside)
{
    NuTo::Prism shape;

    Eigen::Vector3d pIn(0.1, 0.1, 0.1);
    Eigen::Vector3d pOut(0.1, 0.5, 0.6);
    Eigen::Vector3d pBoundary(-1., -1., 0.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(!shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(PyramidIsInside)
{
    NuTo::Pyramid shape;

    Eigen::Vector3d pIn(0.2, 0.2, 0.2);
    Eigen::Vector3d pOut(0.2, 0.2, 0.9);
    Eigen::Vector3d pBoundary(0.5, 0.5, 0.5);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(!shape.IsWithinShape(pBoundary));
}
