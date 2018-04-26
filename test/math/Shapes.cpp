#include "BoostUnitTest.h"
#include "nuto/math/shapes/Line.h"
#include "nuto/math/shapes/Triangle.h"
#include "nuto/math/shapes/Quadrilateral.h"
#include "nuto/math/shapes/Tetrahedron.h"
#include "nuto/math/shapes/Hexahedron.h"
#include "nuto/math/shapes/Prism.h"
#include "nuto/math/shapes/Pyramid.h"

#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronLinear.h"
#include "nuto/mechanics/interpolation/InterpolationPrismLinear.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(LineIsInside)
{
    Line shape;

    Eigen::VectorXd pIn = Eigen::VectorXd::Constant(1, -0.1);
    Eigen::VectorXd pOut = Eigen::VectorXd::Constant(1, -5.1);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(InterpolationTrussLinear::NodeCoordinatesTrussOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationTrussLinear::NodeCoordinatesTrussOrder1(1)));

    BoostUnitTest::CheckOutstream(shape, "Line");
}

BOOST_AUTO_TEST_CASE(TriangleIsInside)
{
    Triangle shape;

    Eigen::Vector2d pIn(0.1, 0.1);
    Eigen::Vector2d pOut(1., 1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(InterpolationTriangleLinear::NodeCoordinatesTriangleOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationTriangleLinear::NodeCoordinatesTriangleOrder1(1)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationTriangleLinear::NodeCoordinatesTriangleOrder1(2)));

    BoostUnitTest::CheckOutstream(shape, "Triangle");
}

BOOST_AUTO_TEST_CASE(QuadIsInside)
{
    Quadrilateral shape;

    Eigen::Vector2d pIn(0.1, -0.1);
    Eigen::Vector2d pOut(1.2, 1.);
    Eigen::Vector2d pBoundary(-1., -1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(shape.IsWithinShape(pBoundary));

    BoostUnitTest::CheckOutstream(shape, "Quadrilateral");
}

BOOST_AUTO_TEST_CASE(TetrahedronIsInside)
{
    Tetrahedron shape;

    Eigen::Vector3d pIn(0.1, 0.1, 0.1);
    Eigen::Vector3d pOut(0.1, 0.1, -0.1);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(InterpolationTetrahedronLinear::NodeCoordinatesTetrahedronOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationTetrahedronLinear::NodeCoordinatesTetrahedronOrder1(1)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationTetrahedronLinear::NodeCoordinatesTetrahedronOrder1(2)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationTetrahedronLinear::NodeCoordinatesTetrahedronOrder1(3)));

    BoostUnitTest::CheckOutstream(shape, "Tetrahedron");
}

BOOST_AUTO_TEST_CASE(HexahedronIsInside)
{
    Hexahedron shape;

    Eigen::Vector3d pIn(0.1, 0.1, -0.1);
    Eigen::Vector3d pOut(-1.1, 0.1, -0.1);
    Eigen::Vector3d pBoundary(-1., 1., 1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(shape.IsWithinShape(pBoundary));

    BoostUnitTest::CheckOutstream(shape, "Hexahedron");
}

BOOST_AUTO_TEST_CASE(PrismIsInside)
{
    Prism shape;

    Eigen::Vector3d pIn(0.1, 0.1, -0.1);
    Eigen::Vector3d pOut(0.6, 0.6, -0.1);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(InterpolationPrismLinear::NodeCoordinatesPrismOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationPrismLinear::NodeCoordinatesPrismOrder1(1)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationPrismLinear::NodeCoordinatesPrismOrder1(2)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationPrismLinear::NodeCoordinatesPrismOrder1(3)));
    BOOST_CHECK(shape.IsWithinShape(InterpolationPrismLinear::NodeCoordinatesPrismOrder1(4)));

    BoostUnitTest::CheckOutstream(shape, "Prism");
}

BOOST_AUTO_TEST_CASE(PyramidIsInside)
{
    Pyramid shape;

    Eigen::Vector3d pIn(0.2, 0.2, 0.2);
    Eigen::Vector3d pOut(0.2, 0.2, 0.9);
    Eigen::Vector3d pBoundary(0.5, 0.5, 0.5);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(shape.IsWithinShape(pBoundary));

    BoostUnitTest::CheckOutstream(shape, "Pyramid");
}
