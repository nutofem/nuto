#include "BoostUnitTest.h"
#include "nuto/math/shapes/Line.h"
#include "nuto/math/shapes/Triangle.h"
#include "nuto/math/shapes/Quadrilateral.h"
#include "nuto/math/shapes/Tetrahedron.h"
#include "nuto/math/shapes/Hexahedron.h"
#include "nuto/math/shapes/Prism.h"
#include "nuto/math/shapes/Pyramid.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"


BOOST_AUTO_TEST_CASE(LineIsInside)
{
    NuTo::Line shape;

    Eigen::VectorXd pIn = Eigen::VectorXd::Constant(1, -0.1);
    Eigen::VectorXd pOut = Eigen::VectorXd::Constant(1, -5.1);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions1D::NodeCoordinatesTrussOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions1D::NodeCoordinatesTrussOrder1(1)));
}

BOOST_AUTO_TEST_CASE(TriangleIsInside)
{
    NuTo::Triangle shape;

    Eigen::Vector2d pIn(0.1, 0.1);
    Eigen::Vector2d pOut(1., 1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions2D::NodeCoordinatesTriangleOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions2D::NodeCoordinatesTriangleOrder1(1)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions2D::NodeCoordinatesTriangleOrder1(2)));
}

BOOST_AUTO_TEST_CASE(QuadIsInside)
{
    NuTo::Quadrilateral shape;

    Eigen::Vector2d pIn(0.1, -0.1);
    Eigen::Vector2d pOut(1.2, 1.);
    Eigen::Vector2d pBoundary(-1., -1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(TetrahedronIsInside)
{
    NuTo::Tetrahedron shape;

    Eigen::Vector3d pIn(0.1, 0.1, 0.1);
    Eigen::Vector3d pOut(0.1, 0.1, -0.1);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesTetrahedronOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesTetrahedronOrder1(1)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesTetrahedronOrder1(2)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesTetrahedronOrder1(3)));
}

BOOST_AUTO_TEST_CASE(HexahedronIsInside)
{
    NuTo::Hexahedron shape;

    Eigen::Vector3d pIn(0.1, 0.1, -0.1);
    Eigen::Vector3d pOut(-1.1, 0.1, -0.1);
    Eigen::Vector3d pBoundary(-1., 1., 1.);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(shape.IsWithinShape(pBoundary));
}

BOOST_AUTO_TEST_CASE(PrismIsInside)
{
    NuTo::Prism shape;

    Eigen::Vector3d pIn(0.1, 0.1, -0.1);
    Eigen::Vector3d pOut(0.6, 0.6, -0.1);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));

    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesPrismOrder1(0)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesPrismOrder1(1)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesPrismOrder1(2)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesPrismOrder1(3)));
    BOOST_CHECK(shape.IsWithinShape(NuTo::ShapeFunctions3D::NodeCoordinatesPrismOrder1(4)));
}

BOOST_AUTO_TEST_CASE(PyramidIsInside)
{
    NuTo::Pyramid shape;

    Eigen::Vector3d pIn(0.2, 0.2, 0.2);
    Eigen::Vector3d pOut(0.2, 0.2, 0.9);
    Eigen::Vector3d pBoundary(0.5, 0.5, 0.5);

    BOOST_CHECK(shape.IsWithinShape(pIn));
    BOOST_CHECK(!shape.IsWithinShape(pOut));
    BOOST_CHECK(shape.IsWithinShape(pBoundary));
}
