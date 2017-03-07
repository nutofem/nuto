#include "BoostUnitTest.h"
#include "mechanics/cell/Jacobian.h"
#include "mechanics/elements/ElementShapeFunctions.h"

BOOST_AUTO_TEST_CASE(Jacobian1DDet)
{
    Eigen::MatrixXd B           = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1(Eigen::VectorXd());
    Eigen::VectorXd coordinates = Eigen::Vector2d({12., 17});
    NuTo::Jacobian<1> jacobian(coordinates, B);
    BOOST_CHECK_CLOSE(jacobian.Det(), 2.5, 1.e-10);
}

BOOST_AUTO_TEST_CASE(Jacobian2DDet)
{
    Eigen::MatrixXd B           = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(Eigen::Vector2d({0, 0}));
    Eigen::VectorXd coordinates = Eigen::VectorXd(8);
    coordinates << 0, 0, 2, 0, 2, 8, 0, 8; // rectangle from (0,0) to (2,8)


    NuTo::Jacobian<2> jacobian(coordinates, B);
    BOOST_CHECK_CLOSE(jacobian.Det(), 4, 1.e-10); // reference area = 4, this area = 16;
}

BOOST_AUTO_TEST_CASE(Jacobian3DDet)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder1(Eigen::Vector3d({0, 0, 0}));
    Eigen::VectorXd coordinates = Eigen::VectorXd(12);
    coordinates << 0, 0, 0, 10, 0, 0, 0, 5, 0, 0, 0, 42;

    NuTo::Jacobian<3> jacobian(coordinates, B);
    BOOST_CHECK_CLOSE(jacobian.Det(), 10 * 5 * 42, 1.e-10);
}

BOOST_AUTO_TEST_CASE(JacobianTransform)
{
    // TODO! How?
}
