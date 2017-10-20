#include "BoostUnitTest.h"
#include "mechanics/cell/Jacobian.h"
#include "mechanics/elements/ElementShapeFunctions.h"

BOOST_AUTO_TEST_CASE(Jacobian1DDet)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1();
    Eigen::VectorXd coordinates = Eigen::Vector2d({12., 17});
    NuTo::Jacobian jacobian(coordinates, B, 1);
    BOOST_CHECK_CLOSE(jacobian.Det(), 2.5, 1.e-10);
}

BOOST_AUTO_TEST_CASE(Jacobian2DDet)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(Eigen::Vector2d({0, 0}));
    Eigen::VectorXd coordinates = Eigen::VectorXd(8);
    coordinates << 0, 0, 2, 0, 2, 8, 0, 8; // rectangle from (0,0) to (2,8)


    NuTo::Jacobian jacobian(coordinates, B, 2);
    BOOST_CHECK_CLOSE(jacobian.Det(), 4, 1.e-10); // reference area = 4, this area = 16;
}

BOOST_AUTO_TEST_CASE(Jacobian3DDet)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder1();
    Eigen::VectorXd coordinates = Eigen::VectorXd(12);
    coordinates << 0, 0, 0, 10, 0, 0, 0, 5, 0, 0, 0, 42;

    NuTo::Jacobian jacobian(coordinates, B, 3);
    BOOST_CHECK_CLOSE(jacobian.Det(), 10 * 5 * 42, 1.e-10);
}

BOOST_AUTO_TEST_CASE(Jacobian1Din2D)
{
    {
        Eigen::MatrixXd B = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1();
        Eigen::VectorXd coordinates = Eigen::VectorXd(4);
        coordinates << 0, 0, 3, 4;

        NuTo::Jacobian jacobian(coordinates, B, 2);
        BOOST_CHECK_CLOSE(jacobian.Det() * 2, 5, 1.e-10);
        //
        // total length of the element is 5 = sqrt(3*3 + 4*4)
        // factor 2 in Det() comes from the integration point weigth
    }
    {
        Eigen::VectorXd coordinates = Eigen::VectorXd(6);
        coordinates << 0, 0, 1.5, 2, 3, 4;

        Eigen::MatrixXd B0 = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(Eigen::VectorXd::Constant(1, -std::sqrt(1./3.)));
        Eigen::MatrixXd B1 = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(Eigen::VectorXd::Constant(1, std::sqrt(1./3.)));
        NuTo::Jacobian jacobian0(coordinates, B0, 2);
        NuTo::Jacobian jacobian1(coordinates, B1, 2);
        BOOST_CHECK_CLOSE(jacobian0.Det() + jacobian1.Det(), 5, 1.e-10);
        //
        // total length of the element is 5 = sqrt(3*3 + 4*4)
        // each integration point has weight 1
    }
}

BOOST_AUTO_TEST_CASE(JacobianTransform)
{
    // TODO! How?
}
