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

        Eigen::MatrixXd B0 = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(
                Eigen::VectorXd::Constant(1, -std::sqrt(1. / 3.)));
        Eigen::MatrixXd B1 = NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(
                Eigen::VectorXd::Constant(1, std::sqrt(1. / 3.)));
        NuTo::Jacobian jacobian0(coordinates, B0, 2);
        NuTo::Jacobian jacobian1(coordinates, B1, 2);
        BOOST_CHECK_CLOSE(jacobian0.Det() + jacobian1.Det(), 5, 1.e-10);
        //
        // total length of the element is 5 = sqrt(3*3 + 4*4)
        // each integration point has weight 1
    }
}

//! @brief calculate the area of a triangle (points a,b,c) element based on ||J|| and check vs. `correctArea`
void CheckJacobian2Din3D(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, double correctArea)
{
    Eigen::MatrixXd B = NuTo::ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder1(Eigen::Vector2d::Zero());
    Eigen::VectorXd coordinates = Eigen::VectorXd(9);
    coordinates << a.x(), a.y(), a.z(), b.x(), b.y(), b.z(), c.x(), c.y(), c.z();

    NuTo::Jacobian jacobian(coordinates, B, 3);

    double ratioAreaRealToAreaNatural = jacobian.Det();
    double areaNatural = 0.5; // triangle (0,0) ; (1,0) ; (0,1)
    double calculatedArea = ratioAreaRealToAreaNatural * areaNatural;

    BOOST_CHECK_CLOSE(correctArea, calculatedArea, 1.e-10);
}

BOOST_AUTO_TEST_CASE(Jacobian2Din3D)
{
    // some manual test cases
    CheckJacobian2Din3D(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), 0.5);
    CheckJacobian2Din3D(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(5, 0, 0), Eigen::Vector3d(0, 0, 7), 5. * 7 / 2);
    CheckJacobian2Din3D(Eigen::Vector3d(0, 4, 0), Eigen::Vector3d(0, 0, 3), Eigen::Vector3d(0, 4, 3), 3. * 4 / 2);

    // some random test cases
    for (int i = 0; i < 100; ++i)
    {
        // arbitrary 3d points A, B, C define a triangle (A,B,C)
        Eigen::Vector3d a = Eigen::Vector3d::Random();
        Eigen::Vector3d b = Eigen::Vector3d::Random();
        Eigen::Vector3d c = Eigen::Vector3d::Random();

        double analyticalArea = 0.5 * (b - a).cross(c - a).norm(); // wiki...
        CheckJacobian2Din3D(a, b, c, analyticalArea);
    }
}

BOOST_AUTO_TEST_CASE(JacobianTransform)
{
    // TODO! How?
}
