#include "BoostUnitTest.h"
#include "mechanics/cell/GradientOperators.h"

constexpr double dN0dX = 1;
constexpr double dN0dY = 2;
constexpr double dN0dZ = 3;
constexpr double dN1dX = 11;
constexpr double dN1dY = 12;
constexpr double dN1dZ = 13;
constexpr double dN2dX = 21;
constexpr double dN2dY = 22;
constexpr double dN2dZ = 23;

BOOST_AUTO_TEST_CASE(Gradients1D)
{
    Eigen::MatrixXd dNdX(3, 1);
    dNdX(0, 0) = dN0dX;
    dNdX(1, 0) = dN1dX;
    dNdX(2, 0) = dN2dX;

    BoostUnitTest::CheckEigenMatrix(dNdX.transpose(), NuTo::Nabla::Gradient()(dNdX));
    BoostUnitTest::CheckEigenMatrix(dNdX.transpose(), NuTo::Nabla::Strain()(dNdX));
}

BOOST_AUTO_TEST_CASE(Gradients2D)
{
    Eigen::MatrixXd dNdX(3, 2);
    dNdX(0, 0) = dN0dX;
    dNdX(0, 1) = dN0dY;
    dNdX(1, 0) = dN1dX;
    dNdX(1, 1) = dN1dY;
    dNdX(2, 0) = dN2dX;
    dNdX(2, 1) = dN2dY;

    Eigen::MatrixXd correctGradient(2, 3);
    correctGradient(0, 0) = dN0dX;
    correctGradient(0, 1) = dN1dX;
    correctGradient(0, 2) = dN2dX;
    correctGradient(1, 0) = dN0dY;
    correctGradient(1, 1) = dN1dY;
    correctGradient(1, 2) = dN2dY;
    BoostUnitTest::CheckEigenMatrix(correctGradient, NuTo::Nabla::Gradient()(dNdX));


    Eigen::MatrixXd correctStrain = Eigen::MatrixXd::Zero(3, 6);
    correctStrain(0, 0) = dN0dX;
    correctStrain(1, 1) = dN0dY;
    correctStrain(2, 0) = dN0dY;
    correctStrain(2, 1) = dN0dX;

    correctStrain(0, 2) = dN1dX;
    correctStrain(1, 3) = dN1dY;
    correctStrain(2, 2) = dN1dY;
    correctStrain(2, 3) = dN1dX;

    correctStrain(0, 4) = dN2dX;
    correctStrain(1, 5) = dN2dY;
    correctStrain(2, 4) = dN2dY;
    correctStrain(2, 5) = dN2dX;

    BoostUnitTest::CheckEigenMatrix(correctStrain, NuTo::Nabla::Strain()(dNdX));
}

BOOST_AUTO_TEST_CASE(Gradients3D)
{
    Eigen::MatrixXd dNdX(3, 3);
    dNdX(0, 0) = dN0dX;
    dNdX(0, 1) = dN0dY;
    dNdX(0, 2) = dN0dZ;
    dNdX(1, 0) = dN1dX;
    dNdX(1, 1) = dN1dY;
    dNdX(1, 2) = dN1dZ;
    dNdX(2, 0) = dN2dX;
    dNdX(2, 1) = dN2dY;
    dNdX(2, 2) = dN2dZ;

    Eigen::MatrixXd correctStrain = Eigen::MatrixXd::Zero(6, 9);
    correctStrain(0, 0) = dN0dX;
    correctStrain(1, 1) = dN0dY;
    correctStrain(2, 2) = dN0dZ;
    correctStrain(3, 1) = dN0dZ;
    correctStrain(3, 2) = dN0dY;
    correctStrain(4, 0) = dN0dZ;
    correctStrain(4, 2) = dN0dX;
    correctStrain(5, 0) = dN0dY;
    correctStrain(5, 1) = dN0dX;

    correctStrain(0, 3) = dN1dX;
    correctStrain(1, 4) = dN1dY;
    correctStrain(2, 5) = dN1dZ;
    correctStrain(3, 4) = dN1dZ;
    correctStrain(3, 5) = dN1dY;
    correctStrain(4, 3) = dN1dZ;
    correctStrain(4, 5) = dN1dX;
    correctStrain(5, 3) = dN1dY;
    correctStrain(5, 4) = dN1dX;

    correctStrain(0, 6) = dN2dX;
    correctStrain(1, 7) = dN2dY;
    correctStrain(2, 8) = dN2dZ;
    correctStrain(3, 7) = dN2dZ;
    correctStrain(3, 8) = dN2dY;
    correctStrain(4, 6) = dN2dZ;
    correctStrain(4, 8) = dN2dX;
    correctStrain(5, 6) = dN2dY;
    correctStrain(5, 7) = dN2dX;

    BoostUnitTest::CheckEigenMatrix(correctStrain, NuTo::Nabla::Strain()(dNdX));
}
