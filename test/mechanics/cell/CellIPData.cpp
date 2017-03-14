#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/cell/CellIPData.h"
#include <iostream>


constexpr double dN0dX = 1;
constexpr double dN0dY = 2;
constexpr double dN0dZ = 3;
constexpr double dN1dX = 11;
constexpr double dN1dY = 12;
constexpr double dN1dZ = 13;
constexpr double dN2dX = 21;
constexpr double dN2dY = 22;
constexpr double dN2dZ = 23;

NuTo::DerivativeShapeFunctionsNatural MockDerivatives1D()
{
    NuTo::DerivativeShapeFunctionsNatural m(3, 1);
    m(0, 0) = dN0dX;
    m(1, 0) = dN1dX;
    m(2, 0) = dN2dX;
    return m;
}

NuTo::DerivativeShapeFunctionsNatural MockDerivatives2D()
{
    NuTo::DerivativeShapeFunctionsNatural m(3, 2);
    m(0, 0) = dN0dX;
    m(0, 1) = dN0dY;
    m(1, 0) = dN1dX;
    m(1, 1) = dN1dY;
    m(2, 0) = dN2dX;
    m(2, 1) = dN2dY;
    return m;
}

NuTo::DerivativeShapeFunctionsNatural MockDerivatives3D()
{
    NuTo::DerivativeShapeFunctionsNatural m(3, 2);
    m(0, 0) = dN0dX;
    m(0, 1) = dN0dY;
    m(0, 2) = dN0dZ;
    m(1, 0) = dN1dX;
    m(1, 1) = dN1dY;
    m(1, 2) = dN1dZ;
    m(2, 0) = dN2dX;
    m(2, 1) = dN2dY;
    m(2, 2) = dN2dZ;
    return m;
}

BOOST_AUTO_TEST_CASE(CellIPData2D)
{
    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({1, 0}));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({0, 1}));
    NuTo::InterpolationTriangleLinear interpolation0(22);
    NuTo::ElementSimple e0({&n0, &n1, &n2}, interpolation0);
    NuTo::DofType d0("dof0", 2, 0);

    NuTo::NodeSimple n3 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(1));
    NuTo::NodeSimple n4 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(3));
    NuTo::NodeSimple n5 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(7));
    fakeit::Mock<NuTo::InterpolationSimple> interpolation1;
    Method(interpolation1, GetDerivativeShapeFunctions) = MockDerivatives2D();
    NuTo::ElementSimple e1({&n3, &n4, &n5}, interpolation1.get());
    NuTo::DofType d1("dof0", 1, 1);

    NuTo::DofContainer<NuTo::ElementSimple*> elements;
    elements[d0] = &e0;
    elements[d1] = &e1;

    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({1. / 3., 1. / 3.});

    NuTo::Jacobian<2> jac(e0.ExtractNodeValues(), e0.GetInterpolation().GetDerivativeShapeFunctions(ipCoords));

    NuTo::CellIPData<2> ipData(elements, jac, ipCoords);

    BoostUnitTest::CheckEigenMatrix(jac.Inv(), Eigen::Matrix2d::Identity());

    Eigen::MatrixXd correctGradient(2, 3);
    correctGradient(0, 0) = dN0dX;
    correctGradient(0, 1) = dN1dX;
    correctGradient(0, 2) = dN2dX;
    correctGradient(1, 0) = dN0dY;
    correctGradient(1, 1) = dN1dY;
    correctGradient(1, 2) = dN2dY;
    BoostUnitTest::CheckEigenMatrix(ipData.GetBMatrixGradient(d1), correctGradient);


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
    BoostUnitTest::CheckEigenMatrix(ipData.GetBMatrixStrain(d1), correctStrain);
}

// BOOST_AUTO_TEST_CASE(InterpolationBStrain1D)
//{
//    MockInterpolation1D p(1);
//    BoostUnitTest::CheckEigenMatrix(p.GetBStrain(Eigen::VectorXd()), p.GetBGradient(Eigen::VectorXd()));
//}
//
//
// BOOST_AUTO_TEST_CASE(InterpolationBStrain3D)
//{
//    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(6, 9);
//    expected(0, 0) = dN0dX;
//    expected(1, 1) = dN0dY;
//    expected(2, 2) = dN0dZ;
//    expected(3, 1) = dN0dZ;
//    expected(3, 2) = dN0dY;
//    expected(4, 0) = dN0dZ;
//    expected(4, 2) = dN0dX;
//    expected(5, 0) = dN0dY;
//    expected(5, 1) = dN0dX;
//
//    expected(0, 3) = dN1dX;
//    expected(1, 7) = dN1dY;
//    expected(2, 5) = dN1dZ;
//    expected(3, 4) = dN1dZ;
//    expected(3, 5) = dN1dY;
//    expected(4, 3) = dN1dZ;
//    expected(4, 5) = dN1dX;
//    expected(5, 3) = dN1dY;
//    expected(5, 4) = dN1dX;
//
//    expected(0, 6) = dN2dX;
//    expected(1, 7) = dN2dY;
//    expected(2, 8) = dN2dZ;
//    expected(3, 7) = dN2dZ;
//    expected(3, 8) = dN2dY;
//    expected(4, 6) = dN2dZ;
//    expected(4, 8) = dN2dX;
//    expected(5, 6) = dN2dY;
//    expected(5, 7) = dN2dX;
//
//    MockInterpolation3D p(3);
//    BoostUnitTest::CheckEigenMatrix(p.GetBStrain(Eigen::VectorXd()), expected);
//}
//
