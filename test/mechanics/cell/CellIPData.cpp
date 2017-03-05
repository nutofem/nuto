#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/interpolation/InterpolationTriangle.h"
#include "mechanics/cell/CellIPData.h"
#include <iostream>

NuTo::DerivativeShapeFunctionsNatural MockDerivatives1D()
{
    NuTo::DerivativeShapeFunctionsNatural m(3, 1);
    m(0, 0) = 0;
    m(1, 0) = 10;
    m(2, 0) = 20;
    return m;
}

NuTo::DerivativeShapeFunctionsNatural MockDerivatives2D()
{
    NuTo::DerivativeShapeFunctionsNatural m(3, 2);
    m(0, 0) = 0;
    m(0, 1) = 1;
    m(1, 0) = 10;
    m(1, 1) = 11;
    m(2, 0) = 20;
    m(2, 1) = 21;
    return m;
}

NuTo::DerivativeShapeFunctionsNatural MockDerivatives3D()
{
    NuTo::DerivativeShapeFunctionsNatural m(3, 2);
    m(0, 0) = 0;
    m(0, 1) = 1;
    m(0, 2) = 2;
    m(1, 0) = 10;
    m(1, 1) = 11;
    m(1, 2) = 12;
    m(2, 0) = 20;
    m(2, 1) = 21;
    m(2, 2) = 22;
    return m;
}

BOOST_AUTO_TEST_CASE(CellIPData2D)
{
    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({1, 0}));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({0, 1}));
    NuTo::InterpolationTriangle interpolation0(NuTo::eInterpolation::GAUSS, 1, 2);
    NuTo::Element e0({&n0, &n1, &n2}, interpolation0);
    NuTo::DofType d0("dof0", 2, 0);

    NuTo::NodeSimple n3 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(1));
    NuTo::NodeSimple n4 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(3));
    NuTo::NodeSimple n5 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(7));
    fakeit::Mock<NuTo::Interpolation> interpolation1;
    Method(interpolation1, GetDerivativeShapeFunctions) = MockDerivatives2D();
    NuTo::Element e1({&n3, &n4, &n5}, interpolation1.get());
    NuTo::DofType d1("dof0", 1, 1);

    NuTo::DofContainer<const NuTo::Element*> elements;
    elements[d0] = &e0;
    elements[d1] = &e1;

    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({1. / 3., 1. / 3.});

    NuTo::Jacobian jac(e0.ExtractNodeValues(), e0.GetInterpolation().GetDerivativeShapeFunctions(ipCoords));

    NuTo::CellIPData ipData(elements, jac, ipCoords);

    BoostUnitTest::CheckEigenMatrix(jac.Inv(), Eigen::Matrix2d::Identity());

    Eigen::MatrixXd correctGradient(2, 3);
    correctGradient(0, 0) = 0;
    correctGradient(0, 1) = 10;
    correctGradient(0, 2) = 20;
    correctGradient(1, 0) = 1;
    correctGradient(1, 1) = 11;
    correctGradient(1, 2) = 21;
    BoostUnitTest::CheckEigenMatrix(ipData.GetBMatrixGradient(d1), correctGradient);


    Eigen::MatrixXd correctStrain = Eigen::MatrixXd::Zero(3, 6);
    correctStrain(0, 0) = 0;
    correctStrain(1, 1) = 1;
    correctStrain(2, 0) = 1;
    correctStrain(2, 1) = 0;

    correctStrain(0, 2) = 10;
    correctStrain(1, 3) = 11;
    correctStrain(2, 2) = 11;
    correctStrain(2, 3) = 10;

    correctStrain(0, 4) = 20;
    correctStrain(1, 5) = 21;
    correctStrain(2, 4) = 21;
    correctStrain(2, 5) = 20;
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
//    expected(0, 0) = 0;
//    expected(1, 1) = 1;
//    expected(2, 2) = 2;
//    expected(3, 1) = 2;
//    expected(3, 2) = 1;
//    expected(4, 0) = 2;
//    expected(4, 2) = 0;
//    expected(5, 0) = 1;
//    expected(5, 1) = 0;
//
//    expected(0, 3) = 10;
//    expected(1, 4) = 11;
//    expected(2, 5) = 12;
//    expected(3, 4) = 12;
//    expected(3, 5) = 11;
//    expected(4, 3) = 12;
//    expected(4, 5) = 10;
//    expected(5, 3) = 11;
//    expected(5, 4) = 10;
//
//    expected(0, 6) = 20;
//    expected(1, 7) = 21;
//    expected(2, 8) = 22;
//    expected(3, 7) = 22;
//    expected(3, 8) = 21;
//    expected(4, 6) = 22;
//    expected(4, 8) = 20;
//    expected(5, 6) = 21;
//    expected(5, 7) = 20;
//
//    MockInterpolation3D p(3);
//    BoostUnitTest::CheckEigenMatrix(p.GetBStrain(Eigen::VectorXd()), expected);
//}
//
