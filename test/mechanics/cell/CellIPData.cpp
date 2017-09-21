#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/interpolation/CellInterpolationIGA.h"
#include "mechanics/cell/CellIPData.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/CellInterpolationFEM.h"
#include "mechanics/IGA/NURBS.h"

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
    NuTo::CellInterpolationFEM e0({&n0, &n1, &n2}, interpolation0);
    NuTo::DofType d0("dof0", 2, 0);

    NuTo::NodeSimple n3 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(1));
    NuTo::NodeSimple n4 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(3));
    NuTo::NodeSimple n5 = NuTo::NodeSimple(Eigen::Matrix<double, 1, 1>::Constant(7));
    fakeit::Mock<NuTo::InterpolationSimple> interpolation1;
    Method(interpolation1, GetDerivativeShapeFunctions) = MockDerivatives2D();
    NuTo::CellInterpolationFEM e1({&n3, &n4, &n5}, interpolation1.get());
    NuTo::DofType d1("dof0", 1, 1);

    NuTo::DofContainer<NuTo::CellInterpolationBase*> elements;
    elements[d0] = &e0;
    elements[d1] = &e1;

    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({1. / 3., 1. / 3.});

    NuTo::Jacobian<2> jac(e0.ExtractNodeValues(), e0.GetDerivativeShapeFunctions(ipCoords));

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

BOOST_AUTO_TEST_CASE(CellIGA)
{
    std::vector<std::vector<NuTo::NodeSimple*>> controlPoints;

    std::vector<NuTo::NodeSimple*> row1;
    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({0, 0}));
    row1.push_back(&n0);
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({1, 0}));
    row1.push_back(&n1);
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({2, 0}));
    row1.push_back(&n2);
    controlPoints.push_back(row1);

    std::vector<NuTo::NodeSimple*> row2;
    NuTo::NodeSimple n3 = NuTo::NodeSimple(Eigen::Vector2d({0, 1}));
    row2.push_back(&n3);
    NuTo::NodeSimple n4 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    row2.push_back(&n4);
    NuTo::NodeSimple n5 = NuTo::NodeSimple(Eigen::Vector2d({2, 1}));
    row2.push_back(&n5);
    controlPoints.push_back(row2);

    std::vector<NuTo::NodeSimple*> row3;
    NuTo::NodeSimple n6 = NuTo::NodeSimple(Eigen::Vector2d({0, 2}));
    row3.push_back(&n6);
    NuTo::NodeSimple n7 = NuTo::NodeSimple(Eigen::Vector2d({1, 2}));
    row3.push_back(&n7);
    NuTo::NodeSimple n8 = NuTo::NodeSimple(Eigen::Vector2d({2, 2}));
    row3.push_back(&n8);
    controlPoints.push_back(row3);

    std::vector<double> knots1D = {0, 0, 0, 1, 1, 1};
    std::array<std::vector<double>, 2> knots = {knots1D, knots1D};

    std::vector<double> weights1D = {1, 1, 1};
    std::vector<std::vector<double>> weights = {weights1D, weights1D, weights1D};

    std::array<int, 2> degree = {2, 2};

    NuTo::NURBS<2> surface(knots, controlPoints, weights, degree);

    std::array<int, 2> knotIDsCell = {2, 2};
    NuTo::CellInterpolationIGA<2> iga(knotIDsCell, surface);

    NuTo::DofType d0("dof0", 2, 0);

    NuTo::DofContainer<NuTo::CellInterpolationBase*> elements;
    elements[d0] = &iga;

    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({-1. / 3., 1. / 3.});

    NuTo::Jacobian<2> jac(iga.ExtractNodeValues(), iga.GetDerivativeShapeFunctions(ipCoords));

    BoostUnitTest::CheckEigenMatrix(jac.Inv(), 0.5 * Eigen::Matrix2d::Identity());
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
