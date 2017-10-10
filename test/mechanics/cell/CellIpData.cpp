#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/cell/CellIpData.h"
#include "mechanics/nodes/NodeSimple.h"

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
    NuTo::DerivativeShapeFunctionsNatural m(3, 3);
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
    fakeit::Mock<NuTo::ElementInterface> mockElement;
    Method(mockElement, GetDerivativeShapeFunctions) = MockDerivatives2D();

    fakeit::Mock<NuTo::ElementCollection> elements;
    Method(elements, DofElement) = mockElement.get(); 


    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({1. / 3., 1. / 3.});

    NuTo::NodeValues nodalValues(6);
    nodalValues << 0, 0, 1, 0, 0, 1;

    NuTo::DerivativeShapeFunctionsNatural derivativeForJacobian =
            NuTo::ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder1(ipCoords);

    NuTo::Jacobian jac(nodalValues, derivativeForJacobian);
    NuTo::CellIpData ipData(elements.get(), jac, ipCoords);

    BoostUnitTest::CheckEigenMatrix(jac.Inv(), Eigen::Matrix2d::Identity());

    Eigen::MatrixXd correctGradient(2, 3);
    correctGradient(0, 0) = dN0dX;
    correctGradient(0, 1) = dN1dX;
    correctGradient(0, 2) = dN2dX;
    correctGradient(1, 0) = dN0dY;
    correctGradient(1, 1) = dN1dY;
    correctGradient(1, 2) = dN2dY;
    
    NuTo::DofType d0("dof0", 2);
    BoostUnitTest::CheckEigenMatrix(ipData.GetBMatrixGradient(d0), correctGradient);


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

    BoostUnitTest::CheckEigenMatrix(ipData.GetBMatrixStrain(d0), correctStrain);
}

BOOST_AUTO_TEST_CASE(InterpolationBStrain3D)
{
    fakeit::Mock<NuTo::ElementInterface> mockElement;
    Method(mockElement, GetDerivativeShapeFunctions) = MockDerivatives3D();

    fakeit::Mock<NuTo::ElementCollection> elements;
    Method(elements, DofElement) = mockElement.get(); 

    NuTo::NaturalCoords ipCoords = Eigen::Vector3d({1. / 3., 1. / 3., 1. / 3.});

    NuTo::NodeValues nodalValues(12);
    nodalValues << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;

    NuTo::DerivativeShapeFunctionsNatural derivativeForJacobian =
            NuTo::ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder1();

    NuTo::Jacobian jac(nodalValues, derivativeForJacobian);
    NuTo::CellIpData ipData(elements.get(), jac, ipCoords);

    BoostUnitTest::CheckEigenMatrix(jac.Inv(), Eigen::Matrix3d::Identity());
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(6, 9);
    expected(0, 0) = dN0dX;
    expected(1, 1) = dN0dY;
    expected(2, 2) = dN0dZ;
    expected(3, 1) = dN0dZ;
    expected(3, 2) = dN0dY;
    expected(4, 0) = dN0dZ;
    expected(4, 2) = dN0dX;
    expected(5, 0) = dN0dY;
    expected(5, 1) = dN0dX;

    expected(0, 3) = dN1dX;
    expected(1, 4) = dN1dY;
    expected(2, 5) = dN1dZ;
    expected(3, 4) = dN1dZ;
    expected(3, 5) = dN1dY;
    expected(4, 3) = dN1dZ;
    expected(4, 5) = dN1dX;
    expected(5, 3) = dN1dY;
    expected(5, 4) = dN1dX;

    expected(0, 6) = dN2dX;
    expected(1, 7) = dN2dY;
    expected(2, 8) = dN2dZ;
    expected(3, 7) = dN2dZ;
    expected(3, 8) = dN2dY;
    expected(4, 6) = dN2dZ;
    expected(4, 8) = dN2dX;
    expected(5, 6) = dN2dY;
    expected(5, 7) = dN2dX;

    NuTo::DofType d0("dof0", 3);
    BoostUnitTest::CheckEigenMatrix(ipData.GetBMatrixStrain(d0), expected);
}
