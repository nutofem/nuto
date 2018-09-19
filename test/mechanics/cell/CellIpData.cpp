#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "nuto/mechanics/cell/CellIpData.h"
#include <type_traits>

BOOST_AUTO_TEST_CASE(CellIpDataMemoizationB)
{
    Eigen::MatrixXd dNdXi = Eigen::MatrixXd::Random(3, 2);
    Eigen::VectorXd nodalValues(6);
    nodalValues << 0, 0, 1, 0, 0, 1;

    fakeit::Mock<NuTo::CoordinateElementInterface> mockElement;
    Method(mockElement, GetDerivativeShapeFunctions) = dNdXi;
    Method(mockElement, GetNMatrix) = Eigen::MatrixXd::Ones(2, 6);
    Method(mockElement, ExtractCoordinates) = nodalValues;

    fakeit::Mock<NuTo::DofElementInterface> mockDofElement;
    Method(mockDofElement, ExtractNodeValues) = nodalValues;
    Method(mockDofElement, GetDerivativeShapeFunctions) = dNdXi;
    Method(mockDofElement, GetNMatrix) = Eigen::MatrixXd::Ones(2, 6);

    fakeit::Mock<NuTo::ElementCollection> elements;
    Method(elements, DofElement) = mockDofElement.get();
    Method(elements, CoordinateElement) = mockElement.get();

    fakeit::Mock<NuTo::Nabla::Interface> mockGradientOperator;
    Method(mockGradientOperator, operator()) = Eigen::Matrix3d::Ones();

    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({1. / 3., 1. / 3.});

    NuTo::Jacobian jac(nodalValues, dNdXi);
    NuTo::CellData cellData(elements.get(), 0);
    NuTo::CellIpData ipData(cellData, jac, ipCoords, 0);

    NuTo::DofType d0("dof0", 2);
    for (int i = 0; i < 10; ++i)
    {
        auto B = ipData.B(d0, mockGradientOperator.get());
        BoostUnitTest::CheckEigenMatrix(B, Eigen::Matrix3d::Ones());
    }

    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(mockGradientOperator, operator())).Exactly(1));

    NuTo::ScalarDofType d1("Scalar");
    BOOST_CHECK((std::is_same<decltype(ipData.Value(d0)), Eigen::VectorXd>::value));
    BOOST_CHECK((std::is_same<decltype(ipData.Value(d1)), double>::value));
}
