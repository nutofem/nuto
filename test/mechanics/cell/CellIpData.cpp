#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellIpData.h"

void AcceptDouble(double, NuTo::DofType d)
{
    BOOST_CHECK(d.GetName() == "Scalar");
}
void AcceptVector(Eigen::VectorXd, NuTo::DofType d)
{
    BOOST_CHECK(d.GetName() != "Scalar");
}


BOOST_AUTO_TEST_CASE(CellIpDataMemoizationB)
{
    NuTo::DerivativeShapeFunctionsNatural dNdXi = Eigen::MatrixXd::Random(3, 2);
    NuTo::NodeValues nodalValues(6);
    nodalValues << 0, 0, 1, 0, 0, 1;

    fakeit::Mock<NuTo::ElementInterface> mockElement;
    Method(mockElement, GetDerivativeShapeFunctions) = dNdXi;
    Method(mockElement, GetNMatrix) = Eigen::MatrixXd::Ones(2, 6);
    Method(mockElement, ExtractNodeValues) = nodalValues;

    fakeit::Mock<NuTo::ElementCollection> elements;
    Method(elements, DofElement) = mockElement.get();
    Method(elements, CoordinateElement) = mockElement.get();

    fakeit::Mock<NuTo::Nabla::Interface> mockGradientOperator;
    Method(mockGradientOperator, operator()) = Eigen::Matrix3d::Ones();

    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({1. / 3., 1. / 3.});

    NuTo::Jacobian jac(nodalValues, dNdXi, 2);
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
    AcceptDouble(ipData.Value(d1), d1);
    AcceptVector(ipData.Value(d0), d0);
}
