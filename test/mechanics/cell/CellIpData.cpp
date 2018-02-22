#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellIpData.h"

BOOST_AUTO_TEST_CASE(CellIpDataMemoizationB)
{
    NuTo::DerivativeShapeFunctionsNatural dNdXi = Eigen::MatrixXd::Random(3, 2);

    constexpr int ipId = 27;

    fakeit::Mock<NuTo::ElementInterface> mockElement;
    Method(mockElement, GetDerivativeShapeFunctions) = dNdXi;

    fakeit::Mock<NuTo::ElementCollection> elements;
    Method(elements, DofElement) = mockElement.get();

    fakeit::Mock<NuTo::B::Interface> mockGradientOperator;
    Method(mockGradientOperator, operator()) = Eigen::Matrix3d::Ones();

    NuTo::NaturalCoords ipCoords = Eigen::Vector2d({1. / 3., 1. / 3.});

    NuTo::NodeValues nodalValues(6);
    nodalValues << 0, 0, 1, 0, 0, 1;


    NuTo::Jacobian jac(nodalValues, dNdXi, 2);
    NuTo::CellIpData ipData(elements.get(), jac, ipCoords, ipId);
    BOOST_CHECK(ipData.GetIpId() == ipId);

    NuTo::DofType d0("dof0", 2);
    for (int i = 0; i < 10; ++i)
    {
        auto B = ipData.B(d0, mockGradientOperator.get());
        BoostUnitTest::CheckEigenMatrix(B, Eigen::Matrix3d::Ones());
    }

    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(mockGradientOperator, operator())).Exactly(1));
}
