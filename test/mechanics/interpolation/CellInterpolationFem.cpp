#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include <type_traits>

#include "mechanics/interpolation/CellInterpolationFem.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"


BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::CellInterpolationFem>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::CellInterpolationFem>::value);
}


struct TestElement : public NuTo::CellInterpolationFem
{
    TestElement()
        : NuTo::CellInterpolationFem({&n0, &n1, &n2}, interpolation)
    {
    }

    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));
    NuTo::InterpolationTriangleLinear interpolation = NuTo::InterpolationTriangleLinear(2);
};

BOOST_AUTO_TEST_CASE(ExtractNodeValues)
{
    NuTo::NodeValues nodeValues = TestElement().ExtractNodeValues();
    BoostUnitTest::CheckVector(nodeValues, std::vector<double>{1, 1, 5, 1, 1, 7}, 6);
}

BOOST_AUTO_TEST_CASE(Interpolation)
{
    auto nodeValues = TestElement().ExtractNodeValues();
    Eigen::VectorXd ipCoord(2);
    ipCoord << 0.5, 0.5;
    //    auto N = TestElement().GetNMatrix(Eigen::Vector2d(0.5, 0.5));
    auto N = TestElement().GetNMatrix(ipCoord);
    BoostUnitTest::CheckVector(N * nodeValues, std::vector<double>{3, 4}, 2);
}

BOOST_AUTO_TEST_CASE(CacheN)
{
    fakeit::Mock<NuTo::InterpolationSimple> interpolation;
    Method(interpolation, GetShapeFunctions) = Eigen::Vector2d({42, 6174});
    Method(interpolation, GetDofDimension) = 1;
    Method(interpolation, GetNumNodes) = 2;

    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::CellInterpolationFem interpolationFem({&n0}, interpolation.get());

    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        interpolationFem.GetNMatrix(Eigen::Vector2d(0, 0));
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(interpolation, GetShapeFunctions)).Exactly(1));
}
