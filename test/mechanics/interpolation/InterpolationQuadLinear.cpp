#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "TypeTraits.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"


BOOST_AUTO_TEST_CASE(InterpolationQuadLinearCopyMove)
{
    NuTo::Test::Copy<NuTo::InterpolationQuadLinear>();
    NuTo::Test::Move<NuTo::InterpolationQuadLinear>();
}

std::vector<Eigen::VectorXd> GetTestPoints()
{
    constexpr double a = 0.577350269189626;
    return {Eigen::Vector2d({-a, -a}), Eigen::Vector2d({-a, a}),
            Eigen::Vector2d({a, a}), Eigen::Vector2d({a, -a})};
}


BOOST_AUTO_TEST_CASE(InterpolationQuadLinearN)
{
    NuTo::InterpolationQuadLinear interpolation(1);
    NuTo::Test::CheckShapeFunctionsAndNodePositions(interpolation);
}


BOOST_AUTO_TEST_CASE(InterpolationQuadLinearUnity)
{
    NuTo::InterpolationQuadLinear interpolation(1);
    NuTo::Test::CheckPartitionOfUnity(interpolation, GetTestPoints());
}

BOOST_AUTO_TEST_CASE(InterpolationQuadLinearB)
{
    NuTo::InterpolationQuadLinear interpolation(1);
    NuTo::Test::CheckDerivativeShapeFunctionsCDF(interpolation, GetTestPoints());
}
