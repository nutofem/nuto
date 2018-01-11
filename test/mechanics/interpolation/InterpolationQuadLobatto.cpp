#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "mechanics/interpolation/InterpolationQuadLobatto.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    constexpr double a = 0.577350269189626;
    return {Eigen::Vector2d({-a, -a}), Eigen::Vector2d({-a, a}), Eigen::Vector2d({a, a}), Eigen::Vector2d({a, -a})};
}

BOOST_AUTO_TEST_CASE(InterpolationQuadLobatto)
{
    int order = 5;
    NuTo::InterpolationQuadLobatto interpolation(order);

    std::vector<Eigen::VectorXd> rPoints = GetTestPoints();

    BOOST_TEST_MESSAGE("Checking copy move...");
    NuTo::Test::CheckCopyMove<NuTo::InterpolationQuadLobatto>();
    BOOST_TEST_MESSAGE("Checking match of shape functions and node positions...");
    NuTo::Test::CheckShapeFunctionsAndNodePositions(interpolation);
    BOOST_TEST_MESSAGE("Checking partition of unity...");
    NuTo::Test::CheckPartitionOfUnity(interpolation, rPoints);
    BOOST_TEST_MESSAGE("Checking shape function derivatives via CDF...");
    NuTo::Test::CheckDerivativeShapeFunctionsCDF(interpolation, rPoints);
}
