#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "mechanics/interpolation/InterpolationTrussLobatto.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    return {Eigen::VectorXd::Constant(1, 1.), Eigen::VectorXd::Constant(1, 0.), Eigen::VectorXd::Constant(1, 1.)};
}

BOOST_AUTO_TEST_CASE(InterpolationTrussLobatto)
{
    int order = 5;
    NuTo::InterpolationTrussLobatto interpolation(1, order);

    std::vector<Eigen::VectorXd> rPoints = GetTestPoints();

    BOOST_TEST_MESSAGE("Checking copy move...");
    NuTo::Test::CheckCopyMove<NuTo::InterpolationTrussLobatto>();
    BOOST_TEST_MESSAGE("Checking match of shape functions and node positions...");
    NuTo::Test::CheckShapeFunctionsAndNodePositions(interpolation);
    BOOST_TEST_MESSAGE("Checking patrition of unity...");
    NuTo::Test::CheckPartitionOfUnity(interpolation, rPoints);
    BOOST_TEST_MESSAGE("Checking shape function derivatives via CDF...");
    NuTo::Test::CheckDerivativeShapeFunctionsCDF(interpolation, rPoints);
}
